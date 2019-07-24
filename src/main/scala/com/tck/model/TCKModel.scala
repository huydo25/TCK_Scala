package com.tck.model

import breeze.linalg.{ Axis, DenseMatrix, DenseVector, det, inv, sum}
import breeze.numerics.{exp, pow, sqrt}

import scala.math._
import org.apache.spark.rdd.RDD
import com.tck.utils.Utils._

import scala.util.Random

class TCKModel(minNi: Double, minVi: Int, maxVi: Int,
                 minTi: Int , maxTi: Int, ii: Int, gi: Int) extends Serializable {
  /** INPUTS
    * X: data array of size N x T x V
    * C: number of mixture components (optional)
    * minN: min percentage of sub-sample (optional)
    * minV: min number of dimensions (optional)
    * maxV: max number of dimensions (optional)
    * minT: min length of time segments (optional)
    * maxT: max length of time segments (optional)
    * I: number of iterations (optional)
    * missing: binary indicator. 1 if there is missing data and 0 if not
    */

  var minN: Double = minNi
  var minV: Int = minVi
  var maxV: Int = maxVi
  var minT: Int = minTi
  var maxT: Int = maxTi
  var I: Int = ii
  var G: Int = gi

  def map_em(x: Array[Array[Array[Double]]], C:Int,  missing: Int): Cell= {
    /** MAP_EM - fit a GMM to time series data with missing values using MAP-EM
      * OUTPUTS
      * Q: cluster posterior probabilities
      * mu: cluster means (time dependant + variable dependant)
      * s2: cluster variances (variable dependant)
      * theta: cluster priors
      * dim_idx: indexes of the subset of dimension considered
      * time_idx: indexes of the subset of time intervals considered
      */
    val N = x.length
    val T = x(0).length
    val V = x(0)(0).length

    // Hyper parameters for mean prior (a0, b0) and the std dev prior (n0) of the mixture components
    val a0: Double = (1.0 - 0.001) * random + 0.001
    val b0: Double = (0.2 - 0.005) * random + 0.005
    val n0: Double = (0.2 - 0.001) * random + 0.001
    /*val a0 = 0.8409
    val b0 = 0.0546 // used for debugs
    val n0 = 0.1630*/

    //Randomly subsample dimensions, time intervals and samples
    var sN: Int = 0
    if (N > 100) {
      sN = math.round(minN*N).toInt + Random.nextInt(N - math.round(minN*N).toInt +1)
      //sN = 180
    } else {
      sN = math.round(0.9 * N).toInt
    }
    var sub_idx: Array[Int] = Array.ofDim(sN)
    sub_idx = Random.shuffle(0 to N-1).take(sN).sortWith(_<_).toArray
    //sub_idx = (0 until sN).toArray // for debugs

    val sV = minV + Random.nextInt(maxV - minV + 1)
    val dim_idx = Random.shuffle(0 to V - 1).take(sV).sortWith(_ < _) //generate sV (sorted) integers between 1 and V

    val t1 = 1 + Random.nextInt(T - minT + 1)
    val t2 = (t1 + minT - 1) + Random.nextInt(min(T, t1 + maxT - 1) - (t1 + minT - 1) + 1) - 1
    val sT = t2 - t1 +1
    val time_idx = t1 to t2 // generate sT continuous integers from t1 to t2

//    val sT = 15
//    val time_idx = (0 until sT)

    val sX: Array[Array[Array[Double]]] = Array.ofDim(sub_idx.length, time_idx.length, dim_idx.length)
    sX.indices.map( n => sX(n).indices.map(t=> sX(n)(t).indices.map( v => sX(n)(t)(v) = x(sub_idx(n))(time_idx(t))(dim_idx(v)))))


    // initialize model parameters
    var theta: Array[Double] = Array.fill(C)(1 / C) // cluster priors        (1 x C)
    var mu: Array[Array[Array[Double]]] = Array.fill(sT, sV, C)(0.0) // cluster means         (sT x sV x C)
    var s2: Array[Array[Double]] = Array.fill(sV, C)(0.0) // cluster variances     (sV x C)
    var Q: Array[Array[Double]] = Array.fill(sN, C)(0.0)

    if (missing == 0) {
      //Calculate empirical moments
      // init variable

      val rsX: Array[Array[Array[Double]]] = Array.ofDim(sV, sT, sN)

      rsX.indices.map(x => {
        rsX(x).indices.map(y => rsX(x)(y) = sX.map(_ (y)).map(_ (x)))
      })

      val mu_0: Array[Array[Double]] = Array.fill(sT, sV)(0) // prior mean over time and variables (sT x sV)
      mu_0.indices.map(x => mu_0(x) = nanMean(rsX.map(_ (x)), 1))

      val rsX1: Array[Array[Double]] = Array.ofDim(sV, sT * sN)
      rsX1.indices.map(x => {
        var t: Array[Double] = Array()
        rsX(x).indices.map(y => t = t ++ sX.map(_ (y)).map(_ (x)))
        rsX1(x) = t
      })

      val s_0: Array[Double] = nanStd(rsX1.transpose)
      val s2_0 = s_0.map(x => x * x)

      /*val newX = sX.context.parallelize(sXa)
      newX.foreach(x => println(x.deep.mkString(" ")))*/

      val invS_0: Array[DenseMatrix[Double]] = Array.ofDim(2)
      val T1: DenseMatrix[Double] = DenseMatrix(Array.fill(sT)(1 to sT toArray).map(x => x.map(y => y.toDouble)).transpose: _*)
      val T2: DenseMatrix[Double] = DenseMatrix(Array.fill(sT)(1 to sT toArray).map(x => x.map(y => y.toDouble)): _*)

      invS_0.zipWithIndex.map(v => {
        var t = s_0(v._2) * b0 * exp(-a0 * pow(T1 - T2, 2))
        if (det(t) == 0) { // check if the matrix can be inverted
          t = t + 0.1 * t(1, 1) * DenseMatrix.eye[Double](sT) // add a small number to the diagonal if S_0 is not invertible
        }
        invS_0(v._2) = inv(t)
      })

      // initialize model parameters
      theta = Array.fill(C)(1.0 / C) // cluster priors        (1 x C)
      mu = Array.fill(sT, sV, C)(0.0) // cluster means         (sT x sV x C)
      s2 = Array.fill(sV, C)(0.0) // cluster variances     (sV x C)
      Q = Array.fill(sN, C)(0.0) // cluster assignments   (sN x C)x: RDD[DenseMatrix[Double]]

      for (i <- 0 until I) {
      //Vector.range(0, G*(C-1)).par.foreach( i => {
        if (i == 0) {
          val cluster: Array[Int] = Array.fill(sN)(Random.nextInt(C))
          //val cluster = Array(9, 18, 0, 14, 16, 29, 36, 18, 0, 36, 8, 38, 21, 20, 1, 10, 10, 25, 21, 23, 15, 35, 39, 37, 39, 22, 24, 6, 25, 8, 19, 39, 1, 36, 18, 7, 19, 11, 2, 0, 19, 10, 0, 35, 27, 30, 22, 0, 10, 12, 22, 3, 37, 5, 18, 9, 21, 7, 38, 10, 30, 23, 3, 25, 31, 1, 16, 6, 18, 18, 39, 12, 9, 5, 25, 38, 11, 26, 39, 37, 26, 29, 22, 29, 18, 26, 14, 5, 25, 6, 15, 13, 21, 0, 33, 18, 14, 39, 6, 37, 10, 21, 6, 18, 14, 9, 23, 27, 24, 26, 7, 12, 6, 39, 8, 18, 18, 36, 38, 29, 37, 39, 16, 37, 12, 1, 32, 32, 28, 32, 9, 17, 6, 32, 20, 5, 6, 33, 23, 6, 30, 12, 19, 26, 26, 16, 0, 30, 0, 17, 23, 28, 30, 32, 11, 36, 31, 13, 29, 39, 35, 17, 3, 6, 34, 16, 5, 25, 21, 11, 19, 36, 18, 28, 6, 3, 34, 26, 26, 6)
          val temp: Array[Int] = (1 to C).toArray
          Q.indices.map(a => Q(a).indices.map(j => if (cluster(a) == temp(j)) Q(a)(j) = 1 else Q(a)(j) = 0))
        } else {
          // start
          val distr_c: Array[Array[Array[Double]]] = Array.ofDim(sN, sT, sV)
          for (c <- 0 until C) {
            val temp1: Array[Array[Double]] = Array.ofDim(sN, sV * sT)
            distr_c.indices.map(n =>
              distr_c(n).indices.map(t =>
                distr_c(n)(t).indices.map(v => {
                  distr_c(n)(t)(v) = normpdf(sX(n)(t)(v), mu(t)(v)(c), sqrt(s2(v)(c)))
                  if (distr_c(n)(t)(v) < normpdf(3))
                    distr_c(n)(t)(v) = normpdf(3)
                  temp1(n)(v * sT + t) = distr_c(n)(t)(v)
                })))

            val prod_distrc_c: Array[Double] = Array.fill(sN)(1)
            temp1.indices.map(n => temp1(n).indices.map(vt => prod_distrc_c(n) *= temp1(n)(vt)))
            prod_distrc_c.indices.map(n => Q(n)(c) = prod_distrc_c(n) * theta(c))
          }
          val sumtQ = Array.fill(C)(Q.map(_.sum)).transpose
          Q.indices.map(i => Q(i).indices.map(j => Q(i)(j) /= sumtQ(i)(j)))
        }
        // update mu, s2 and theta
        for (c <- 0 until C) {
          val sumQ = Q.map(_(c)).sum

          theta(c) = sumQ / sN
          //println(theta.deep.mkString(" "))
          for (v <- 0 until sV) {
            val var2 = sT * sumQ
            val temp_var2 = (DenseMatrix(sX.map(_.map(_ (v))): _*) - DenseMatrix(Array.fill(sN)(mu.map(_.map(_ (c))).map(_ (v))): _*)).map(x => x * x)
            val var1 = DenseVector(Q.map(_ (c)): _*).t dot sum(temp_var2, Axis._1).t
            s2(v)(c) = (n0 * s2_0(v) + var1) / (n0 + var2)

            val A = invS_0(v) + (sumQ / s2(v)(c)) * DenseMatrix.eye[Double](sT)
            val b = invS_0(v) * DenseMatrix(mu_0.map(_ (v)): _*) + DenseMatrix(sX.map(_.map(_ (v))): _*).t * DenseMatrix(Q.map(_ (c)): _*) / s2(v)(c)
            val temp_r: Array[Double] = (A \ b).toArray
            for (t <- 0 until sT) {
              mu(t)(v)(c) = temp_r(t)
            }
          }
        }
      }
      Q = posterior(x, C, mu, s2, theta, dim_idx, time_idx, missing)
      //println(Q.deep.mkString("\n"))
    }
    new Cell(Q , mu, s2, theta, dim_idx, time_idx)
  }

  def posterior(x: Array[Array[Array[Double]]], C:Int, mu: Array[Array[Array[Double]]], s2: Array[Array[Double]],
                theta: Array[Double], dim_idx: Seq[Int], time_idx: Seq[Int], missing: Int): Array[Array[Double]] = {

    /* GMMposterior - Evaluate the posterior for the data X of the GMM described by C, mu, s2 and theta
     * INPUTS
     *  X: data array of size N x T x V
     *  C: number of mixture components (optional)
     *  mu: cluster means over time and variables (V x T)
     *  s2: cluster stds over variables (sV x 1)
     *  theta: cluster priors
     *  dim_idx: subset of variables to be used in the clustering
     *  time_idx: subset of time intervals to be used in the clustering
     *  missing: binary indicator. 1 if there is missing data and 0 if not
     * OUTPUTS
     *  Q: posterior
     **/
    val N = x.length
    val Q: Array[Array[Double]] = Array.fill(N, C)(0.0)
    val sV = dim_idx.length
    val sT = time_idx.length

    val sX: Array[Array[Array[Double]]] = Array.ofDim(N, time_idx.length, dim_idx.length)

    sX.indices.map( n => sX(n).indices.map(t=> sX(n)(t).indices.map( v => sX(n)(t)(v) = x(n)(time_idx(t))(dim_idx(v)))))

    if (missing == 0) {

      // Compute GMM posterior
      val distr_c: Array[Array[Array[Double]]] = Array.ofDim(N, sT, sV)
      for (c <- 0 until C) {
        var temp1: Array[Array[Double]] = Array.ofDim(N, sV * sT)
        distr_c.indices.map(n =>
          distr_c(n).indices.map(t =>
            distr_c(n)(t).indices.map(v => {
              distr_c(n)(t)(v) = normpdf(sX(n)(t)(v), mu(t)(v)(c), sqrt(s2(v)(c)))
              if (distr_c(n)(t)(v) < normpdf(3))
                distr_c(n)(t)(v) = normpdf(3)
              temp1(n)(v * sT + t) = distr_c(n)(t)(v)
            })))
        var prod_distrc_c: Array[Double] = Array.fill(N)(1)
        temp1.indices.map(n => temp1(n).indices.map(vt => prod_distrc_c(n) *= temp1(n)(vt)))
        prod_distrc_c.indices.map(n => Q(n)(c) = prod_distrc_c(n) * theta(c))
      }
      val sumtQ = Array.fill(C)(Q.map(_.sum)).transpose
      Q.indices.map(i => Q(i).indices.map(j => Q(i)(j) /= sumtQ(i)(j)))
    }
    Q
  }

  def train(x: Array[Array[Array[Double]]] ) : (Array[Cell],Int,Int) = {
    val N = x.length
    val T = x(0).length
    val V = x(0)(0).length
    var C = 10
    if(N < 100) {
      C = 10
    } else {
      C = 40
    }

    if (V == 1){
      val minV = 1
    } else{
      val minV = 2
    }
    var check:Double = 0.0
    var missing = 0
    x.map( a => a.map(b => b.map(c => if(c.isNaN) check += 1.0) ))
    if (check == 0.0){
      missing = 0
      println("The dataset does not contain missing data\n\n")
    } else {
      missing = 1
      println("The dataset contains missing data\n\n")
    }

    println(s"Training the TCK using the following parameters:\n C = $C, G =$G\n " +
      s"Number of MTS for each GMM: ${floor(minN*N)} - $N (${floor(minN*100)} - 100 percent)\n " +
      s"Number of attributes sampled from [$minV, $maxV]\n Length of time segments sampled from [$minT, $maxT]\n\n ")

    val res : Array[Cell] = Array.ofDim((G*(C-1)))
    //Vector.range(0, G*(C-1)).par.foreach( i => {
      //Thread.sleep(100)
    for (i <- 0 until (G*(C-1))) {
      println(i)
      val c: Int = (floor((i)/G) + 2).toInt
      res(i) = map_em(x, c, missing)
    }//)
    (res, C, G)
  }

  def computeKM(x: Array[Cell], C: Int, G: Int, marking: Int = 0, xte: Array[Array[Array[Double]]]): Array[Array[Double]] ={
     /*function [ K ] = TCK(GMM, C, G, Xte)
      *TCK -  compute TCK kernel matrix between training data and test data Xte
      *INPUTS
      *  GMM : Cell output from the function trainTCK
      *  C: Second output from trainTCK
      *  G: Third output from trainTCK
      *  marking: mark to compute in-sample kernel matrix
      *  Xte: data array of size Nte x T x V, where Nte is the number of
      *  multivariate time series, T the length and V the number of attributes.
      *OUTPUTS
      *  K: kernel matrix
      **/
    if (marking == 0 ){
      // Check if the dataset contains mising elements
      var check:Double = 0.0
      var missing = 0
      xte.map( a => a.map(b => b.map(c => if(c.isNaN) check += 1.0) ))

      if (check == 0.0){
        missing = 0
      } else {
        missing = 1
      }
      val K : Array[Array[Double]] = Array.fill(x(0).Q.length,xte.length)(0.0)
      var result: DenseMatrix[Double] = DenseMatrix.zeros(x(0).Q.length,xte.length)
      //Vector.range(0, G*(C-1)).par.foreach( i => {
      for (i <- 0 until G*(C-1)){
        //println(i)
        val c = floor(i/G).toInt + 2
        val q : Array[Array[Double]] = posterior( xte, c, x(i).mean, x(i).std, x(i).Th, x(i).D, x(i).T, missing)
        result = result + DenseMatrix(x(i).Q:_*) * DenseMatrix(q:_*).t
      }//)
      result = result.t
      for (k <- 0 until K.length) {
        K(k) = result(::,k).toArray
      }
      K
    }  else { // in-sample kernel matrix
      val K : Array[Array[Double]] = Array.fill(x(0).Q.length,x(0).Q.length)(0)
      var result: DenseMatrix[Double] = DenseMatrix.zeros(x(0).Q.length,x(0).Q.length)

      for (i <- 0 until G*(C-1)){
        result = result + DenseMatrix(x(i).Q:_*) * DenseMatrix(x(i).Q:_*).t
      }
      for (k <- 0 until result.rows ){
        K(k) = result(::,k).toArray
      }
      K
    }
  }


}
