package com.tck.models

import com.tck.utils.Utils._
import scala.math._
import com.tck.models.GMM_MAP_EM._

object trainTCK{
  def trainTCK(x: Array[Array[Array[Double]]],
               C: Int = 10, minN: Double = 0.8,
               minV: Int = 2, maxV: Int = 2,
               minT : Int =6, maxT: Int = 25,
               I: Int = 20,  G: Int = 5): (List[ (Array[Array[Double]], Array[Array[Array[Double]]], Array[Array[Double]],Array[Double], Array[Int], Array[Int])], Int, Int) = {

    // trainTCK - Train the TCK
    //
    //INPUTS
    //   x: data array of size N x T x V, where N is the number of multivariate time series,
    //      T the length and V the number of attributes.
    //   minN: min percentage of sub-sample (optional)
    //   minV: min number of attributes for each GMM (optional)
    //   maxV: max number of attributes for each GMM  (optional)
    //   minT: min length of time segments for each GMM (optional)
    //   maxT: max length of time segments for each GMM (optional)
    //   C: max number of mixture components for each GMM (optional)
    //   G: number of randomization for each number of components (optional)
    //   I: number of iterations (optional)
    //OUTPUTS
    //   res: A cell of size ((C-1)*G,6) that for each q = 1:(C-1)*G contain
    //   Q: cluster posterior probabilities
    //   mu: cluster means (time dependant + variable dependant)
    //   s2: cluster variances (variable dependant)
    //   theta: cluster priors
    //   dim_idx: indexes of the subset of dimension considered
    //   time_idx: indexes of the subset of time intervals considered
    //  C
    //  G

    val N = x.length
    val T = x(0).length
    val V = x(0)(0).length
    // optional parameters
    assert(C >= 2, "C must be larger than 1" )
    if(N < 100) {
      val C = 10
    } else {
      val C = 40
    }

    if (V == 1){
      val minV = 1
    } else{
      val minV = 2
    }
//    println(maxV)
//    println(V)
    assert(minN > 0 && minN <= 1, "The minimum percentage of subsample must be in (0,1]" )
    assert(minV >= 1 && minV <= V, "The minimum number of variables must be in [1,V]")
    assert(maxV >= 1 && maxV <= V, "The maximum number of variables must be in [1,V]")
    assert(minT >= 1 && minT <= T, "The minimum number of variables must be in [1,T]")
    assert(maxT >= 1 && maxT <= T, "The maximum number of variables must be in [1,T]")


    // Check if there is missing data in the dataset.
    val nan_idx = isNaN(x)
    var sum = 0
    for (i <- 0 until nan_idx(0).length){
      for (j <- 0 until nan_idx(0)(0).length ){
        sum += nan_idx.map(_(i)).map(_(j)).sum.toInt
      }
    }
    var missing = 0
    if (sum > 0 ){
      missing = 1
      println("The dataset contains missing data\n\n")
    } else {
      missing = 0
      println("The dataset does not contain missing data\n\n")
    }
    println(s"Training the TCK using the following parameters:\n C = $C, G =$G\n " +
            s"Number of MTS for each GMM: ${floor(minN*N)} - $N (${floor(minN*100)} - 100 percent)\n " +
            s"Number of attributes sampled from [$minV, $maxV]\n Length of time segments sampled from [$minT, $maxT]\n\n ")

    var res : List[ (Array[Array[Double]], Array[Array[Array[Double]]], Array[Array[Double]],Array[Double], Array[Int], Array[Int])]  = List()

    for (i <- 0 until (G*(C-1))){
      val c: Int = (floor((i-1)/G) + 2).toInt
      res = res :+ GMM_MAP_EM.GMM_MAP_EM(x,c,minN,minV,V,minT,maxT,I,missing)
    }
    //println(GMM_MAP_EM.GMM_MAP_EM(x,(floor((5-1)/G) + 2).toInt,minN,minV,V,minT,maxT,I,missing)._1.deep.mkString("\n"))
    (res, C, G)
  }
}
