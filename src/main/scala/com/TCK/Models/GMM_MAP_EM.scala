package com.TCK.Model

import breeze.linalg._

import scala.math._
import breeze.numerics._
import com.TCK.Ultils.Ultils.isNaN

import scala.util.Random

object  GMM_MAP_EM{
  def GMM_MAP_EM(x: Array[Array[Array[Double]]],
                  C: Int = 40, minN: Double = 0.8,
                  minV: Int = 2, maxV: Int = 100,
                  minT : Int =6, maxT: Int = 25,
                  I: Int = 20,  missing: Int = 2): Unit ={
//    MAP_EM - fit a GMM to time series data with missing values using MAP-EM
//
//    INPUTS
//      X: data array of size N x T x V
//      C: number of mixture components (optional)
//      minN: min percentage of subsample (optional)
//      minV: min number of dimensions (optional)
//      maxV: max number of dimensions (optional)
//      minT: min length of time segments (optional)
//      maxT: max length of time segments (optional)
//      I: number of iterations (optional)
//      missing: binary indicator. 1 if there is missing data and 0 if not
//
//    OUTPUTS
//      Q: cluster posterior probabilities
//      mu: cluster means (time dependant + variable dependant)
//      s2: cluster variances (variable dependant)
//      theta: cluster priors
//      dim_idx: indexes of the subset of dimension considered
//      time_idx: indexes of the subset of time intervals considered
//
//      Reference: "Time Series Cluster Kernel for Learning Similarities between Multivariate Time Series with Missing Data", 2017 Pattern Recognition, Elsevier.

    val N = x.length
    val T = x(0).length
    val V = x(0)(0).length
    // optional parameters
    assert(minN > 0 && minN <= 1, "The minimum percentage of subsample must be in (0,1]" )
    if (V == 1){
      var minV = 1
    } else{
      var minV = 2
    }
    assert( minV >= 1 && minV <= V, "The minimum number of variables must be in [1,V]")
    assert(maxV >= 1 && maxV <= V, "The maximum number of variables must be in [1,V]")
    assert(minT >= 1 && minT <= T, "The minimum number of variables must be in [1,T]")
    assert(maxT >= 1 && maxT <= T, "The maximum number of variables must be in [1,T]")
    var maxT_new = 0
    maxT_new = min(floor(0.8*T), maxT)

    // Hyperparameters for mean prior (a0, b0) and the std dev prior (n0) of the mixture components
    var a0 : Double = ( 1   - 0.001 ) * random + 0.001
    var b0 : Double = ( 0.2 - 0.005 ) * random + 0.001
    var n0 : Double = ( 0.2 - 0.001 ) * random + 0.001

    //Randomly subsample dimensions, time intervals and samples
    var sN: Int = 0
    if(N > 100){
      sN = round(minN*N).toInt + Random.nextInt(N - round(minN*N).toInt +1)
    } else {
      sN = round(0.9*N).toInt
    }
    var sub_idx: Array[Int] = Array.ofDim(sN)
    sub_idx = Random.shuffle(1 to N).take(sN).sortWith(_<_).toArray

    val sV =  minV + Random.nextInt(maxV - minV +1)
    val dim_idx = Random.shuffle(1 to V).take(sV).sortWith(_<_).toArray//generate sV (sorted) integers between 1 and V

    val t1 =  1 + Random.nextInt(T-minT+1)
    val t2 = (t1 + minT - 1) + Random.nextInt(min(T,(t1 + maxT -1)) - (t1 + minT -1) +1)
    val sT = t2 - t1 +1
    val time_idx: Array[Double] = (t1 to t2).toArray // generate sT continuous integers from t1 to t2

    val sX: Array[Array[Array[Double]]] = Array.ofDim(sub_idx.length, time_idx.length, dim_idx.length)
    for (i <- 0 until sub_idx.length){
      for (j <- 0 until time_idx.length){
        for (k <- 0 until dim_idx.length ){
          sX(i)(j)(k) = x(sub_idx(i))(time_idx(j))(dim_idx(k))
        }
      }
    }

    if (missing == 1 ){
      // handle missing data
      val nan_idx = isNaN(sX)
      val R: Array[Array[Array[Double]]] = Array.fill(N,sV, sT)(1)
      for (i <- 0 until x.length){
        for (j <- 0 until x(i).length){
          for (k <- 0 until x(i)(j).length){
            if (nan_idx(i)(j)(k) == 1 ){
              R(i)(j)(k) = 0
//              sX(i)(j)(k) = -100000
            }
          }
        }
      }

    //Calculate empirical moments
      var mu_0 : Array[Array[Double]] = Array.fill(sT,sV)(0) // prior mean over time and variables (sT x sV)
      for (v <- 0 until sV){

      }
//      for v = 1:sV
//      mu_0(:,v) = nanmean(sX(:,:,v),1);
//      end
//      s_0 = zeros(sV,1); % prior std over variables (sV x 1)
//      tempX = reshape(sX,[sN*sT,sV]);
//      for v = 1:sV
//      s_0(v) = nanstd(tempX(:,v),0,1);
//      end
//      s2_0 = s_0.^2;


    }  else if (missing == 0 ){

    } else {

    }


  }
}