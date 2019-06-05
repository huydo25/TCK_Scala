package com.TCK.Models

import breeze.linalg._
import com.TCK.Ultils.Ultils._

object GMMposterior{

  def GMM_posterior(x : Array[Array[Array[Double]]], C: Int, mu: Array[Array[Array[Double]]],
                    // needed update s2
                    s2: Array[Array[Double]], theta: Array[Double],
                    dim_idx: Array[Int], time_idx: Array[Int], missing: Int): Array[Array[Double]]={
    //    GMMposterior - Evaluate the posterior for the data X of the GMM described by C, mu, s2 and theta
    //    INPUTS
    //      X: data array of size N x T x V
    //      C: number of mixture components (optional)
    //      mu: cluster means over time and variables (V x T)
    //      s2: cluster stds over variables (sV x 1)
    //      theta: cluster priors \
    //      dim_idx: subset of variables to be used in the clustering
    //      time_idx: subset of time intervals to be used in the clustering
    //      missing: binary indicator. 1 if there is missing data and 0 if not
    //    OUTPUTS
    //      Q: posterior
    //    Reference: "Time Series Cluster Kernel for Learning Similarities between Multivariate Time Series with Missing Data", 2017 Pattern Recognition, Elsevier.

    //initialize variables
    val N = x.length
    var Q : Array[Array[Double]]= Array.fill(N,C)(0)
    val sV = dim_idx.length
    val sT = time_idx.length
    //sample of X
    val sX: Array[Array[Array[Double]]] = Array.ofDim(N,sT,sV)
    for (i <- 0 until N){
      for (j <- 0 until sT ){
        for (k <- 0 until sV ){
          sX(i)(j)(k) = x(i)(time_idx(j))(dim_idx(k))
        }
      }
    }
    //end sample
    if (missing == 1){
      // handle missing data
      val nan_idx = isNaN(sX)
      val R: Array[Array[Array[Double]]] = Array.fill(N,sT,sV)(1)
      for (i <- 0 until N){
        for (j <- 0 until sT){
          for (k <- 0 until sV){
            if (nan_idx(i)(j)(k) == 1 ){
              R(i)(j)(k) = 0
              sX(i)(j)(k) = -100000
            }
          }
        }
      }
      // Compute GMM posterior
      var distr_c : Array[Array[Array[Double]]] = Array.ofDim(N, sT, sV)
      for (i <- 0 until  C){
        var temp: Double = 0
        for (j <- 0 until N){
          for (k <- 0 until sT) {
            for (l <- 0 until sV){
              // need to update in term of dimension array
              temp = normpdf(sX(j)(k)(l), mu(k)(l)(i), s2(l)(i))
              if (temp < normpdf(3)) {
                temp = normpdf(3)
              }
              distr_c(j)(k)(l) = temp
            }
          }
        }
        //reshape distribution vector
        var temp1: Array[Array[Double]] = Array.ofDim(N, sV*sT)
        for (l <- 0 until N){
          for (j <- 0 until sT){
            for (k <- 0 until sV){
              temp1(l)(k*sT+j) = distr_c(l)(j)(k)
            }
          }
        }
        // product of distribution function
        var prod_distrc_c : Array[Double] = Array.fill(N)(1)
        for (j <- 0 until N){
          for (k <- 0 until sV*sT){
            prod_distrc_c(j) *= temp1(j)(k)
          }
        }
        for (j <- 0 until N){
          Q(j)(i) = prod_distrc_c(j) * theta(i)
        }
        //end
      }
      val sum = Array.fill(C)(Q.map(_.sum)).transpose
      Q.indices.map(i =>  Q(i).indices.map(j => Q(i)(j)/= sum(i)(j)))

    } else if (missing == 0){
      // Compute GMM posterior
      var distr_c : Array[Array[Array[Double]]] = Array.ofDim(N, sT, sV)
      for (i <- 0 until C){
        var temp: Double = 0
        for (j <- 0 until N){
          for (k <- 0 until sT) {
            for (l <- 0 until sV){
              // need to update in term of dimension array
              temp = normpdf(sX(j)(k)(l), mu(k)(l)(i), s2(l)(i))
              if (temp < normpdf(3)) {
                temp = normpdf(3)
              }
              distr_c(j)(k)(l) = temp
            }
          }
        }
        //reshape distribution vector
        var temp1: Array[Array[Double]] = Array.ofDim(N, sV*sT)
        for (l <- 0 until N){
          for (j <- 0 until sT){
            for (k <- 0 until sV){
              temp1(l)(k*sT+j) = distr_c(l)(j)(k)
            }
          }
        }
        // product of distribution function
        var prod_distrc_c : Array[Double] = Array.fill(N)(1)
        for (j <- 0 until N){
          for (k <- 0 until sV*sT){
            prod_distrc_c(j) *= temp1(j)(k)
          }
        }
        for (j <- 0 until N){
          Q(j)(i) = prod_distrc_c(j) * theta(i)
        }
        //end
      }
      //reshape Q
      val sum = Array.fill(C)(Q.map(_.sum)).transpose
      Q.indices.map(i =>  Q(i).indices.map(j => Q(i)(j)/= sum(i)(j)))

    } else
      sys.error("The value of the variable missing is not 0 or 1")
    Q
  }

}


