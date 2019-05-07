package com.TCK.Model

import com.TCK.Ultils.Ultils
import breeze.linalg._
import com.TCK.Ultils.Ultils._

object GMM_posterior{

  def GMM_posterior(X : Array[Array[Array[Double]]], C: Int, mu: Array[Array[Array[Double]]],
                    // needed update s2
                    s2: Array[Array[Array[Double]]], theta: Array[Double],
                    dim_idx: Array[Double], time_idx: Array[Double], missing: Int): Unit ={
    //    GMMposterior - Evaluate the posterior for the data X of the GMM described by C, mu, s2 and theta
    //    INPUTS
    //      X: data array of size N x V x T
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
    val N = X.length
    var Q : Array[Array[Double]]= Array.fill(N,C)(0)
    val sV = dim_idx.length
    val sT = time_idx.length
    //sample of X
    var sX: Array[Array[Array[Double]]] = Array.ofDim(N,sV, sT)
    for (i <- 0 to sX.length){
      for (j <- 0 to sX(i).length){
        for (k <- 0 to sX(i)(j).length){
          sX(i)(j)(k) = X(i, time_idx(j), dim_idx(k))
        }
      }
    }
    //end sample

    if (missing == 1 ){
      // handle missing data
      var nan_idx = isNaN(sX)
      var R: Array[Array[Array[Double]]] = Array.fill(N,sV, sT)(1)
      for (i <- 0 to X.length-1){
        for (j <- 0 to X(i).length-1){
          for (k <- 0 to X(i)(j).length-1){
            if (nan_idx(i)(j)(k) == 1 ){
              R(i)(j)(k) = 0
              sX(i)(j)(k) = -100000
            }
          }
        }
      }
      // Compute GMM posterior
      var distr_c : Array[Array[Array[Double]]] = Array.ofDim(N, sV, sT)
      for (i <- 0 to C-1){
        var temp: Double = 0
        for (j <- 0 to N){
          for (k <- 0 to sV) {
            for (l <- 0 to sT){
              // need to update in term of dimension array
              temp = normpdf(sX(j)(k)(l), mu(j)(k)(l), s2(j)(k)(l))
              if (temp < normpdf(3)) {
                temp = normpdf(3)
              }
              distr_c(j)(k)(l) = temp
            }
          }
        }
        //reshape distribution vector
        var temp1: Array[Array[Double]] = Array.ofDim(N, sV*sT)
        for (l <- 0 to N){
          for (j <- 0 to sV){
            for (k <- 0 to sT){
              temp1(l)(j*sT+k) = distr_c(l)(j)(k)
            }
          }
        }
        // product of distribution function
        var prod_distrc_c : Array[Double] = Array.fill(N)(1)
        for (j <- 0 to N){
          for (k <- 0 to sV*sT){
            prod_distrc_c(j) *= temp1(j)(k)
          }
        }
        for (j <- 0 to N){
          Q(j)(i) = prod_distrc_c(j) * theta(i)
        }
        //end
      }
      //reshape Q
      // Q = Q./repmat(sum(Q,2),[1,C]);
    } else if (missing == 0){
      // Compute GMM posterior
      var distr_c : Array[Array[Array[Double]]] = Array.ofDim(N, sV, sT)
      for (i <- 0 to C-1){
        var temp: Double = 0
        for (j <- 0 to N){
          for (k <- 0 to sV) {
            for (l <- 0 to sT){
              // need to update in term of dimension array
              temp = normpdf(sX(j)(k)(l), mu(j)(k)(l), s2(j)(k)(l))
              if (temp < normpdf(3)) {
                temp = normpdf(3)
              }
              distr_c(j)(k)(l) = temp
            }
          }
        }
        //reshape distribution vector
        var temp1: Array[Array[Double]] = Array.ofDim(N, sV*sT)
        for (l <- 0 to N){
          for (j <- 0 to sV){
            for (k <- 0 to sT){
              temp1(l)(j*sT+k) = distr_c(l)(j)(k)
            }
          }
        }
        // product of distribution function
        var prod_distrc_c : Array[Double] = Array.fill(N)(1)
        for (j <- 0 to N){
          for (k <- 0 to sV*sT){
            prod_distrc_c(j) *= temp1(j)(k)
          }
        }
        for (j <- 0 to N){
          Q(j)(i) = prod_distrc_c(j) * theta(i)
        }
        //end
      }
      //reshape Q
      // Q = Q./repmat(sum(Q,2),[1,C]);

    } else
      sys.error("The value of the variable missing is not 0 or 1")
  }
}


