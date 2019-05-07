package com.TCK.Model

import com.TCK.Ultils.Ultils
import breeze.linalg._

object GMM_posterior{

  def GMM_posterior(X : Array[Array[Array[Double]]], C: Int, mu: Array[Array[Array[Double]]],
                    s2: Array[Array[Double]], theta: Array[Double],
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
    var Q = Array.fill(N,C)(0)
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
      var nan_idx = Ultils.isNaN(sX)
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
      for (i <- 0 to C-1){

        //        distr_c = normpdf(sX, permute(repmat(mu(:,:,c),[1,1,N]),[3,1,2]), permute(repmat(sqrt(s2(:,c)),[1,N,sT]),[2,3,1]) ).^R;
        //        distr_c(distr_c < normpdf(3)) = normpdf(3);
        //        distr_c = reshape(distr_c,[N,sV*sT]);
        //        Q(:,c) = theta(c)*prod(distr_c,2);

      }
    }
  }
}


