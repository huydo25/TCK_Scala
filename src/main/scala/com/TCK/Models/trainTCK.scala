package com.TCK.Models

object trainTCK{
  def trainTCK(): Unit = {

    // trainTCK - Train the TCK
    //
    //INPUTS
    //   X: data array of size N x T x V, where N is the number of multivariate time series, T the length and V the number of attributes.
    //   minN: min percentage of subsample (optional)
    //   minV: min number of attributes for each GMM (optional)
    //   maxV: max number of attributes for each GMM  (optional)
    //   minT: min length of time segments for each GMM (optional)
    //   maxT: max length of time segments for each GMM (optional)
    //   C: max number of mixture components for each GMM (optional)
    //   G: number of randomizations for each number of components (optional)
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

  }
}
