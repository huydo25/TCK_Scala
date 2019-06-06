package com.tck.models

import breeze.linalg._
import breeze.numerics.floor
import com.tck.utils.Utils._
import com.tck.models.GMMposterior._

object TCK {
  def TCK( gmmParameter : List[ (Array[Array[Double]], Array[Array[Array[Double]]], Array[Array[Double]],Array[Double], Array[Int], Array[Int])] ,
           C: Int, G: Int, marking: Int = 0 , xte: Array[Array[Array[Double]]] = Array()): Array[Array[Double]] ={
  // function [ K ] = TCK(GMM, C, G, Xte)
  // TCK -  compute TCK kernel matrix between training data and test data Xte
  // INPUTS
  //   GMM : Cell output from the function trainTCK
  //    C: Second output from trainTCK
  //    G: Third output from trainTCK
  //    marking: mark to compute in-sample kernel matrix
  //    Xte: data array of size Nte x T x V, where Nte is the number of
  //      multivariate time series, T the length and V the number of attributes.
  // OUTPUTS
  //    K: kernel matrix

  if (marking == 0 ){
    // Check if the dataset contains mising elements
    val nan_idx = isNaN(xte)
    var sum = 0
    for (i <- 0 until nan_idx(0).length){
      for (j <- 0 until nan_idx(0)(0).length ){
        sum += nan_idx.map(_(i)).map(_(j)).sum.toInt
      }
    }
    var missing = 0
    if (sum > 0 ){
      missing = 1
      //println("The dataset contains missing data\n\n")
    } else {
      missing = 0
      //println("The dataset does not contain missing data\n\n")
    }
    var K : Array[Array[Double]] = Array.fill(gmmParameter(0)._1.length,xte.length)(0)
    var result: DenseMatrix[Double] = DenseMatrix.zeros(gmmParameter(0)._1.length,xte.length)
    for (i <- 0 until G*(C-1)){
      val c = floor((i-1)/G) +2
      val q : Array[Array[Double]] = GMMposterior.GMM_posterior( xte, c, gmmParameter(i)._2, gmmParameter(i)._3,
                                              gmmParameter(i)._4, gmmParameter(i)._5, gmmParameter(i)._6, missing)
      result = result + DenseMatrix(gmmParameter(i)._1:_*) * DenseMatrix(q:_*).t
    }
    for (k <- 0 until result.rows ){
      K(k) = result(::,k).toArray
    }
    K
  }  else { // in-sample kernel matrix
    var K : Array[Array[Double]] = Array.fill(gmmParameter(0)._1.length,gmmParameter(0)._1.length)(0)
    var result: DenseMatrix[Double] = DenseMatrix.zeros(gmmParameter(0)._1.length,gmmParameter(0)._1.length)
    for (i <- 0 until G*(C-1)){
      result = result + DenseMatrix(gmmParameter(i)._1:_*) * DenseMatrix(gmmParameter(i)._1:_*).t
    }
    for (k <- 0 until result.rows ){
      K(k) = result(::,k).toArray
    }
    K
  }


  }
}
