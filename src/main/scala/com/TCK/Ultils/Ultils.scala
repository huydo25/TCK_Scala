package com.TCK.Ultils

import scala.math._

object Ultils{

  def isNaN(x: Array[Array[Array[Double]]]) : Array[Array[Array[Double]]] =  {
    var binary_vector: Array[Array[Array[Double]]] = Array.ofDim(x.length, x(0).length, x(0)(1).length)
    for (i <- 0 until x.length){
      for (j <- 0 until  x(i).length){
        for (k <- 0 until  x(i)(j).length){
          if (x(i)(j)(k).isNaN){
            binary_vector(i)(j)(k) = 1
          } else{
            binary_vector(i)(j)(k) = 0
          }
        }
      }
    }
    binary_vector
  }

  def normpdf(x: Double, mu: Double = 0, sigma: Double =1 ): Double ={
    //  NORMPDF Normal probability density function (pdf).
    //    Y = NORMPDF(X,MU,SIGMA) returns the pdf of the normal distribution with
    //    mean MU and standard deviation SIGMA, evaluated at the values in X.
    //    The size of Y is the common size of the input arguments.  A scalar
    //    input functions as a constant matrix of the same size as the other
    //    inputs.
    //  Default values for MU and SIGMA are 0 and 1 respectively.
    val y: Double = exp(-0.5 * pow(((x - mu) / sigma),2 ) ) / (sqrt(2*Pi) * sigma)

    y

  }
}

