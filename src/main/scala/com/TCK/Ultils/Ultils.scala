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
    //    input functions as a constant matrix of the same size as the other inputs.
    //  Default values for MU and SIGMA are 0 and 1 respectively.
    val y: Double = exp(-0.5 * pow(((x - mu) / sigma),2 ) ) / (sqrt(2*Pi) * sigma)
    y
  }

  def nanmean2D(x : Array[Array[Double]], dim: Int = 1): Array[Double]={
    var mean: Array[Double] = Array()
    if (dim == 2 ){
      mean = Array.ofDim(x.length)
      var sum = 0.0
      var idx = 0.0
      for (i <- 0 until x.length){
        sum = 0.0 ; idx = 0.0
        for (j <- 0 until x(i).length){
          if (x(i)(j).isNaN){
            sum += 0
            idx += 0
          } else {
            sum += x(i)(j)
            idx += 1
          }
        }
        //        println (sum, idx, mean.length)
        mean(i) = sum / idx
      }
    } else if (dim == 1){
      mean = Array.ofDim(x(0).length)
      var sum = Array.fill(x(0).length)(0.0)
      var idx = Array.fill(x(0).length)(0)
      for (i <- 0 until x.length){
        for (j <- 0 until x(i).length){
          if (x(i)(j).isNaN){
            sum(j) += 0.0
            idx(j) += 0
          } else {
            sum(j) += x(i)(j)
            idx(j) += 1
          }
        }
      }
      //      println( sum.deep.mkString("\n"))
      //      println( idx.deep.mkString("\n"))
      for (i <- 0 until mean.length){
        mean(i) = sum(i)/idx(i)
      }
    } else {
      sys.error("dimension is between 0 and 1")
    }

    mean
  }

}

