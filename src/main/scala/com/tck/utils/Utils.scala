package com.tck.utils

import breeze.linalg._
import breeze.stats._

import scala.math.{exp, pow, sqrt, Pi}

object Utils{

  def isNaN(x: Array[Array[Array[Double]]]): Array[Array[Array[Double]]] = {
    var bv: Array[Array[Array[Double]]] = Array.ofDim(x.length, x(0).length, x(0)(0).length)
    bv.indices.map(a => bv(a).indices.map(b => bv(a)(b).indices.map(c => {if (x(a)(b)(c).isNaN) bv(a)(b)(c) = 1
                                                                          else bv(a)(b)(c) = 0 } ) ) )
    bv
  }

  def normpdf(x: Double, mu: Double = 0.0, sigma: Double = 1.0 ): Double ={
    //  NORMPDF Normal probability density function (pdf).
    //    Y = NORMPDF(X,MU,SIGMA) returns the pdf of the normal distribution with
    //    mean MU and standard deviation SIGMA, evaluated at the values in X.
    //    The size of Y is the common size of the input arguments.  A scalar
    //    input functions as a constant matrix of the same size as the other inputs.
    //  Default values for MU and SIGMA are 0 and 1 respectively.
    val y: Double = exp(-0.5 * pow(((x - mu) / sigma),2 ) ) / (sqrt(2*Pi) * sigma)
    y
  }
  def nanMean(x : Array[Array[Double]], dim: Int =1 ):  Array[Double] = {
    assert(dim == 0 || dim == 1 , "Dimension is between 0 and 1")
    var m: Array[Double] = Array()
    var check:Double = 0.0
    x.map( a => a.map(b => if(b.isNaN) check += 1.0 ))
    if (check == 0.0){
      val temp = DenseMatrix(x:_*)
      //println(temp)
      if (dim == 0){
        m = mean(temp(::,*)).t.toArray
      } else {
        m = mean(temp(*,::)).toArray
      }
      m
    } else {
      if (dim == 1){
        m = Array.ofDim(x.length)
        var sum = 0.0
        var idx = 0.0
        x.indices.map( a => {sum = 0.0 ; idx = 0.0; x(a).map(b => if (!b.isNaN) {sum += b; idx += 1});
                                  m(a) = sum/idx;})
      } else{
        m = Array.ofDim(x(0).length)
        val sum = Array.fill(x(0).length)(0.0)
        val idx = Array.fill(x(0).length)(0.0)
        x.indices.map( a => {x(a).indices.map(b => {if (!x(a)(b).isNaN) {sum(b) += x(a)(b);idx(b) += 1;}})})
        m.indices.map( a => m(a) = sum(a)/idx(a))
      }
      m
    }
  }

  def nanStd(x : Array[Array[Double]] ,dim: Int = 0): Array[Double]= {
    assert(dim == 0 || dim == 1 , "Dimension is between 0 and 1")
    var m: Array[Double] = Array()
    var std: Array[Double] = Array()
    var check:Double = 0.0
    x.map( a => a.map(b => if(b.isNaN) check += 1.0 ))
    if (check == 0.0){
      val temp = DenseMatrix(x:_*)
      //println(temp)
      if (dim == 0){
        m = mean(temp(::,*)).t.toArray
        std = stddev(temp(::,*)).t.toArray
      } else {
        m = mean(temp(*,::)).toArray
        std = stddev(temp(*,::)).toArray
      }
      std
    } else {
      if (dim == 1){
        m = Array.ofDim(x.length)
        std = Array.ofDim(x.length)
        var sum: Array[Double] = Array.fill(x.length)(0.0)
        var idx: Array[Double] = Array.fill(x.length)(0.0)
        x.indices.map( a => {x(a).map(b => if (!b.isNaN) {sum(a) += b; idx(a) += 1}); m(a) = sum(a)/idx(a);})
        x.indices.map( a => std(a) = sqrt(x(a).map(b => if (b.isNaN) 0 else pow(b - m(a), 2)).sum / (idx(a) - 1)))
      } else{
        m = Array.ofDim(x(0).length)
        std = Array.ofDim(x(0).length)
        val sum = Array.fill(x(0).length)(0.0)
        val idx = Array.fill(x(0).length)(0.0)
        x.indices.map( a => {x(a).indices.map(b => {if (!x(a)(b).isNaN) {sum(b) += x(a)(b);idx(b) += 1;}})})
        m.indices.map( a => m(a) = sum(a)/idx(a))
        std.indices.map(a  => std(a) = sqrt(x.map(_(a)).map(b => if (b.isNaN) 0 else pow(b - m(a), 2)).sum / {if (idx(a) > 1) (idx(a) - 1)  else 1} ) )
      }
      std
    }
  }

  /*def nanMeanStd(x : Array[Array[Double]] ,dim: Int = 0): (Array[Double], Array[Double])= {
    assert(dim == 0 || dim == 1 , "Dimension is between 0 and 1")
    var m: Array[Double] = Array()
    var std: Array[Double] = Array()
    var check:Double = 0.0
    x.map( a => a.map(b => if(b.isNaN) check += 1.0 ))
    if (check == 0.0){
      val temp = DenseMatrix(x:_*)
      //println(temp)
      if (dim == 0){
        m = mean(temp(::,*)).t.toArray
        std = stddev(temp(::,*)).t.toArray
      } else {
        m = mean(temp(*,::)).toArray
        std = stddev(temp(*,::)).toArray
      }
      (m,std)
    } else {
      if (dim == 1){
        m = Array.ofDim(x.length)
        std = Array.ofDim(x.length)
        var sum: Array[Double] = Array.fill(x.length)(0.0)
        var idx: Array[Double] = Array.fill(x.length)(0.0)
        x.indices.map( a => {x(a).map(b => if (!b.isNaN) {sum(a) += b; idx(a) += 1}); m(a) = sum(a)/idx(a);})
        x.indices.map( a => std(a) = sqrt(x(a).map(b => if (b.isNaN) 0 else pow(b - m(a), 2)).sum / (idx(a) - 1)))
      } else{
        m = Array.ofDim(x(0).length)
        std = Array.ofDim(x(0).length)
        val sum = Array.fill(x(0).length)(0.0)
        val idx = Array.fill(x(0).length)(0.0)
        x.indices.map( a => {x(a).indices.map(b => {if (!x(a)(b).isNaN) {sum(b) += x(a)(b);idx(b) += 1;}})})
        m.indices.map( a => m(a) = sum(a)/idx(a))
        std.indices.map(a  => std(a) = sqrt(x.map(_(a)).map(b => if (b.isNaN) 0 else pow(b - m(a), 2)).sum / {if (idx(a) > 1) (idx(a) - 1)  else 1} ) )
      }
      (m,std)
    }
  }*/




}
