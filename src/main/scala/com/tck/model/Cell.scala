package com.tck.model

class Cell(q: Array[Array[Double]], mu: Array[Array[Array[Double]]],
           s2: Array[Array[Double]], theta: Array[Double], dim_idx: Seq[Int], time_id: Seq[Int]) {
  //(Array[Array[Double]], Array[Array[Array[Double]]], Array[Array[Double]],Array[Double], Array[Int], Array[Int])
  val Q = q
  val mean = mu
  val std = s2
  val Th = theta
  val D = dim_idx
  val T = time_id
}
