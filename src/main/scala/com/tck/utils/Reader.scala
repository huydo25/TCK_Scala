package com.tck.utils


import scala.io.Source

object Reader{
  def readCSV(filename:String) : Array[Array[Double]] = {
    Source.fromFile(filename)
      .getLines()
      .map(_.split(",").map(_.trim.toDouble))
      .toArray
  }
}
