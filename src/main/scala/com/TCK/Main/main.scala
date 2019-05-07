package com.TCK.main

import scala.io.Source


object Main{

  def readCSV(filename:String) : Array[Array[Double]] = {
    Source.fromFile(filename)
      .getLines()
      .map(_.split(",").map(_.trim.toDouble))
      .toArray
  }
  def main(args: Array[String]): Unit = {
    val resourcesPath = getClass.getResource("/x_VAR.csv")
    val x = readCSV(resourcesPath.getPath)
    val resourcesPath1 = getClass.getResource("/xte_VAR.csv")
    val xte = readCSV(resourcesPath1.getPath)
    //val x_new = DenseMatrix(x:_*)
    // println(x.length, x(0).length)

    // Reshape raw data into MTS
    var X : Array[Array[Array[Double]]] = Array.ofDim[Double](200,50,2)
    for(i <- 0 to x.length -1 ){
      for (j <- 0 to x(i).length-1){
        if ( j < (x(i).length/2)){
          X(i)(j)(0) = x(i)(j)
          //println(j, x(i)(j))
        }
        else{
          X(i)(j-x(i).length/2)(1) = x(i)(j)
          //println(j-(x(i).length/2),x(i)(j))
        }
      }
    }
    // label data
    var Y = Array.fill(200)(1)
    for (i <- 100 to Y.length -1 ){
      Y(i) += 1
    }
    // Reshape xte raw data into MTS
    var Xte : Array[Array[Array[Double]]] = Array.ofDim[Double](200,50,2)
    for(i <- 0 to xte.length -1 ){
      for (j <- 0 to xte(i).length-1){
        if ( j < (xte(i).length/2)){
          Xte(i)(j)(0) = xte(i)(j)
        }
        else{
          Xte(i)(j-x(i).length/2)(1) = xte(i)(j)
        }
      }
    }
    val Yte =Y
    print(Xte(0).deep.mkString("\n"))
  }

}
