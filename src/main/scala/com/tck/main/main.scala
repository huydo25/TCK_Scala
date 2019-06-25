package com.tck.main

import breeze.linalg.{Axis, DenseMatrix, sum}

import scala.io.Source
import com.tck.models.trainTCK._
import com.tck.models.computeKM_TCK._

object Main{

  def readCSV(filename:String) : Array[Array[Double]] = {
    Source.fromFile(filename)
      .getLines()
      .map(_.split(",").map(_.trim.toDouble))
      .toArray
  }
  def main(args: Array[String]): Unit = {
//    val resourcesPath = getClass.getResource("/x_VAR.csv")
//    val x = readCSV(resourcesPath.getPath)
//    val resourcesPath1 = getClass.getResource("/xte_VAR.csv")
//    val xte = readCSV(resourcesPath1.getPath)
//
//    // Reshape raw data into MTS
//    var X : Array[Array[Array[Double]]] = Array.ofDim[Double](200,50,2)
//    for(i <- 0 until x.length ){
//      for (j <- 0 until x(i).length){
//        if ( j < (x(i).length/2)){
//          X(i)(j)(0) = x(i)(j)
//        }
//        else{
//          X(i)(j-x(i).length/2)(1) = x(i)(j)
//        }
//      }
//    }
//    // label data
//    var Y = Array.fill(200)(1)
//    for (i <- 100 until  Y.length  ){
//      Y(i) += 1
//    }
//
//    //println(sum((DenseMatrix(X.map(_.map(_(0))):_*) - DenseMatrix(Array.fill(200)(X.map(_.map(_(0))).map(_(0)).slice(0,50)):_*)).map(x => x*x), Axis._1).t)
//    // Reshape xte data into MTS
//    var Xte : Array[Array[Array[Double]]] = Array.ofDim[Double](200,50,2)
//    for(i <- 0 until xte.length  ){
//      for (j <- 0 until xte(i).length){
//        if ( j < (xte(i).length/2)){
//          Xte(i)(j)(0) = xte(i)(j)
//        }
//        else{
//          Xte(i)(j-x(i).length/2)(1) = xte(i)(j)
//        }
//      }
//    }
//    val Yte =Y

    val resourcesPath = getClass.getResource("/penDigits/pendigits_train.csv")
    val rawx : Array[Array[Double]]= readCSV(resourcesPath.getPath)

    val x = rawx.map(_.slice(0,rawx(0).length-1))
    //println(x.length,x(0).length)
    // Reshape raw data into MTS
    var X : Array[Array[Array[Double]]] = Array.ofDim[Double](x.length,x(0).length/2,2)
    for(i <- 0 until x.length ){
      for (j <- 0 until x(i).length){
        if ( j < (x(i).length/2)){
          X(i)(j)(0) = x(i)(j)
        }
        else{
          X(i)(j-x(i).length/2)(1) = x(i)(j)
        }
      }
    }
    // label data
    var Y :Array[Int] = rawx.map(_(rawx(0).length-1)).map(x => x.toInt)
    //println(X.deep.mkString("\n"))

    // test data
    val resourcesPath1 = getClass.getResource("/penDigits/pendigits_test.csv")
    val rawXte = readCSV(resourcesPath1.getPath)
    val xte = rawXte.map(_.slice(0,rawXte(0).length-1))
    //println(xte.length, xte(0).length)
    // Reshape xte data into MTS
    var Xte : Array[Array[Array[Double]]] = Array.ofDim[Double](xte.length,xte(0).length/2,2)
    for(i <- 0 until xte.length  ){
      for (j <- 0 until xte(i).length){
        if ( j < (xte(i).length/2)){
          Xte(i)(j)(0) = xte(i)(j)
        }
        else{
          Xte(i)(j-xte(i).length/2)(1) = xte(i)(j)
        }
      }
    }
    val Yte :Array[Int] = rawXte.map(_(rawXte(0).length-1)).map(x => x.toInt)
    //println(Yte.length)
    //println(Xte.length, Xte(0).length, Xte(0)(0).length)

    var gmmParameter : List[(Array[Array[Double]], Array[Array[Array[Double]]], Array[Array[Double]],Array[Double], Array[Int], Array[Int])]  = List()
    var C: Int = 0
    var G: Int = 0

    var temp_r = trainTCK(X)
    gmmParameter = temp_r._1
    C = temp_r._2
    G = temp_r._3

    // Compute in-sample kernel matrix
    var K = TCK(gmmParameter,C,G,1)
    println("Done training TCK")
    // Compute similarity between Xte and the training points
    var Kte = TCK(gmmParameter,C,G,0,Xte)

    //println(Kte.length, Kte(0).length)


    // 1NN -classifier
    val Nte = Yte.length
    var I : Array[Int] = Array.ofDim(Kte(0).length)
    for (i <- 0 until Kte(0).length){
      val temp = Kte.map(_(i))
      I(i) = temp.indexOf(temp.max)
    }
    var predY : Array[Int] = Array.fill(I.length)(1)
    I.indices.map(i => predY(i) = Y(I(i)))

    //println(predY.deep.mkString(" "))
    var sum : Int = 0
    for (i <- 0 until Yte.length ){
      if (predY(i) == Yte(i)){
        sum  = sum + 1
      }
    }
    val accuracy : Double = sum.toDouble / Nte.toDouble * 100.0
    println()
    print("Cluster accuracy: ")
    println(accuracy)
    println("Done!!")
  }


}
