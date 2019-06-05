package com.TCK.main

import scala.io.Source
import com.TCK.Models.trainTCK._
import com.TCK.Models.TCK._

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
    for(i <- 0 until x.length ){
      for (j <- 0 until x(i).length){
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
    for (i <- 100 until  Y.length  ){
      Y(i) += 1
    }
    // Reshape xte data into MTS
    var Xte : Array[Array[Array[Double]]] = Array.ofDim[Double](200,50,2)
    for(i <- 0 until xte.length -1 ){
      for (j <- 0 until xte(i).length){
        if ( j < (xte(i).length/2)){
          Xte(i)(j)(0) = xte(i)(j)
        }
        else{
          Xte(i)(j-x(i).length/2)(1) = xte(i)(j)
        }
      }
    }
    val Yte =Y
//    print(Xte(0).deep.mkString("\n"))
    val gmmParameter : List[(Array[Array[Double]], Array[Array[Array[Double]]], Array[Array[Double]],Array[Double], Array[Int], Array[Int])]  = List()
    var C: Int = 0
    var G: Int = 0
    (gmmParameter, C, G) = trainTCK(X)

    // Compute in-sample kernel matrix
    var K = TCK(gmmParameter,C,G,1)

    // % Compute similarity between Xte and the training points
    var Kte = TCK(gmmParameter,C,G,0,Xte);

    // 1NN -classifier
//    Nte = length(Yte);
//    [C,I] = max(Kte);    % find training series with maximum similarity
//      pred_Y = Y(I);     % 1NN classification
//    accuracy=sum(pred_Y==Yte)/Nte



    println("Done!!")
  }


}
