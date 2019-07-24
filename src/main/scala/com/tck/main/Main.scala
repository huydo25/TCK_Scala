package com.tck.main

import com.tck.model.TCKModel
import com.tck.utils._


object Main{

  def main(args: Array[String]): Unit = {
    val resourcesPath = getClass.getResource("/x_VAR.csv")
    val x = Reader.readCSV(resourcesPath.getPath)
    val resourcesPath1 = getClass.getResource("/xte_VAR.csv")
    val xte = Reader.readCSV(resourcesPath1.getPath)

    // Reshape raw data into MTS
    var X : Array[Array[Array[Double]]] = Array.ofDim[Double](200,50,2)
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
    var Y = Array.fill(200)(1)
    for (i <- 100 until  Y.length  ){
      Y(i) += 1
    }

    //println(sum((DenseMatrix(X.map(_.map(_(0))):_*) - DenseMatrix(Array.fill(200)(X.map(_.map(_(0))).map(_(0)).slice(0,50)):_*)).map(x => x*x), Axis._1).t)
    // Reshape xte data into MTS
    var Xte : Array[Array[Array[Double]]] = Array.ofDim[Double](200,50,2)
    for(i <- 0 until xte.length  ){
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

    val t1 = System.nanoTime
    val model = new TCKModel(0.8, 2, 2, 6, 18, 30, 30 )
    val train = model.train(X)
    //model.map_em(X,40,0)
    var duration = (System.nanoTime - t1) / 1e9d
    println(duration)
    val t2 = System.nanoTime
    val KM = model.computeKM(train._1, train._2, train._3, 0, Xte )
    duration = (System.nanoTime - t2) / 1e9d
    println(duration)

    /*val Y = Array.fill(200)(1)
    for (i <- 100 until  Y.length  ){
      Y(i) += 1
    }
    val Yte =Y*/

    // 1NN -classifier
    val Nte = Yte.length
    val I : Array[Int] = Array.ofDim(KM(0).length)
    for (i <- 0 until KM(0).length){
      val temp = KM.map(_(i))
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
