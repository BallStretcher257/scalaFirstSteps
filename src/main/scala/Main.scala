import math._
import breeze.linalg._
import breeze.integrate._
object Main {
  def main(args: Array[String]): Unit = {
    val m = 7
    def phi(s: Double, derivativeOrder: Int = 0) = {
      derivativeOrder match {
        case 0 => DenseVector.tabulate(m){i => pow(s, i)}
        case 1 => DenseVector.tabulate(m){i => i * pow(s, i-1)}
        case 2 => DenseVector.tabulate(m){i => i * (i-1) * pow(s, i-2)}
      }
    }
    def tfunc(c: DenseVector[Double], s: Double, derivativeOrder: Int = 0) = {
      c.t * phi(s, derivativeOrder)
    }
    val int_steps = 50 // важный параметр
    val s0 = -0.01
    val S = 0.85
    val bmat = DenseMatrix.tabulate(m, m){case (i, j) => trapezoid(s => phi(s)(i) * phi(s)(j), s0, S, int_steps)}
    val x0 = 2.3 // важный параметр
    val chi = (x0*x0 - 1)/(x0*x0 + 1)
    val epsilon0 = 1 // почему? важный параметр?
    val epsilon   = epsilon0/sqrt(1 - chi*chi)
    val d = 2 // важный параметр
    val sigma0 = sqrt(2) * d
    def sigma(t: Double) = {
      sigma0/sqrt(1 - (t*chi)/(1 - chi*chi))
    }
    def afunc(t: Double, r: Double) = {
      sigma0/(r-sigma(t)+sigma0)
    }
    def rcrit(t : Double) = (pow(2, 1/6) - 1) * d + sigma(t)
    def uappr(t: Double, r: Double) = {
      if (r <= rcrit(t))
        4 * epsilon * (pow(afunc(t, r), 12) - pow(afunc(t, r), 6)) + epsilon
        else 0
    }
    def bigFfunc(t: Double, r: Double) = {
      exp(-uappr(t, r))/0.86
    }
    val n = 10
    val rdisc = DenseVector.tabulate(n){i => 0.9 * d + 1.2 * d * i / n}
    def jacmat(c: DenseVector[Double]) = {
      DenseMatrix.tabulate(n, m) { case (i, j) =>
        def temp_func(s: Double) = bigFfunc(tfunc(c, s), rdisc(i)) * phi(s)(j)/tfunc(c, s, 1)
        temp_func(S) - temp_func(s0) - trapezoid(
          s =>
          bigFfunc(tfunc(c, s), rdisc(i)) *
            (phi(s, 1)(j) * tfunc(c, s, 1)
              - phi(s)(j) * tfunc(c, s, 2))/
            (tfunc(c, s, 1)*tfunc(c, s, 1))
          , s0, S, int_steps)
      }
    }
    def avec(c: DenseVector[Double]) = {
      rdisc.map(r => trapezoid(s => bigFfunc(tfunc(c, s), r), s0, S, int_steps))
    }
    def gfunc(t: Double, z: Double, f: Double) = {
      (1 + cos(t) * cos(t) * (1 - 2 * chi)) * cos(z) * cos(z) +
        (1.0/2 * (1 - chi) * sin(2 * t) * sin(2 * z) +
          sin(t)*sin(t)*sin(z)*sin(z)*cos(f))*cos(f)
    }
    def sigmaOriginal(t: Double, z: Double, f: Double) = {
      d / sqrt(1 - chi * gfunc(t: Double, z: Double, f: Double) / (1 - chi * chi * cos(t) * cos(t)))
    }
    def afuncOriginal(t: Double, z: Double, f: Double, r: Double) = {
      d / (r - sigmaOriginal(t: Double, z: Double, f: Double) + d)
    }
    def epsilonOriginal(t: Double) = {
      epsilon0 / sqrt(1 - chi*chi*cos(t)*cos(t))
    }
    def rcritOriginal(t: Double, z: Double, f: Double) = {
      (pow(2, 1/6) - 1) * d + sigmaOriginal(t: Double, z: Double, f: Double)
    }
    def uOriginal(t: Double, z: Double, f: Double, r: Double) = {
      if (r <= rcritOriginal(t, z, f))
        4 * epsilonOriginal(t) * (pow(afuncOriginal(t, z, f, r), 12) - pow(afuncOriginal(t, z, f, r), 6)) + epsilonOriginal(t)
      else 0
    }
    def bigFfuncOriginal(t: Double, z: Double, f: Double, r: Double) = {
      exp(- uOriginal(t, z, f, r)) * sin(z) * sin (t) / (8 * Pi)
    }
    val gvec = rdisc.map(r =>
      trapezoid(t =>
        trapezoid(z =>
          trapezoid(f =>
            bigFfuncOriginal(t, z, f, r),
            0, Pi, int_steps),
          0, Pi, int_steps),
        0, 2 * Pi, int_steps))
    val alpha = 0.05 // важный параметр
    def step(c: DenseVector[Double]) = {
      inv(jacmat(c).t * jacmat(c) + alpha * bmat) * (jacmat(c).t * (avec(c) - gvec) + alpha * bmat * c)
    }
    var c = DenseVector.ones[Double](m)
    for (i <- 1 to 100) {
      c = c - step(c)
      println(c)
    }
  }
}