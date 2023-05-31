import math._
import breeze.linalg._
import breeze.integrate._
import breeze.plot._
object Main {
  def main(args: Array[String]): Unit = {
    val m = 7
    def phiOld(s: Double, derivativeOrder: Int = 0) = {
      DenseVector.tabulate(m){
        i =>
          derivativeOrder match {
            case 0 => pow(s, i)
            case 1 => i * pow(s, i-1)
            case 2 => i * (i-1) * pow(s, i-2)
          }
      }
    }
    def binom(n: Int, k: Int): Double = {
      require(0 <= k && k <= n)
      @annotation.tailrec
      def binomtail(nIter: Int, kIter: Int, ac: Double): Double = {
        if (kIter > k) ac
        else binomtail(nIter + 1, kIter + 1, (nIter * ac) / kIter)
      }
      if (k == 0 || k == n) 1
      else binomtail(n - k + 1, 1, 1.0)
    }
    val s0 = -0.01
    val S = 0.85
    def phi(s: Double, derivativeOrder: Int = 0) = {
      DenseVector.tabulate(m){
        i => {
          val x = 2*(s - s0)/(S - s0) - 1
          (0 to i).map(k => binom(i, k) * binom(k + i, k) * pow((x - 1)/2, k)).sum
        }
      }
    }

    def tfunc(c: DenseVector[Double], s: Double, derivativeOrder: Int = 0) = {
      c.t * phi(s, derivativeOrder)
    }
    val int_steps = 50 // важный параметр
    //val bmat = DenseMatrix.tabulate(m, m){case (i, j) => trapezoid(s => phi(s)(i) * phi(s)(j), s0, S, int_steps)}
    val bmat = DenseMatrix.tabulate(m, m){(i, j) => if (i == j) 0.86/(2*i + 1) else 0}
    val x0 = 1.6 // важный параметр
    val chi = (x0*x0 - 1)/(x0*x0 + 1)
    val epsilon0 = 1 // почему? важный параметр?
    val epsilon   = epsilon0/sqrt(1 - chi*chi)
    val d = 2 // важный параметр
    val sigma0 = sqrt(2) * d
    def sigma(t: Double) = {DenseMatrix.eye[Double](m)
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
    val n = 50
    val rdisc = DenseVector.tabulate(n){i => 0.9 * d + 1.2 * d * i / n}
    def jacmatByParts(c: DenseVector[Double]) = {
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
    var stepnum = 0
    def nanprobe(arg: Double) = {
      stepnum = stepnum + 1
      if (arg.isNaN) throw new Exception("nan found at step " + stepnum) else println(arg)
      arg
    }
    def jacmatAnal(c: DenseVector[Double]) = {
      DenseMatrix.tabulate(n, m) { case (i, j) =>
        trapezoid(
          s => {
            val t = tfunc(c, s)
            val r = rdisc(i)
            if (r <= rcrit(t))
              bigFfunc(t, r) * 4 * epsilon * (12*pow(afunc(t, r), 11) - 6 * pow(afunc(t, r), 5)) *
              sigma0 / ((r-sigma(t)+sigma0) * (r-sigma(t)+sigma0)) *
              sigma0 * pow(1 - (t*chi)/(1 - chi*chi), - 3 / 2) * chi / (chi*chi - 1)
            else 0
          }, s0, S, int_steps)
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
      trapezoid(f =>
        trapezoid(z =>
          trapezoid(t =>
            bigFfuncOriginal(t, z, f, r),
            0, Pi, int_steps),
          0, Pi, int_steps),
        0, 2 * Pi, int_steps))
    val alpha = 0.05  // важный параметр
    def step(c: DenseVector[Double], jacmat: DenseVector[Double] => DenseMatrix[Double]) = {
      inv(jacmat(c).t * jacmat(c) + alpha * bmat) * (jacmat(c).t * (avec(c) - gvec) + alpha * bmat * c)
    }
    val c = DenseVector.ones[Double](m)
    for (_ <- 1 to 100) {
      c -= step(c, jacmatAnal)
    }
    val f = Figure()
    val p = f.subplot(0)
    p += plot(rdisc, gvec, '.')
    p += plot(rdisc, avec(c))
  }
}