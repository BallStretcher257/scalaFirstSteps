import breeze.linalg._
import breeze.integrate._

object Main {
  def main(args: Array[String]): Unit = {
    val m = 7
    def phi(s: Double, derivativeOrder: Int = 0) = {
      derivativeOrder match {
        case 0 => DenseVector.tabulate(m){i => math.pow(s, i)}
        case 1 => DenseVector.tabulate(m){i => i * math.pow(s, i-1)}
        case 2 => DenseVector.tabulate(m){i => i * (i-1) * math.pow(s, i-2)}
      }
    }
    def tfunc(c: DenseVector[Double], s: Double, derivativeOrder: Int = 0) = {
      c.t * phi(s, derivativeOrder)
    }
    val int_steps = 1000 // важный параметр
    val s0 = -0.01
    val S = 0.85
    def bmat(c:DenseVector[Double]) = {
      DenseMatrix.tabulate(m, m){case (i, j) => trapezoid(s => phi(s)(i) * phi(s)(j), s0, S, int_steps)}
    }
    val x0 = 1.1 // важный параметр
    val chi = (x0*x0 - 1)/(x0*x0 + 1)
    val epsilon0 = 1 // почему? важный параметр?
    val epsilon   = epsilon0/math.sqrt(1 - chi*chi)
    val d = 2 // важный параметр
    val sigma0 = math.sqrt(2) * d
    def sigma(t: Double) = {
      sigma0/math.sqrt(1 - (t*chi)/(1 - chi*chi))
    }
    def afunc(t: Double, r: Double) = {
      sigma0/(r-sigma(t)+sigma0)
    }
    def rcrit(t : Double) = (math.pow(2, 1/6) - 1) * d + sigma(t)
    def uappr(t: Double, r: Double) = {
      if (r <= rcrit(t))
        4 * epsilon * (math.pow(afunc(t, r), 12) - math.pow(afunc(t, r), 6)) + epsilon
        else 0
    }
    def bigFfunc(t: Double, r: Double) = {
      math.exp(-uappr(t, r))/0.86
    }
    val n = 10
    val rdisc = (BigDecimal(0.9 * d) to BigDecimal(2.1 * d) by 1.2 * d/n).map(_.doubleValue)
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
    val alpha = 0.05 // важный параметр
    val gvec = rdisc.map(math.pow(_, 2))
    def step(c: DenseVector[Double]) = {
      inv(jacmat(c).t * jacmat(c) + alpha * bmat(c))
        * (jacmat(c).t * (avec(c) - gvec) + alpha * bmat(c) * c)
    }
  }
}