#ifndef OPENMC_BRENT_H
#define OPENMC_BRENT_H

// ============================================================
// Brent's method for bracketed root finding.
//
// Original source: John Burkardt, C++ version of:  (not available in pdf, but this is the original source)
//   Richard Brent, "Algorithms for Minimization Without Derivatives",
//   Dover, 2002, ISBN: 0-486-41998-3.
//   Function: zero(), distributed under GNU LGPL license.
//   URL: https://people.sc.fsu.edu/~jburkardt/classes/isc_2011/week11/brent.cpp

#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

// brent function defined as template!! like it was done with bisect() and toms748()
template<class F, class Tol>
std::pair<double, double> brent(F f, double min, double max, Tol tol, std::uintmax_t& max_iter)
{
  // sa, sb = current bracket endpoints (sa=left, sb=right=best estimate)
  // c      = previous sb, such that f(sb) and f(c) have opposite signs
  // d      = last step
  // e      = step before last
  double sa = min;
  double sb = max;
  double fa = f(sa);
  double fb = f(sb);

  double c  = sa;
  double fc = fa;
  double e  = sb - sa;
  double d  = e;

  // Exact roots at endpoints cases
  if (fa == 0.0) { max_iter = 1; return {sa, sa}; }
  if (fb == 0.0) { max_iter = 1; return {sb, sb}; }

  // Machine tolerance ---> 2.22*10^(-16)
  const double macheps = std::numeric_limits<double>::epsilon();

  // Use while cycle instead of for to loop on max number of iterations
  while (max_iter > 0 && !tol(sa, sb)) {

    // sb must always be the best estimate: if c is better, swap!!
    if (std::abs(fc) < std::abs(fb)) {
      sa = sb;
      sb = c;
      c  = sa;
      fa = fb;
      fb = fc;
      fc = fa;
    }

    // half-interval and machine-precision tolerance 
    double t    = 1.e-8;                // tolerance
    double tol1 = 2.0 * macheps * std::abs(sb) + t;
    double m    = 0.5 * (c - sb);        // half interval thickness

    // check on machine-precision convergence  
    if (std::abs(m) <= tol1 || fb == 0.0)
      break;

    // try interpolation if previous step was large enough and sb is better than sa
    if (std::abs(e) < tol1 || std::abs(fa) <= std::abs(fb)) {
      // convergence too slow: bisection
      e = m;
      d = e;
    } 
    else {
      double s = fb / fa;
      double p, q, r;

      if (sa == c) {       // only two distinct points: secant method
        p = 2.0 * m * s;
        q = 1.0 - s;
      } 
      else {               // three distinct points: inverse quadratic interpolation
        q = fa / fc;
        r = fb / fc;
        p = s * (2.0 * m * q * (q - r) - (sb - sa) * (r - 1.0));
        q = (q - 1.0) * (r - 1.0) * (s - 1.0);
      }

      // p must be positive for the acceptance check
      if ( p > 0.0 )
      {
        q = - q;
      }
      else
      {
        p = - p;
      }

      s = e;
      e = d;

      // accept interpolated step only if within bracket and faster than last step
      if (2.0 * p < 3.0 * m * q - std::abs(tol1 * q) && p < std::abs(0.5 * s * q)) {    // accept interpolation
        d = p / q;    
      } 
      else {          // fall back to bisection
        e = m;        
        d = e;
      }
    }

    // advance sb by the chosen step
    sa = sb;
    fa = fb;

    if (std::abs(d) > tol1) 
    {
      sb = sb + d;
    }
    else if ( 0.0 < m )
    {
      sb = sb + tol1;
    }
    else
    {
      sb = sb - tol1;
    }
    
    fb = f(sb);

    // ensure [sb, c] still brackets the root; if not, replace c with old sa
    if ((fb > 0.0 && fc > 0.0) || (fb <= 0.0 && fc <= 0.0)) {
      c  = sa;
      fc = fa;
      e  = sb - sa;
      d  = e;
    }

    --max_iter;
  }

  // return bracket ordered as (smaller, larger) (modification 4)
  return (sb < c) ? std::make_pair(sb, c) : std::make_pair(c, sb);
}

#endif



