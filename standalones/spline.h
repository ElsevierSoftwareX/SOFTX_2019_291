#ifndef SPLINE_HPP
#define SPLINE_HPP

#include <cstdio>
#include <cstdlib>
#include <limits>
#include "mt_modes_base.h"

inline mt_real b3_spline(const mt_real x)
{
    mt_real absx = std::abs(x);
    if (absx < 1)
    {
        mt_real y = 2 - absx;
        mt_real z = 1 - absx;
        return (1. / 6.) *(y*y*y - 4*z*z*z);
    }
    if (absx < 2)
    {
        mt_real y = 2 - absx;
        return (1. / 6.) * y*y*y;
    }
    return (mt_real)(0);
}

inline mt_real b3_spline_prime(const mt_real x)
{
    if (x < 0)
    {
        return -b3_spline_prime(-x);
    }

    if (x < 1)
    {
        return x*(1.5*x - 2);
    }
    if (x < 2)
    {
        return -0.5*(2 - x)*(2 - x);
    }
    return (mt_real)(0.0);
}

// spline interpolation
class spline
{
private:
    mt_vector<mt_real> m_beta;
    mt_real m_h_inv;
    mt_real m_a;
    mt_real m_avg;
public:
    // set default boundary condition to be zero curvature at both ends
    spline() : m_beta(), m_h_inv(0.0), m_a(0.0), m_avg(0.0)
    {
        ;
    }
    template<class InputIterator>
    spline(InputIterator s, InputIterator e, const mt_real left_endpoint,
            const mt_real step_size, const mt_real left_endpoint_derivative = std::numeric_limits<mt_real>::quiet_NaN(),
            const mt_real right_endpoint_derivative = std::numeric_limits<mt_real>::quiet_NaN());

    mt_real operator() (const mt_real x) const;
    mt_real prime(const mt_real x) const;
};


template<class InputIterator>
spline::spline(InputIterator s, InputIterator e, const mt_real left_endpoint,
            const mt_real step_size, const mt_real left_endpoint_derivative,
            const mt_real right_endpoint_derivative) : m_a(left_endpoint), m_avg(0)
{
  const mt_real third = 1. / 3.;
  size_t length = std::distance(s,e);
  if (length < 5)
  {
      if (std::isnan(left_endpoint_derivative) || std::isnan(right_endpoint_derivative))
      {
          fprintf(stderr, "ERROR in %s: Interpolation using a cubic b spline with derivatives estimated at the endpoints requires at least 5 points.\n", __func__);
          exit(EXIT_FAILURE);
      }
      if (length < 3)
      {
          fprintf(stderr, "ERROR in %s: Interpolation using a cubic b spline requires at least 3 points.\n", __func__);
          exit(EXIT_FAILURE);
      }
  }
  if (std::isnan(left_endpoint))
  {
      fprintf(stderr, "ERROR in %s: Left endpoint is NAN; this is disallowed.\n", __func__);
      exit(EXIT_FAILURE);
  }
  if (left_endpoint + length*step_size >= (std::numeric_limits<mt_real>::max)())
  {
      fprintf(stderr, "ERROR in %s: Right endpoint overflows the maximum representable number of the specified precision.\n", __func__);
      exit(EXIT_FAILURE);
  }
  if (step_size <= 0)
  {
      fprintf(stderr, "ERROR in %s: The step size must be strictly > 0.\n", __func__);
      exit(EXIT_FAILURE);
  }

  m_h_inv = 1 / step_size;
  // Following Kress's notation, s'(a) = a1, s'(b) = b1
  mt_real a1 = left_endpoint_derivative;
  // See the finite-difference table on Wikipedia for reference on how
  // to construct high-order estimates for one-sided derivatives:
  // https://en.wikipedia.org/wiki/Finite_difference_coefficient#Forward_and_backward_finite_difference
  // Here, we estimate then to O(h^4), as that is the maximum accuracy we could obtain from this method.
  if (std::isnan(a1))
  {
      // For simple functions (linear, quadratic, so on)
      // almost all the error comes from derivative estimation.
      // This does pairwise summation which gives us another digit of accuracy over naive summation.
      mt_real t0 = 4*(s[1] + third*s[3]);
      mt_real t1 = -(25*third*s[0] + s[4])/4  - 3*s[2];
      a1 = m_h_inv*(t0 + t1);
  }
  mt_real b1 = right_endpoint_derivative;
  if (std::isnan(b1))
  {
      size_t n = length - 1;
      mt_real t0 = 4*(s[n-3] + third*s[n - 1]);
      mt_real t1 = -(25*third*s[n - 4] + s[n])/4  - 3*s[n - 2];
      b1 = m_h_inv*(t0 + t1);
  }
  // s(x) = \sum \alpha_i B_{3}( (x- x_i - a)/h )
  // Of course we must reindex from Kress's notation, since he uses negative indices which make C++ unhappy.
  m_beta.assign(length + 2, std::numeric_limits<mt_real>::quiet_NaN());
  // Since the splines have compact support, they decay to zero very fast outside the endpoints.
  // This is often very annoying; we'd like to evaluate the interpolant a little bit outside the
  // boundary [a,b] without massive error.
  // A simple way to deal with this is just to subtract the DC component off the signal, so we need the average.
  // This algorithm for computing the average is recommended in
  // http://www.heikohoffmann.de/htmlthesis/node134.html
  mt_real t = 1;
  for (size_t i = 0; i < length; ++i)
  {
      if (std::isnan(s[i]))
      {
          fprintf(stderr, "ERROR in %s: This function you are trying to interpolate is a nan at index %zd\n", __func__, i);
          exit(EXIT_FAILURE);
      }
      m_avg += (s[i] - m_avg) / t;
      t += 1;
  }
  // Now we must solve an almost-tridiagonal system, which requires O(N) operations.
  // There are, in fact 5 diagonals, but they only differ from zero on the first and last row,
  // so we can patch up the tridiagonal row reduction algorithm to deal with two special rows.
  // See Kress, equations 8.41
  // The the "tridiagonal" matrix is:
  // 1  0 -1
  // 1  4  1
  //    1  4  1
  //       1  4  1
  //          ....
  //          1  4  1
  //          1  0 -1
  // Numerical estimate indicate that as N->Infinity, cond(A) -> 6.9, so this matrix is good.
  mt_vector<mt_real> rhs(length + 2, std::numeric_limits<mt_real>::quiet_NaN());
  mt_vector<mt_real> super_diagonal(length + 2, std::numeric_limits<mt_real>::quiet_NaN());

  rhs[0] = -2*step_size*a1;
  rhs[rhs.size() - 1] = -2*step_size*b1;

  super_diagonal[0] = 0;

  for(size_t i = 1; i < rhs.size() - 1; ++i)
  {
      rhs[i] = 6*(s[i - 1] - m_avg);
      super_diagonal[i] = 1;
  }


  // One step of row reduction on the first row to patch up the 5-diagonal problem:
  // 1 0 -1 | r0
  // 1 4 1  | r1
  // mapsto:
  // 1 0 -1 | r0
  // 0 4 2  | r1 - r0
  // mapsto
  // 1 0 -1 | r0
  // 0 1 1/2| (r1 - r0)/4
  super_diagonal[1] = 0.5;
  rhs[1] = (rhs[1] - rhs[0])/4;

  // Now do a tridiagonal row reduction the standard way, until just before the last row:
  for (size_t i = 2; i < rhs.size() - 1; ++i)
  {
      mt_real diagonal = 4 - super_diagonal[i - 1];
      rhs[i] = (rhs[i] - rhs[i - 1])/diagonal;
      super_diagonal[i] /= diagonal;
  }

  // Now the last row, which is in the form
  // 1 sd[n-3] 0      | rhs[n-3]
  // 0  1     sd[n-2] | rhs[n-2]
  // 1  0     -1      | rhs[n-1]
  mt_real final_subdiag = -super_diagonal[rhs.size() - 3];
  rhs[rhs.size() - 1] = (rhs[rhs.size() - 1] - rhs[rhs.size() - 3])/final_subdiag;
  mt_real final_diag = -1/final_subdiag;
  // Now we're here:
  // 1 sd[n-3] 0         | rhs[n-3]
  // 0  1     sd[n-2]    | rhs[n-2]
  // 0  1     final_diag | (rhs[n-1] - rhs[n-3])/diag

  final_diag = final_diag - super_diagonal[rhs.size() - 2];
  rhs[rhs.size() - 1] = rhs[rhs.size() - 1] - rhs[rhs.size() - 2];


  // Back substitutions:
  m_beta[rhs.size() - 1] = rhs[rhs.size() - 1]/final_diag;
  for(size_t i = rhs.size() - 2; i > 0; --i)
  {
      m_beta[i] = rhs[i] - super_diagonal[i]*m_beta[i + 1];
  }
  m_beta[0] = m_beta[2] + rhs[0];
}

mt_real spline::operator()(const mt_real x) const
{
  // See Kress, 8.40: Since B3 has compact support, we don't have to sum over all terms,
  // just the (at most 5) whose support overlaps the argument.
  mt_real z = m_avg;
  mt_real t = m_h_inv*(x - m_a) + 1;

  size_t k_min = (size_t) (std::max)(static_cast<long>(0), static_cast<long>(std::trunc(std::ceil(t - 2))));
  size_t k_max = (size_t) (std::max)((std::min)(static_cast<long>(m_beta.size() - 1), static_cast<long>(std::trunc(std::floor(t + 2)))), (long) 0);
  for (size_t k = k_min; k <= k_max; ++k)
  {
      z += m_beta[k]*b3_spline(t - k);
  }
  return z;
}

mt_real spline::prime(const mt_real x) const
{
    mt_real z = 0;
    mt_real t = m_h_inv*(x - m_a) + 1;

    size_t k_min = (size_t) (std::max)(static_cast<long>(0), static_cast<long>(std::trunc(std::ceil(t - 2))));
    size_t k_max = (size_t) (std::min)(static_cast<long>(m_beta.size() - 1), static_cast<long>(std::trunc(std::floor(t + 2))));

    for (size_t k = k_min; k <= k_max; ++k)
    {
        z += m_beta[k]*b3_spline_prime(t - k);
    }
    return z*m_h_inv;
}

#endif /* SPLINE_HPP */
