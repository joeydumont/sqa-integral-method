/*******************************************************-/
 * This source code is subject to the terms of the GNU  -/
 * Public License. If a copy of the GPL was not         -/
 * distributed with this file, you can obtain one at    -/
 * https://www.gnu.org/licenses/gpl.html.               -/
 ********************************************************/

#ifndef COMMON_FUNCTIONS_H
#define COMMON_FUNCTIONS_H

/** \file 2dsw_commonFunctions.h
 *  \author Joey Dumont <joey.dumont@gmail.com>
 *  \since 2014-07-24
 *  \date 2014-07-31
 *  \brief Defines some common utility functions for the library.
 *
 * This file defines some functions that are used inside the whole
 * library. They might concern error handling, short numerical algorithms
 * and other types of functions. 
 *
 * \copyright GPL
*/

#include <cmath>
#include <complex>
#include <armadillo>

#include <2d-sword_bits/2dsw_swordConstants.h>

/// Namespace for the whole library.
namespace BD_SWORD {

/*! Proper definition of modular arithmetic.
 *    @param[in] value Value of which we wish to take the modulo.
 *    @param[in] modulo Module of the arithmetic.
 *    @retval user_mod value mod modulo. */
inline double user_mod(double value, double modulo)
{
    return value-modulo*floor(value/modulo);
} // user_mod()

/*! Proper definition of the atan2 function. Returns a value in [0, 2*pi).
 *    @param[in] x x-coordinate.
 *    @param[in] y y-coordinate.
 *    @retval atan2_pos Value of the angle between x and y in [0,2*pi).*/  
inline double atan2_pos(double y, double x)
{
  return user_mod(atan2(y, x), 2.0*datum<double>::pi);
} // atan2_pos()

/*! This function computes the derivative of matrices using
 * the Romberg method (which is based on Richardson's deferred
 * approach to the limit).
 *  @param[in] func Functor or function pointer that evaluates the matrix function.
 *  @param[in] x Value at which we want to evaluate the derivative.
 *  @param[in] h Initial value of the stepsize.
 *  @param[in] b Factor by which we reduce the stepsize each iteration.
 *  @retval matrixDerivative Returns the value of the matrix derivative. */
template <class mat_type, class func_type>
arma::Mat<mat_type> matrixComplexDerivative(func_type func, const std::complex<double> x, double h, double b=1.4, double tol = 1.0e-5)
{
  // Compute the square of the stepsize reducing factor.
  double b2 = b*b;
  double fac;

  // We prepare the error variables.
  const double safe = 2.0;
  double err = std::numeric_limits<double>::max();
  double errt;

  // We prepare a cube to store the matrices.
  std::complex<double> I(0.0,1.0);
  arma::Mat<mat_type> ans = (func(x+h)-func(x-h)-I*(func(x+I*h)-func(x-I*h)))/(4.0*h);
  arma::Cube<mat_type> tableauOld(ans.n_rows, ans.n_cols, 5);
  arma::Cube<mat_type> tableauNew(ans.n_rows, ans.n_cols, 6);

  // We store the first function evaluation in the first element of the
  // tableau.
  tableauOld.slice(0) = ans;

  // At each iteration, we compute the first element in the ith row.
  // We limit the size of the tableau to 10 to limit the number of function
  // evaluations.
  unsigned int i = 1;
  while(i<10)
  {
    // We resize the cubes only if necessary.
    if (i>5)
    {
      tableauOld.resize(ans.n_rows, ans.n_cols, i+1);
      tableauNew.set_size(ans.n_rows, ans.n_cols, i+2);
    }

    // Reduce the stepsize and evaluate the function.
    h /= b;
    tableauNew.slice(0) = (func(x+h)-func(x-h)-I*(func(x+I*h)-func(x-I*h)))/(4.0*h);
    fac = b2;

    for (int j=1;j<=i;j++)
    {
      // Compute extrapolations of higher order.
      tableauNew.slice(j) = (tableauNew.slice(j-1)*fac-tableauOld.slice(j-1))/(fac-1.0);
      fac = b2*fac;

      // We verify the error made.
      errt = std::max(arma::norm(tableauNew.slice(j)-tableauNew.slice(j-1)),
                    arma::norm(tableauNew.slice(j)-tableauOld.slice(j-1)));

      // If the error has decreased, we accept the answer.
      if (errt <= err)
      {
        err=errt;
        ans = tableauNew.slice(j);
        }
      }

      // If higher order is worse by a significant factor, we quit.
      if (arma::norm(tableauNew.slice(i)-tableauOld.slice(i-1)) >= safe*err
        ||err <= tol) break;


      // Bookkeeping for next iteration.
      i++;
      tableauOld = tableauNew;
    }

    return ans;
} // matrixComplexDerivative()

/*! This function finds the poles of a matrix. We use Newton's method 
 * to compute the roots of the determinant function.
 *    @param[in] 
 *    @param[in]
 *    @param[in]
 *    @retval poleMatrix A pole of the matrix. */
template <class mat_type, class func_type>
std::complex<double> polesMatrix(func_type matrix, std::complex<double> pole, const double h = 0.1, const double tol = 1.0e-5)
{
  // We compute the error
  double errNew = std::numeric_limits<double>::max();
  double errOld = std::numeric_limits<double>::max();
  double safe = 2.0;
  std::complex<double> poleNew = pole;

  while (errOld>=tol)
  {
    // We compute the matrix inverse and derivative.
    arma::Mat<mat_type> i_mat = arma::inv(matrix(pole));
    arma::Mat<mat_type> d_mat = matrixComplexDerivative<mat_type, func_type>(matrix, pole, h);

    // We compute the new pole.
    poleNew = pole + 1.0/arma::trace(d_mat*i_mat);
    errNew = fabs(poleNew-pole);

    if (errNew > safe*errOld)
    {
      break;
    }

    // Bookkeeping
    errOld = errNew;
    pole = poleNew;
    std::cout << pole << std::endl;
  }

  return pole;
}

} // namespace BD_SWORD

#endif // COMMON_FUNCTIONS_H