/*******************************************************-/
 * This source code is subject to the terms of the GNU  -/
 * Public License. If a copy of the GPL was not         -/
 * distributed with this file, you can obtain one at    -/
 * https://www.gnu.org/licenses/gpl.html.               -/
 ********************************************************/

/** \file 2dsw_testCommonFunctions.cpp
 *  \author Joey Dumont <joey.dumont@gmail.com>
 *  \since 2014-07-29
 *  \date 2014-08-04 
 *  \brief Tests the common functions of the library. 
 *  \copyright GPL
 */

#include <2d-sword.h>
#include <armadillo>
#include <iostream>
#include <complex>

using namespace BD_SWORD;

// We test the complex derivative function.
typedef arma::cx_mat (*f_ptr)(std::complex<double>);
arma::cx_mat CosZ(std::complex<double> z)
{
  arma::cx_mat result(2,2);
  result(0,0) = cos(z);
  result(0,1) = z*z;
  result(1,0) = atan(z);
  result(1,1) = exp(z);
  return result;
}

arma::cx_mat OneOverZ(std::complex<double> z)
{
  arma::cx_mat result(1,1);
  result(0,0) = 1.0/z;
  return result;
}

int main(int argc, char* argv[])
{
  arma::vec x = arma::linspace(-datum<double>::pi, datum<double>::pi, 10);
  arma::mat angles(10,10);

  for (int i=0; i<10; i++)
  {
    for (int j=0; j<10; j++)
    {
      angles(j,i) = atan2_pos(x(j),x(i));
    }
    x(i) = user_mod(x(i), 2.0*datum<double>::pi);
  }

  x.print();
  std::cout << std::endl;
  angles.print();

  std::complex<double> z(1.0,1.0);
  arma::cx_mat testCos = matrixComplexDerivative<std::complex<double>, f_ptr>(&CosZ, z, 0.1);
  testCos.print();

  std::complex<double> pole = polesMatrix<std::complex<double>, f_ptr>(&OneOverZ, 2.0);
  std::cout << pole << std::endl;

  return 0;
}