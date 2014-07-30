/*******************************************************-/
 * This source code is subject to the terms of the GNU  -/
 * Public License. If a copy of the GPL was not         -/
 * distributed with this file, you can obtain one at    -/
 * https://www.gnu.org/licenses/gpl.html.               -/
 ********************************************************/

/** \file 2dsw_testCommonFunctions.cpp
 *  \author Joey Dumont <joey.dumont@gmail.com>
 *  \since 2014-07-29
 *  \date 2014-07-29 
 *  \brief Tests the common functions of the library. 
 *  \copyright GPL
 */

#include <2dsw.h>
#include <armadillo>
#include <iostream>

using namespace BD_SWORD;

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
  return 0;
}