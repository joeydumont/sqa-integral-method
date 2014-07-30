/*******************************************************-/
 * This source code is subject to the terms of the GNU  -/
 * Public License. If a copy of the GPL was not         -/
 * distributed with this file, you can obtain one at    -/
 * https://www.gnu.org/licenses/gpl.html.               -/
 ********************************************************/

/** \file 2dsw_testCavity.cpp
 *  \author Joey Dumont <joey.dumont@gmail.com>
 *  \since 2014-07-29
 *  \date 2014-07-29 
 *  \brief Tests the functions of the Cavity abstract class.
 * 
 * We define some cavity geometries and refractive index profile 
 * and test that the class behaves properly.
 *
 * \copyright GPL
 */

#include <2dsw.h>
#include <armadillo>

using namespace BD_SWORD;

template <class T>
class HomoCircCavity : public Cavity<T>
{

public:
  HomoCircCavity(T _nIn, T _n0, double _rMax) : Cavity<T>(_n0, _rMax){nIn=_nIn;}
  virtual double boundary(double theta)
  {
    return 1.0;
  }

  void RefIndexInfo(){}

protected:
  T nIn;

private:
  virtual T evaluateRefractiveIndex(double r, double theta, double k)
  {
    if (r<=this->rMax)
      return nIn;
    else
      return this->n0;
  }
};

int main(int argc, char* argv[])
{
  HomoCircCavity<double> cav(2.0,1.0,1.0);
  arma::vec x = arma::linspace(-1.0,1.0,200);
  arma::mat refIndex(200,200);

  for (int i=0; i<200; i++)
  {
    for (int j=0; j<200; j++)
    {
      double r = sqrt(x(i)*x(i)+x(j)*x(j));
      double theta = BD_SWORD::atan2_pos(x(j),x(i));
      refIndex(i,j) = cav(r,theta);
    }
  }

  refIndex.save("2dsw_testCavity.dat", arma::raw_ascii); 
  return 0;
}