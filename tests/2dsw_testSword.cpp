/*******************************************************-/
 * This source code is subject to the terms of the GNU  -/
 * Public License. If a copy of the GPL was not         -/
 * distributed with this file, you can obtain one at    -/
 * https://www.gnu.org/licenses/gpl.html.               -/
 ********************************************************/

/** \file 2dsw_testSword.cpp
 *  \author Joey Dumont <joey.dumont@gmail.com>
 *  \since 2014-07-29
 *  \date 2014-07-29 
 *  \brief Tests the TM solver of 2D-SWORD.
 * 
 * We define some cavity geomeries and solve the scattering problem.
 *
 * \copyright GPL
 */

#include <2dsw.h>
#include <armadillo>
#include <complex_bessel.h>

using namespace BD_SWORD;

template <class T>
class HomoCircCavity : public Cavity<T>
{

public:
  HomoCircCavity(T _nIn, T _n0, double _rMax) : Cavity<T>(_n0, _rMax){nIn=_nIn;}
  virtual double boundary(double theta)
  {
    return this->rMax;
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

typedef std::complex<double> (*func_ptr)(double r, double theta);

std::complex<double> besselJ0(double r, double theta)
{
  return sp_bessel::besselJ(0,r)*exp(std::complex<double>(0.0,1.0)*0.0*theta);
}

int main(int argc, char* argv[])
{
  HomoCircCavity<double> cav(2.0,1.0,1.0);
  SurfaceMesh<double> meshCav(75,75,cav);
  meshCav.prepareMesh();
  BDSword_TM<double> bdsword(meshCav, 1.0);
  arma::cx_vec test = bdsword.computeInteriorField<func_ptr>(besselJ0);
  arma::abs(test).eval().save("2dsw_testSword.dat", arma::raw_ascii);

  // We compute the S00 element manually. 
  std::complex<double> sum(0.0,0.0);
  arma::mat centerPositions = meshCav.getCenterPositions();
  centerPositions.save("2dsw_testSwordCenterPos.dat", arma::raw_ascii);
  arma::vec areas = meshCav.getAreaCells();
  for (int i=0; i<centerPositions.n_rows; i++)
  {
    sum+= -3.0
    *areas(i)
    *test(i)
    *besselJ0(sqrt(pow(centerPositions(i,0),2.0)+pow(centerPositions(i,1),2.0)), 0.0);
  }

  std::cout << "S_00: " << 1.0-0.5*std::complex<double>(0.0,1.0)*sum << std::endl;

  return 0;
}
