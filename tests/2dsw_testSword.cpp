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

#include <2d-sword.h>
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
typedef arma::cx_rowvec (*arma_ptr)(double r, double theta);

std::complex<double> besselJ0(double r, double theta)
{
  return sp_bessel::besselJ(0,r)*exp(std::complex<double>(0.0,1.0)*0.0*theta);
}

arma::cx_rowvec gradBesselJ0(double r, double theta)
{
  arma::cx_rowvec grad(3);
  grad(0) = sp_bessel::besselJ(0,r);
  grad(1) = cos(theta)*sp_bessel::besselJ(1,r);
  grad(2) = -sin(theta)*sp_bessel::besselJ(1,r);
  return grad;
}

int main(int argc, char* argv[])
{
  HomoCircCavity<double> cav(2.0,1.0,1.0);
  SurfaceMesh<double> meshCav(25,25,cav);
  meshCav.prepareMesh();

  // TM pol.
  BDSword_TM<double> bdsword(meshCav, 1.0);
  arma::cx_vec test = bdsword.computeInteriorField<func_ptr>(besselJ0);
  arma::abs(test).eval().save("2dsw_testSword.dat", arma::raw_ascii);
  arma::cx_mat scatMat = bdsword.computeScatteringMatrix(0);
  scatMat.print();

  arma::mat centerPositions = meshCav.getCenterPositions();
  centerPositions.save("2dsw_testSwordCenterPos.dat", arma::raw_ascii);

  // TE pol.
  BDSword_TE_const<double> bdswordTE(meshCav, 1.0);
  arma::cx_vec intFieldTE = bdswordTE.computeInteriorField<arma_ptr>(gradBesselJ0);
  arma::abs(intFieldTE).eval().save("2dsw_testSword_te.dat", arma::raw_ascii);
  return 0;
}
