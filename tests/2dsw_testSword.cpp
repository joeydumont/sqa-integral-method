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
using namespace sp_bessel;

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
  virtual T evaluateRefractiveIndex(double r, double theta, std::complex<double> k)
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

arma::cx_mat scatteringMatrixDiff(std::complex<double> k,
                                  std::complex<double> n1,
                                  std::complex<double> n2, 
                                  unsigned int mMax)
{
  arma::cx_vec scatMatDiff(2*mMax+1);
  for (unsigned int i = 0; i < 2*mMax+1; i++)
  {
    int m = i - mMax;

    std::complex<double> A,B,Ap,Bp;
    A  =    n1*besselJp(m,n1*k,1)*hankelH2(m,n2*k)-n2*besselJ(m,n1*k)*hankelH2p(m,n2*k,1);
    B  =    n1*besselJp(m,n1*k,1)*hankelH1(m,n2*k)-n2*besselJ(m,n1*k)*hankelH1p(m,n2*k,1);
    Ap = n1*n1*besselJp(m,n1*k,2)*hankelH2(m,n2*k)-n2*n2*besselJ(m,n1*k)*hankelH2p(m,n2*k,2);
    Bp = n1*n1*besselJp(m,n1*k,2)*hankelH1(m,n2*k)-n2*n2*besselJ(m,n1*k)*hankelH1p(m,n2*k,2);

    scatMatDiff(i) = -(Ap*B-A*Bp)/(B*B);
  }

  return arma::diagmat(scatMatDiff);
}

int main(int argc, char* argv[])
{
  HomoCircCavity<double> cav(2.0,1.0,1.0);
  SurfaceMesh<double> meshCav(100,100,cav);
  meshCav.prepareMesh();

  // TM pol.
  BDSword_TM<double> bdsword(1.0, 25, meshCav);

  /*! Computation of the kernel and its roots. */
  arma::cx_mat kernelP = bdsword(1.0);
  std::complex<double> root = rootsMatrix<std::complex<double>, BDSword_TM<double>&>(bdsword, 1.0, 0.1);

  /*! Computation of the scattering matrix and derivative.

  arma::cx_vec test = bdsword.computeInteriorField<func_ptr>(besselJ0);
  arma::abs(test).eval().save("2dsw_testSword.dat", arma::raw_ascii);
  arma::cx_mat scatMat = bdsword.computeScatteringMatrix(1.0);
  //scatMat.print();

  arma::cx_mat scatMatDiff = matrixComplexDerivative<std::complex<double>, BDSword_TM<double>&>
                            (bdsword, 1.0, 0.1);
  arma::real(scatMatDiff).eval().save("scatMatDiffReal.dat", arma::raw_ascii);
  arma::imag(scatMatDiff).eval().save("scatmatDiffImag.dat", arma::raw_ascii);

  arma::mat diff = arma::abs(scatMatDiff-scatteringMatrixDiff(1.0,2.0,1.0,25));
  diff.save("scatMatDifference.dat", arma::raw_ascii);

  std::cout << arma::max(arma::max(arma::abs(scatMatDiff-scatteringMatrixDiff(1.0,2.0,1.0,25)))) << std::endl;

  arma::mat centerPositions = meshCav.getCenterPositions();
  centerPositions.save("2dsw_testSwordCenterPos.dat", arma::raw_ascii);

  // TE pol.
  //BDSword_TE_const<double> bdswordTE(meshCav, 1.0);
  //arma::cx_vec intFieldTE = bdswordTE.computeInteriorField<arma_ptr>(gradBesselJ0);
  //arma::abs(intFieldTE).eval().save("2dsw_testSword_te.dat", arma::raw_ascii);
  */

  return 0;
}
