/*******************************************************-/
 * This source code is subject to the terms of the GNU  -/
 * Public License. If a copy of the GPL was not         -/
 * distributed with this file, you can obtain one at    -/
 * https://www.gnu.org/licenses/gpl.html.               -/
 ********************************************************/

/** \file 2dsw_testSwordTimingTM.cpp
 *  \author Joey Dumont <joey.dumont@gmail.com>
 *  \since 2014-07-31
 *  \date 2014-07-31
 *  \brief Tests the TM solver of 2D-SWORD. 
 * 
 * We define some cavity geomeries and solve the scattering problem.
 * We determine the accuracy of the method as a function of the area
 * of the cells. 
 *
 * We also determine the takes it takes to compute a fully converged
 * (up to a given tolerance) as a function of the angular momentum 
 * channels. 
 *
 * \copyright GPL
 */

#include <2d-sword.h>
#include <armadillo>
#include <complex_bessel.h>
#include <fstream>
#include <time.h>

using namespace BD_SWORD;

typedef arma::cx_mat (*testMatrix)(std::complex<double>);
arma::cx_mat testMatrixPole(std::complex<double> z)
{
  arma::cx_mat test(3,3);
  test(0,0) = 1.;
  test(0,1) = 4.;
  test(0,2) = 5.;
  test(1,0) = 7.;
  test(1,1) = -1.0/z; 
  test(1,2) = 9.;
  test(2,0) = 1.;
  test(2,1) = 4.;
  test(2,2) = 2.;

  return test;
}

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

// Analytical form of the S-matrix.
arma::cx_mat analSmatrix(int M, std::complex<double> nIn, std::complex<double> k, double rad)
{
  // Define a matrix of zeros.
  arma::cx_mat scatMat = arma::zeros<arma::cx_mat>(2*M+1,2*M+1);

  // Loop to fill the diagonal.
  for (unsigned int i=0;i<2*M+1;i++)
  {
    // Preparation of order and variables.
    int order = -M+i;
    std::complex<double> Zc = nIn*k*rad;
    std::complex<double> Zo = k*rad;

    // Compute numerator and denominator.
    std::complex<double> num = nIn*sp_bessel::hankelH2(order,Zo)*sp_bessel::besselJp(order,Zc)
                -sp_bessel::hankelH2p(order,Zo)*sp_bessel::besselJ(order,Zc);
    std::complex<double> den = nIn*sp_bessel::hankelH1(order,Zo)*sp_bessel::besselJp(order,Zc)
                -sp_bessel::hankelH1p(order,Zo)*sp_bessel::besselJ(order,Zc);

    // Compute diagonal entries.
    scatMat(i,i) = -num/den;
  }

  return scatMat;
}

int main(int argc, char* argv[])
{
  int Nmesh = atoi(argv[1]);

  // We prepare the cavity and the mesh. 
  HomoCircCavity<double> cav(2.0,1.0,1.0);
  SurfaceMesh<double> meshCav(Nmesh,Nmesh,cav);
  meshCav.prepareMesh();

  BDSword_TM<double> bdsword(1.0, 25, meshCav);
  clock_t start = clock();
  arma::cx_mat scatMat = bdsword.computeScatteringMatrix(0.5);
  clock_t end = clock();
  arma::cx_mat analSMat = analSmatrix(25, 
                                      std::complex<double>(2.0,0.0),
                                      std::complex<double>(1.0,0.0),
                                      0.5);

  double err = arma::max(arma::max(arma::abs(scatMat-analSMat)));
  std::cout << err << std::endl;

 // std::cout << polesMatrix<std::complex<double>, BDSword_TM<double>&>(bdsword, 0.4);

  //std::cout << polesMatrix<std::complex<double>, testMatrix>(testMatrixPole, 1.0) << std::endl; 
  // We write this in a file. 
  std::ofstream results;
  results.open("convergence-avx2.dat", std::ios::app);
  results.setf(std::ios::scientific);
  results.precision(15);
  results.put(' ');
  results.width(20);
  results << arma::mean(meshCav.getAreaCells());
  results.put(' ');
  results.width(20);
  results << err;
  results.put(' ');
  results.width(20);
  results << ((double)end-(double)start) / CLOCKS_PER_SEC;
  results.put('\n');
  results.close();

  return 0;
}
