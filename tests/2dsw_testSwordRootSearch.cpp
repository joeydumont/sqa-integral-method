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
 *  \brief Computes some roots of the homogeneous disk. 
 * 
 * We compute some roots of the homogeneous disk and see if they 
 * correspond to the ones computed via the analytical scattering matrix.
 *
 * \copyright GPL
 */

#include <2d-sword.h>
#include <armadillo>
#include <complex_bessel.h>
#include <fstream>
#include <time.h>

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

  // We prepare the file.
  std::ofstream results;
  results.open("polesHomoCirc.dat", std::ios::app);
  results.setf(std::ios::scientific);
  results.precision(15);
  
  // We compute a fixed set of roots.
  int mMax = 2;
  int jMax = 5;
  std::complex<double> pole;

  for (int i = 0; i < mMax; i++)
  {
    for (int j = 0; j < jMax; j++)
    {
      // We compute the pole.
      clock_t start = clock();
      double initReK = datum<double>::pi/2.0*(i/2.0+j+0.25);
      double initImK = log(abs(1.0)/(3.0))/(2.0*2.0);
      pole 
            = polesMatrix<std::complex<double>, BDSword_TM<double>&>
              (bdsword, std::complex<double>(initReK, initImK));
      clock_t end = clock();

      // We write this in a file. 
      results.put(' ');
      results.width(20);
      results << std::real(pole);
      results.put(' ');
      results.width(20);
      results << std::imag(pole);
      results.put(' ');
      results.width(20);
      results << ((double)end-(double)start) / CLOCKS_PER_SEC;
      results.put('\n');
    }
  }
  
  results.close();

  return 0;
}
