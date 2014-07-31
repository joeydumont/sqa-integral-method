/*******************************************************-/
 * This source code is subject to the terms of the GNU  -/
 * Public License. If a copy of the GPL was not         -/
 * distributed with this file, you can obtain one at    -/
 * https://www.gnu.org/licenses/gpl.html.               -/
 ********************************************************/

/** \file 2dsw_testSurfaceMesh.cpp
 *  \author Joey Dumont <joey.dumont@gmail.com>
 *  \since 2014-07-29
 *  \date 2014-07-29 
 *  \brief Tests the meshing capabilities of the library. 
 * 
 * We define some cavity geomeries and inspect the meshes.
 *
 * \copyright GPL
 */

#include <2d-sword.h>
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
  SurfaceMesh<double> meshCav(100,100,cav);
  meshCav.prepareMesh();
  meshCav.getCenterPositions().save("2dsw_testSurfaceMesh_centerPositions.dat", arma::raw_ascii);
  return 0;
}
