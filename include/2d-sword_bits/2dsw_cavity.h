/*******************************************************-/
 * This source code is subject to the terms of the GNU  -/
 * Public License. If a copy of the GPL was not         -/
 * distributed with this file, you can obtain one at    -/
 * https://www.gnu.org/licenses/gpl.html.               -/
 ********************************************************/ 

#ifndef CAVITY_H
#define CAVITY_H

/** \file 2dsw_cavity.h
 *  \author Joey Dumont <joey.dumont@gmail.com>
 *  \since 2014-07-24
 *  \date 2014-07-24
 *  \brief Defines the refractive index function and geometry of the cavity.
 *
 * This class defines both the geometry and the refractive
 * index profile of the cavity. The refractive index can be given 
 * in either polar or Cartesian coordinates and be either real
 * or complex. The class also supports frequency-dependent 
 * refractive indices.
 *
 * \copyright GPL
 */

#include <complex>
#include <string>
#include <iostream>
#include <cmath>

#include <2d-sword_bits/2dsw_commonFunctions.h>

/// Generic namesapce for the whole library.
namespace BD_SWORD {

template <class T>
class Cavity
{
public:
  /// Constructor sets the radius of the LSS.
  Cavity(T _n0, double _rMax) {setExtPotential(_n0);setLSSRadius(_rMax);}
  
  /*! Refractive index sampling. 
   *    @param[in] x1 Value of the first coordinate. 
   *    @param[in] x2 Value of the second coordinate.
   *    @param[in] k Frequency at which to evaluate the refractive index.
   *    @param[in] coord_sys Coordinate system. 
   *    @retval refractiveIndex Real or complex refractive index at the point (x1,y1). */
  T operator()(double x1, double x2, double k = 0.0, std::string coord_sys = "polar")
  {
    T result;
    if (coord_sys == std::string("polar"))
      result = evaluateRefractiveIndex(x1, x2, k);

    else if (coord_sys == std::string("cartesian"))
    {
      double r = sqrt(x1*x1+x2*x2);
      double theta = atan2_pos(x2,x1);
      result = evaluateRefractiveIndex(r, theta, k);
     }

    else
    {
      std::cout << "The only available options are \"polar\" and \"cartesian\"" << std::endl;
      /// \todo Better error handling. 
      //throw;
    }

    return result;
  }

  /// Boundary of the cavity in polar coordinates, i.e. r(theta).
  virtual double boundary(double theta) = 0;

  /// Writes down information about the cavity in a file.
  virtual void RefIndexInfo() = 0;

  /*! @name Accessor Functions */
  ///@{
  void setLSSRadius(double _rMax){rMax=_rMax;}
  double getLSSRadius(){return rMax;}

  void setExtPotential(T _n0){n0=_n0;}
  T getExtPotential(){return n0;}
  ///@}

protected:
  double rMax;
  T n0;

  /*! Private virtual function to be implemented in derived classes. */
  virtual T evaluateRefractiveIndex(double r, double theta, double k) = 0;

};// class Cavity
} // namespace BD_SWORD

#endif // CAVITY_H