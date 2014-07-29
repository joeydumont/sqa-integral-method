/*******************************************************-/
 * This source code is subject to the terms of the GNU  -/
 * Public License. If a copy of the GPL was not         -/
 * distributed with this file, you can obtain one at    -/
 * https://www.gnu.org/licenses/gpl.html.               -/
 ********************************************************/

#ifndef BD_SWORD_H
#define BD_SWORD_H

/** \file 2dsw_sword.h
 *  \author Joey Dumont <joey.dumont@gmail.com>
 *  \since 2014-07-28
 *  \date 2014-07-28
 *  \brief Solves the scattering problem of a given cavity. 
 * 
 * This class implements diverse solution methods for the 
 * scattering problem in bidimensional cavities.
 * 
 * \copyright GPL
 */

#include <2d-sword_bits/2dsw_cavity.h>
#include <2d-sword_bits/2dsw_surfaceMesh.h>

/// Generic namespace for the whole library
namespace BD_SWORD {

template <class T>
class BDSword
{
  /*! Constructor sets the cavity and the mesh.
   *    @param[in] */
  BDSword(Cavity<T>& cav, SurfaceMesh& mesh)
  {

  }

public:
  /*! Computes the interior field due to an incident field described by func. 
   *    @param[in] func_type Incident field.
   *    @retval Returns the interior field due to the incident field. */
  template <class func_type>
  arma::cx_vec computeInteriorField(func_type func);

  /*! Computes the scattering matrix of the cavity. */
  arma::cx_mat computeScatteringMatrix();

protected:
  arma::cx_vec interiorField;
  arma::cx_mat scatMat;
  
};// class BDSword
} // namespace BD_SWORD

#endif // BD_SWORD_H