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
 *  \date 2014-07-31
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
public:
  /*! Constructor sets the cavity and the mesh.
   *    @param[in] _mesh Mesh object. */
  BDSword(std::complex<double> _k, unsigned int _Mmax, SurfaceMesh<T>& _mesh) : mesh(_mesh){k=_k;Mmax=_Mmax;}

  /*! Computes the interior field due to an incident field described by func. 
   *    @param[in] func_type Incident field.
   *    @retval Returns the interior field due to the incident field. */
  template <class func_type>
  arma::cx_vec computeInteriorField(func_type func);

  /*! Computes the scattering matrix associated with a cavity. 
   *    @param[in] Mmax Maximum angular momentum to consider.
   *    @retval scatMat Scattering matrix of truncated size (2*Mmax+1)x(2*Mmax+1).*/
  virtual arma::cx_mat computeScatteringMatrix(std::complex<double> k, unsigned int Mmax) = 0;

  /*! @name Accessor Functions */
  ///@{
  SurfaceMesh<T>& getMesh(){return mesh;}
  void setMesh(SurfaceMesh<T>& _mesh){mesh=_mesh;}
  ///@}

protected:
  std::complex<double> k;
  unsigned int Mmax;
  arma::cx_vec interiorField;
  arma::cx_mat scatMat;
  SurfaceMesh<T>& mesh;
  
};// class BDSword
} // namespace BD_SWORD

#endif // BD_SWORD_H