/*******************************************************-/
 * This source code is subject to the terms of the GNU  -/
 * Public License. If a copy of the GPL was not         -/
 * distributed with this file, you can obtain one at    -/
 * https://www.gnu.org/licenses/gpl.html.               -/
 ********************************************************/

#ifndef SURFACE_MESH_H
#define SURFACE_MESH_H

/** \file 2dsw_surfaceMesh.h
 *  \author Joey Dumont <joey.dumont@gmail.com>
 *  \since 2014-07-25
 *  \date 2014-07-25
 *  \brief Computes the mesh of a given cavity. 
 * 
 * This class computes the mesh of a given cavity object. 
 * Given a Cavity object, we generate a square grid consisting
 * of N_x*N_y cells. The center of each cell is stored in a matrix.
 * We then determine if the center of each cells in is in the interior
 * of the cavity. We store this information in a matrix containing boolean 
 * values. We also store the area of each cells in a third matrix. 
 * We then form a vector with the values in the first matrix
 * that have a corresponding TRUE value in the second matrix. 
 *
 * \copyright GPL
 */

#include <complex>
#include <cmath>

#include <2d-sword_bits/2dsw_cavity.h>

/// Generic namespace for the whole library
namespace 2D-SWORD {

class SurfaceMesh
{
  /*! Constructor sets the mesh size and the cavity. 
   *    @param[in] _Nx Number of cells in the x direction.
   *    @param[in] _Ny Number of cells in the y direction.
   *    @param[in] _cav Reference to a Cavity object.*/
  SurfaceMesh(int _Nx, int _Ny, Cavity& _cav);

public:
  /*! Prepares the mesh used to discretize the cavity. */
  void prepareMesh();

protected:
  /*! @name Construction Variables 
   * Variables necessary to construct the SurfaceMesh. */
  ///@{
  int Nx, Ny;
  Cavity& cav;
  ///}

  /*! @name Mesh Variables 
   * Properties of the constructed mesh. */
  ///@{
  arma::mat centerPositions;
  arma::mat areaCells;
  std::vector<std::vector<bool> > interiorBool;
  ///@}

  /*! Bookkeeping Variables */
  ///@{
  bool meshReady(false);
  ///@}

};// class SurfaceMesh
} // namespace 2D-SWORD

