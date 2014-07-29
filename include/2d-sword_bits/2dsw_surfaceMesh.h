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
#include <vector>
#include <algorithm>
#include <armadillo>

#include <2d-sword_bits/2dsw_cavity.h>

/// Generic namespace for the whole library
namespace BD_SWORD {

template <class T>
class SurfaceMesh
{
  /*! Constructor sets the mesh size and the cavity. 
   *    @param[in] _Nx Number of cells in the x direction.
   *    @param[in] _Ny Number of cells in the y direction.
   *    @param[in] _cav Reference to a Cavity object.*/
  SurfaceMesh(int _Nx, int _Ny, Cavity<T>& _cav)
  {
    setNumberMeshPointsX(_Nx);
    setNumberMeshPointsY(_Ny);
    setCavity(_cav);
  }

public:
  /*! Prepares the mesh used to discretize the cavity. */
  void prepareMesh()
  {
    {
  // We determine the size of the mesh.
  double rMax = cav.getLSSRadius();

  // We generate the boundary points in the x and y directions. 
  arma::vec xPoints = arma::linspace(-rMax, rMax, Nx+1);
  arma::vec yPoints = arma::linspace(-rMax, rMax, Ny+1);
  double dx = xPoints(1)-xPoints(0);
  double dy = yPoints(1)-yPoints(0);

  // We compute and store the values of the associated
  // center points.
  centerPoints.set_size(Nx,Ny);
  for (unsigned int idx = 0; idx < Nx; idx++)
  {
    for (unsigned int idy = 0; idy < Ny; idy++)
    {
      centerPoints(idx,idy).set_size(2);
      centerPoints(idx,idy)(0) = (idx+0.5)*dx;
      centerPoints(idx,idy)(1) = (idy+0.5)*dy;
    }
  }

  // We compute the areas of each cell.
  /// \todo Determine boundary cells and adjust area according to boundary curve.
  areaCells.set_size(Nx,Ny);
  areaCells.fill(dx*dy);

  // We determine if the cell are inside the cavity or not.
  interiorBool.resize(Nx);
  for (unsigned int idx = 0; idx < Nx; idx++)
  {
    interiorBool[idx].resize(Ny);
    for (unsigned int idy = 0; idy < Ny; idy++)
    {
      double x = centerPoints(idx,idy)(0);
      double y = centerPoints(idx,idy)(1);
      double r = sqrt(x*x+y*y);
      double theta = atan2_pos(y,x);
      interiorBool[idx][idy] = ((cav.boundary(theta)-r <= 0.0) ? false : true);
    }
  }

  // Set the mesh to ready.
  meshReady = true;
} // SurfaceMesh::prepareMesh()
  }

protected:
  /*! @name Construction Variables 
   * Variables necessary to construct the SurfaceMesh. */
  ///@{
  unsigned int Nx, Ny;
  Cavity<T>& cav;
  ///}

  /*! @name Mesh Variables 
   * Properties of the constructed mesh. */
  ///@{
  arma::field<arma::vec> centerPoints;
  arma::mat areaCells;
  std::vector<std::vector<bool> > interiorBool;
  ///@}

  /*! Bookkeeping Variables */
  ///@{
  bool meshReady = false;
  ///@}

  /*! Accessor Functions */
  ///@{
  unsigned int getNumberMeshPointsX(){return Nx;}
  void setNumberMeshPointsX(int _Nx){Nx=_Nx;setMeshReady(false);}

  unsigned int getNumberMeshPointsY(){return Ny;}
  void setNumberMeshPointsY(int _Ny){Ny=_Ny;setMeshReady(false);}

  Cavity<T>& getCavity(){return cav;}
  void setCavity(Cavity<T>& _cav){cav=_cav;setMeshReady(false);}

  arma::mat getCenterPositions()
  {
  // We return proper values if the mesh is ready.
  if (!meshReady)
  {
    std::cout << "Mesh is not ready. May yield erreneous results." << std::endl;
    throw;
  }
    
  // We vectorize the centerPositions field into a matrix.
  unsigned int size = 0;
  for (unsigned int i = 0; i < Nx; i++)
  {
    size += std::count(interiorBool[i].begin(), interiorBool[i].end(), true);
  }

  arma::mat centerPositions(size,2);
  unsigned int count = 0;
  for (unsigned int i = 0; i < Nx; i++)
  {
    for (unsigned int j = 0; j < Ny; j++)
    {
      if (interiorBool[i][j] == true)
      {
        centerPositions(count, 0) = centerPoints(i,j)(0);
        centerPositions(count, 1) = centerPoints(i,j)(1);
        count++;
        }
    }
  }

  return centerPositions;
} // SurfaceMesh::getCenterPositions()

  arma::vec getAreaCells()
  {
  if (!meshReady)
  {
    std::cout << "Mesh is not ready. May yield erreneous results." << std::endl;
    throw;
  }

  // We vectorize the areaCells matrix into a vector.
  unsigned int size = 0;
  for (unsigned int i = 0; i < Nx; i++)
  {
    size += std::count(interiorBool[i].begin(), interiorBool[i].end(), true);
  }

  arma::vec areaCellsVec(size);
  unsigned int count = 0;
  for (unsigned int i = 0; i < Nx; i++)
  {
    for (unsigned int j = 0; j < Ny; j++)
    {
      if (interiorBool[i][j] == true)
      {
        areaCellsVec(count) = areaCells(i,j);
        count++;
      }
    }
  }

  return areaCellsVec;
} // SurfaceMesh::getAreaCells()

  bool getMeshReady(){return meshReady;}
  ///@}

};// class SurfaceMesh
} // namespace BD_SWORD

#endif // SURFACE_MESH_H