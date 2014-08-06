/*******************************************************-/
 * This source code is subject to the terms of the GNU  -/
 * Public License. If a copy of the GPL was not         -/
 * distributed with this file, you can obtain one at    -/
 * https://www.gnu.org/licenses/gpl.html.               -/
 ********************************************************/

#ifndef BD_SWORD_TE_CONST_H
#define BD_SWORD_TE_CONST_H

/** \file 2dsw_sword_te_const.h
 *  \author Joey Dumont <joey.dumont@gmail.com>
 *  \since 2014-07-31
 *  \date 2014-07-31
 *  \brief Implements the solution method for the TE polarization. 
 *  \bug Doesn't give good answer.
 * 
 * This class defines the necessary functions and numerical 
 * algorithms for the solution of the scattering problem 
 * in the TE polarization. 
 * 
 * \copyright GPL
 */

#include <complex_bessel.h>

#include <2d-sword_bits/2dsw_sword.h>

/// Generic namesapce for the whole library.
namespace BD_SWORD {

template <class T>
class BDSword_TE_const : public BDSword<T>
{
public:
  /*! Constructor sets the cavity, the mesh and the wavenumber.
   *    @param[in] _mesh Mesh object. 
   *    @param[in] _k Wavenumber for which to solve the scattering problem.*/
   BDSword_TE_const(std::complex<double> _k, unsigned int _Mmax, SurfaceMesh<T>& _mesh) 
                    : BDSword<T>(_k, _Mmax, _mesh){}

   template <class func_type>
   arma::cx_vec computeInteriorField(func_type func)
   {
    // Get the positions at which to evaluate the field.
    arma::mat centerPositions = this->mesh.getCenterPositions();

    // Evaluate the incident field at those points.
    unsigned int size = centerPositions.n_rows;
    arma::cx_mat incField(size,3);
    for (unsigned int i = 0; i < size; i++)
    {
      double r = arma::norm(centerPositions.row(i));
      double theta = atan2_pos(centerPositions(i,1),centerPositions(i,0));
      incField.row(i) = func(r,theta);
    }

    // Evaluate the kernel at the necessary points.
    arma::cx_mat kernel = evaluateKernel(centerPositions, this->mesh.getAreaCells());

    // Prepare the linear system to be solved.
    arma::cx_mat kernelP = arma::zeros<arma::cx_mat>(3*size,3*size);
    kernelP(arma::span(0,size-1), arma::span(0,size-1)) = arma::eye<arma::cx_mat>(size,size);
    kernelP(arma::span(size,3*size-1), arma::span(size,3*size-1)) += 
      (1.0-pow(this->mesh.cav.getExtPotential(),2.0)/pow(this->mesh.cav(0.0,0.0),2.0))
      *kernel;

      return arma::solve(kernelP, arma::vectorise(incField));
   } // computeInteriorField()

   arma::cx_mat computeScatteringMatrix(unsigned int Mmax)
   {
    return arma::zeros<arma::cx_mat>(2*Mmax+1,2*Mmax+1);
   }

   double k;

 protected:
  /*! Computes the kernel of the linear equation. 
   *    @param[in] centerPositions Matrix containing the positions at which to evaluate the kernel.
   *    @param[in] areaCells Area of the cells in the mesh. 
   *    @retval evaluateKernel Returns the value of the kernel at these points.*/
  arma::cx_mat evaluateKernel(arma::mat centerPositions, arma::vec areaCells)
  {
    // We determine the size of the kernel.
    unsigned int size = centerPositions.n_rows;
    arma::cx_mat kernel(2*size,2*size);

    // Gradient in x direction.
    for (unsigned int i = 0; i < size; i++)
    {
      double r1 = arma::norm(centerPositions.row(i));
      double phi1 = atan2_pos(centerPositions(i,1), centerPositions(i,0));
      for (unsigned int j = 0; j < size; j++)
      {
        // Compute diagonal elements.
        if (i==j)
        {
          kernel(i,j) = 0.0;
        }

        // Compute off-diagonal elements.
        else
        {
          // Compute the distance.
          double r2 = arma::norm(centerPositions.row(j));
          double phi2 = atan2_pos(centerPositions(j,1), centerPositions(j,0));
          double d = sqrt(r1*r1+r2*r2-2.0*r1*r2*cos(phi1-phi2));

          // Compute the Hankel function.
          kernel(i,j) = -0.25*datum<double>::i
                        *this->mesh.cav.getExtPotential()*k*(r1*cos(phi1)-r2*cos(phi2))/d
                        *sp_bessel::hankelH1(0,this->mesh.cav.getExtPotential()*k*d)
                        *areaCells(j);

          kernel(arma::span(size,2*size-1), arma::span(size,2*size-1))(i,j) =
                        -0.25*datum<double>::i
                        *this->mesh.cav.getExtPotential()*k*(r1*sin(phi1)-r2*sin(phi2))/d
                        *sp_bessel::hankelH1(0,this->mesh.cav.getExtPotential()*k*d)
                        *areaCells(j);
        }
      }
    }
    return kernel;
  } // evaluateKernel()
};// class BDSWord_TE_const
} // namespace BD_SWORD

#endif // BD_SWORD_TE_CONST_H