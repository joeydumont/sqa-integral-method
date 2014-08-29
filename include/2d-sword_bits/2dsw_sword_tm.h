/*******************************************************-/
 * This source code is subject to the terms of the GNU  -/
 * Public License. If a copy of the GPL was not         -/
 * distributed with this file, you can obtain one at    -/
 * https://www.gnu.org/licenses/gpl.html.               -/
 ********************************************************/

#ifndef BD_SWORD_TM_H
#define BD_SWORD_TM_H

/** \file 2dsw_sword_tm.h
 *  \author Joey Dumont <joey.dumont@gmail.com>
 *  \since 2014-07-29
 *  \date 2014-07-31
 *  \brief Implements the solution method for the TM polarization. 
 * 
 * This class defines the necessary functions and numerical 
 * algorithms for the solution of the scattering problem 
 * in the TM polarization. 
 * 
 * \copyright GPL
 */

#include <complex_bessel.h>

#include <2d-sword_bits/2dsw_sword.h>

/// Generic namespace for the whole library. 
namespace BD_SWORD {

template <class T>
class BDSword_TM : public BDSword<T>
{
public:
  /*! Constructor sets the cavity and the mesh. 
   *    @param[in] _mesh Mesh object. */
  BDSword_TM(std::complex<double> _k, unsigned int _Mmax, SurfaceMesh<T>& _mesh) 
              : BDSword<T>(_k, _Mmax, _mesh){}

  template <class func_type>
  arma::cx_vec computeInteriorField(func_type func)
  {
    // Get the positions at which to evaluate the field.
    arma::mat centerPositions = this->mesh.getCenterPositions();

    // Evaluate the incident field at those points. 
    unsigned int size = centerPositions.n_rows;
    arma::cx_vec incField(size);
    for (unsigned int i = 0; i < size; i++)
    {
      double r = sqrt(pow(centerPositions(i,0),2.0)+pow(centerPositions(i,1),2.0));
      double theta = atan2_pos(centerPositions(i,1),centerPositions(i,0));
      incField(i) = func(r,theta);
    }

    // Evaluate the kernel at the necessary points.
    arma::cx_mat kernel = evaluateKernel(centerPositions, this->mesh.getAreaCells());

    // Solve the linear equation.
    arma::cx_mat kernelP = arma::eye<arma::cx_mat>(size,size)+kernel;
    return arma::solve(kernelP, incField);
  } // computeInteriorField()

  arma::cx_mat operator()(std::complex<double> _k)
  {
    return getTransferMatrix(_k);
  }

  arma::cx_mat computeScatteringMatrix(std::complex<double> _k, unsigned int _Mmax=25)
  {
    // We prepare the variables.  
    this->k=_k;
    this->Mmax = _Mmax;
    this->scatMat.set_size(2*this->Mmax+1,2*this->Mmax+1);

    // We compute the scattering matrix column-by-column.
    eigenFunctions* eig = new eigenFunctions(-this->Mmax);
    for (unsigned int i = 0; i < 2*this->Mmax+1; i++)
    {
      // We set the "incoming" angular momentum. 
      int m = i-this->Mmax;
      eig->sgnExp = 1;
      eig->M = m;
      arma::cx_vec intField =  computeInteriorField<eigenFunctions>(*eig);

      // We compute the S_m'm, the transition probability from momentum m' to m.
      for (unsigned int j = 0; j < 2*this->Mmax+1; j++)
      {
        int mp = j-this->Mmax;
        eig->sgnExp = -1;
        eig->M = mp;
        for (unsigned int l = 0; l < this->mesh.getAreaCells().n_rows; l++)
        {
          double r = sqrt(pow(this->mesh.getCenterPositions()(l,0), 2.0)
                          +pow(this->mesh.getCenterPositions()(l,1),2.0));
          double theta = atan2_pos(this->mesh.getCenterPositions()(l,1),
                                    this->mesh.getCenterPositions()(l,0));

          this->scatMat(j,i) += eig->operator()(r,theta)
                    *(pow(this->mesh.cav(r,theta),2.0)-pow(this->mesh.cav.getExtPotential(),2.0))
                    *intField(l)
                    *this->mesh.getAreaCells()(l);
        }
        this->scatMat(j,i) *= 0.5*datum<double>::i*pow(this->k,2.0);
        this->scatMat(j,i) += ((mp==m) ? 1.0 : 0.0);
      }
    }
    delete eig;
    return this->scatMat;
  } // computeScatteringMatrix()

  arma::cx_mat getTransferMatrix(std::complex<double> _k = 0.0)
  {
    // Verify default value.
    if (_k != 0.0)
    {
      this->k=_k;
    }

    // We evaluate the kernel of the integral equation.
    arma::cx_mat kernel = evaluateKernel(this->mesh.getCenterPositions(),
                                         this->mesh.getAreaCells());

    return arma::eye<arma::cx_mat>(this->mesh.getCenterPositions().n_rows,this->mesh.getCenterPositions().n_rows)+kernel;
  } // getTransferMatrix()

protected:
  /*! Computes the kernel of the linear equation.
   *    @param[in] centerPositions Matrix containing the positions at which to evaluate the kernel.
   *    @param[in] areaCells Area of the cells in the mesh. 
   *    @retval evaluateKernel Returns the value of the kernel at these points. */
  arma::cx_mat evaluateKernel(arma::mat centerPositions, arma::vec areaCells)
  {
    unsigned int size = centerPositions.n_rows;
    arma::cx_mat kernel(size,size);
    for (unsigned int i = 0; i < size; i++)
    {
      double r1 = sqrt(pow(centerPositions(i,0),2.0)+pow(centerPositions(i,1),2.0));
      double phi1 = atan2_pos(centerPositions(i,1),centerPositions(i,0));
      for (unsigned int j = 0; j < size; j++)
      {
        // Compute diagonal elements.
        if (i==j)
        {
          kernel(i,j) = 0.0;//0.5*(pow(this->mesh.cav(r1,phi1,k),2.0)-pow(this->mesh.cav.getExtPotential(),2.0))
                          ///pow(this->mesh.cav.getExtPotential(),2.0);
        }
        // Compute off-diagonal elements.
        if (i!=j)
        {
          // Compute the distance
          double r2 = sqrt(pow(centerPositions(j,0),2.0)+pow(centerPositions(j,1),2.0));
          double phi2 = atan2_pos(centerPositions(j,1),centerPositions(j,0));
          double d = sqrt(r1*r1+r2*r2-2.0*r1*r2*cos(phi1-phi2));

          // Compute the Bessel function.
          kernel(i,j) = -0.25*pow(this->k,2.0)*datum<double>::i
                        *sp_bessel::hankelH1(0,this->k*this->mesh.cav.getExtPotential()*d)
                        *(pow(this->mesh.cav(r2,phi2,this->k),2.0)-pow(this->mesh.cav.getExtPotential(),2.0))
                        *areaCells(j);
        }
      }
    }
    return kernel;
  } // evaluateKernel()

  struct eigenFunctions
  {
    eigenFunctions(int _M=0, int _sgnExp=1){M=_M;sgnExp=_sgnExp;}
    int sgnExp;
    int M;
    std::complex<double> operator()(double r, double theta)
    {
      return sp_bessel::besselJ(M,r)*exp((double)sgnExp*datum<double>::i*(double)M*theta);
    }
  };

};// class BDSword_TM
} // namespace BD_SWORD

#endif // BD_SWORD_TM_H
