/*******************************************************-/
 * This source code is subject to the terms of the GNU  -/
 * Public License. If a copy of the GPL was not         -/
 * distributed with this file, you can obtain one at    -/
 * https://www.gnu.org/licenses/gpl.html.               -/
 ********************************************************/

#ifndef SWORD_CONSTANTS_H
#define SWORD_CONSTANTS_H

/** \file 2dsw_swordConstants.h
 *  \author Joey Dumont <joey.dumont@gmail.com>
 *  \since 2014-07-25
 *  \date 2014-07-25
 *  \brief Defines some static constants.
 *
 * This file defines some mathemetical constants that are of use
 * in the library. 
 *
 * \copyright GPL
 */

#include <complex>

/// Namespace for the whole library.
namespace BD_SWORD {

template <typename eT = double>
class datum
{
public:
  static const eT pi;
  static const std::complex<eT> i;
};

template <typename eT> const eT datum<eT>::pi = eT(3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679);
template <typename eT> const std::complex<eT> datum<eT>::i = std::complex<eT>(0.0,1.0);
} // namespace 2D-SWORD

#endif // SWORD_CONSTANTS_H