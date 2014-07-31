/*******************************************************-/
 * This source code is subject to the terms of the GNU  -/
 * Public License. If a copy of the GPL was not         -/
 * distributed with this file, you can obtain one at    -/
 * https://www.gnu.org/licenses/gpl.html.               -/
 ********************************************************/

#ifndef COMMON_FUNCTIONS_H
#define COMMON_FUNCTIONS_H

/** \file 2dsw_commonFunctions.h
 *  \author Joey Dumont <joey.dumont@gmail.com>
 *  \since 2014-07-24
 *  \date 2014-07-31
 *  \brief Defines some common utility functions for the library.
 *
 * This file defines some functions that are used inside the whole
 * library. They might concern error handling, short numerical algorithms
 * and other types of functions. 
 *
 * \copyright GPL
*/

#include <cmath>
#include <complex>

#include <2d-sword_bits/2dsw_swordConstants.h>

/// Namespace for the whole library.
namespace BD_SWORD {

/*! Proper definition of modular arithmetic.
 *    @param[in] value Value of which we wish to take the modulo.
 *    @param[in] modulo Module of the arithmetic.
 *    @retval user_mod value mod modulo. */
inline double user_mod(double value, double modulo)
{
    return value-modulo*floor(value/modulo);
}

/*! Proper definition of the atan2 function. Returns a value in [0, 2*pi).
 *    @param[in] x x-coordinate.
 *    @param[in] y y-coordinate.
 *    @retval atan2_pos Value of the angle between x and y in [0,2*pi).*/  
inline double atan2_pos(double y, double x)
{
  return user_mod(atan2(y, x), 2.0*datum<double>::pi);
}


} // namespace BD_SWORD

#endif // COMMON_FUNCTIONS_H