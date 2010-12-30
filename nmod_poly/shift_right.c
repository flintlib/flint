/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2010 William Hart

******************************************************************************/

#include <mpir.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

void nmod_poly_shift_right(nmod_poly_t res, nmod_poly_t poly, long k)
{
   const long len = poly->length - k;

   if (k >= poly->length) /* shift all coeffs out */
   {
      res->length = 0;
      return;
   }

   nmod_poly_fit_length(res, len);

   mpn_copyi(res->coeffs, poly->coeffs + k, len);
   
   res->length = len;
}
