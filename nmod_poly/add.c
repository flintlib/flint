/*============================================================================

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

===============================================================================*/
/****************************************************************************

   Copyright (C) 2008, 2009 William Hart
   
*****************************************************************************/

#include <mpir.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

void _nmod_poly_add(mp_ptr res, mp_srcptr poly1, ulong len1, mp_srcptr poly2, ulong len2, nmod_t mod)
{
   ulong longer = FLINT_MAX(len1, len2);
   ulong shorter = FLINT_MIN(len1, len2);
   ulong i;
   
   _nmod_vec_add(res, poly1, poly2, shorter, mod);
   
   if (poly1 != res) // copy any remaining coefficients from poly1
      for (i = shorter; i < len1; i++)
         res[i] = poly1[i];

   if (poly2 != res) // copy any remaining coefficients from poly2
      for (i = shorter; i < len2; i++)
         res[i] = poly2[i];
}

void nmod_poly_add(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2)
{
   ulong longer = FLINT_MAX(poly1->length, poly2->length);

   nmod_poly_fit_length(res, longer);
	
   _nmod_poly_add(res->coeffs, poly1->coeffs, poly1->length, poly2->coeffs, poly2->length, poly1->mod);
    
   res->length = longer;
   _nmod_poly_normalise(res); // there may have been cancellation
}
