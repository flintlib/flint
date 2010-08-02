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

   Copyright (C) 2010 William Hart
   
*****************************************************************************/

#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

void _nmod_poly_sub(mp_ptr res, mp_srcptr poly1, long len1, mp_srcptr poly2, long len2, nmod_t mod)
{
   long longer = FLINT_MAX(len1, len2);
   long shorter = FLINT_MIN(len1, len2);
   long i;
   
   _nmod_vec_sub(res, poly1, poly2, shorter, mod);
   
   if (poly1 != res) // copy any remaining coefficients from poly1
      for (i = shorter; i < len1; i++)
         res[i] = poly1[i];

   // careful, it is *always* necessary to negate coeffs from poly2, even if this is already res
	for (i = shorter; i < len2; i++) 
      res[i] = nmod_neg(poly2[i], mod);
}

void nmod_poly_sub(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2)
{
   long longer = FLINT_MAX(poly1->length, poly2->length);

   nmod_poly_fit_length(res, longer);
   
   _nmod_poly_sub(res->coeffs, poly1->coeffs, poly1->length, poly2->coeffs, poly2->length, poly1->mod);

   res->length = longer;
   _nmod_poly_normalise(res); // there may have been cancellation
}
