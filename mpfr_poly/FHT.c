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

#include <mpir.h>
#include <mpfr.h>
#include "flint.h"
#include "mpfr_vec.h"
#include "mpfr_poly.h"

void _mpfr_poly_FHT_recursive(mpfr * coeffs, ulong n, mp_bitcnt_t prec, 
							          mpfr * costw, mpfr * sintw, ulong stride);

void _mpfr_poly_FHT_recursive(mpfr * coeffs, ulong n, mp_bitcnt_t prec, 
							          mpfr * costw, mpfr * sintw, ulong stride)
{
   ulong i;
   mpfr_t temp, temp2;
   mpfr * cos, * sin;
   mpfr_init2(temp, prec);
   mpfr_init2(temp2, prec);
  
   if (n == 0) return;

   ulong len2 = (1UL<<(n - 1));
   
   for (i = 0; i < len2; i++)
   {
      mpfr_add(temp, coeffs + i, coeffs + len2 + i, GMP_RNDN);
	  mpfr_sub(coeffs + len2 + i, coeffs + i, coeffs + len2 + i, GMP_RNDN);
      mpfr_swap(temp, coeffs + i);
   }
   
   if (len2 != 1)
   {
      for (i = 1; i < len2/4; i++)
	  {
         cos = costw + i*stride;
		 sin = sintw + i*stride;
		 
		 mpfr_mul(temp, coeffs + len2 + i, cos, GMP_RNDN);
         mpfr_mul(temp2, coeffs + 2*len2 - i, sin, GMP_RNDN);
         
         mpfr_mul(coeffs + len2 + i, coeffs + len2 + i, sin, GMP_RNDN);
         mpfr_mul(coeffs + 2*len2 - i, coeffs + 2*len2 - i, cos, GMP_RNDN);
         
		 mpfr_sub(coeffs + 2*len2 - i, coeffs + len2 + i, coeffs + 2*len2 - i, GMP_RNDN);
		 mpfr_add(coeffs + len2 + i, temp, temp2, GMP_RNDN);

         mpfr_mul(temp, coeffs + 3*len2/2 - i, sin, GMP_RNDN);
         mpfr_mul(temp2, coeffs + 3*len2/2 + i, cos, GMP_RNDN);
         
         mpfr_mul(coeffs + 3*len2/2 + i, coeffs + 3*len2/2 + i, sin, GMP_RNDN);
         mpfr_mul(coeffs + 3*len2/2 - i, coeffs + 3*len2/2 - i, cos, GMP_RNDN);
         
         mpfr_sub(coeffs + 3*len2/2 + i, coeffs + 3*len2/2 - i, coeffs + 3*len2/2 + i, GMP_RNDN);
		 mpfr_add(coeffs + 3*len2/2 - i, temp, temp2, GMP_RNDN);
	  }
      
	  if (len2 != 2)
	  {
		  sin = sintw + (len2/4)*stride;
          mpfr_add(temp, coeffs + 5*len2/4, coeffs + 7*len2/4, GMP_RNDN);
	      mpfr_sub(coeffs + 7*len2/4, coeffs + 5*len2/4, coeffs + 7*len2/4, GMP_RNDN);
          mpfr_swap(temp, coeffs + 5*len2/4);
	      mpfr_mul(coeffs + 5*len2/4, coeffs + 5*len2/4, sin, GMP_RNDN);
	      mpfr_mul(coeffs + 7*len2/4, coeffs + 7*len2/4, sin, GMP_RNDN);
	  }
   }

   mpfr_clear(temp);
   mpfr_clear(temp2);
 
   _mpfr_poly_FHT_recursive(coeffs, n - 1, prec, costw, sintw, stride*2);
   _mpfr_poly_FHT_recursive(coeffs + len2, n - 1, prec, costw, sintw, stride*2);
}

void _mpfr_poly_FHT(mpfr * coeffs, ulong n, mp_bitcnt_t prec)
{
	mpfr * cos, * sin;
    ulong length = (1UL<<n), i;
	mpfr_t t1, t2;

	if (n == 0) return;

	mpfr_init2(t1, prec);
	mpfr_init2(t2, prec);

	cos = _mpfr_vec_init(length/8 + 1, prec);
	sin = _mpfr_vec_init(length/8 + 1, prec);

	mpfr_set_ui(cos, 1, GMP_RNDN);
	mpfr_set_ui(sin, 0, GMP_RNDN);

	mpfr_const_pi(t1, GMP_RNDN);
	mpfr_div_2exp(t1, t1, n - 1, GMP_RNDN);
	mpfr_set(t2, t1, GMP_RNDN);

    for (i = 1; i < length/8; i++)
	{
		mpfr_sin_cos(sin + i, cos + i, t2, GMP_RNDN);
		mpfr_add(t2, t2, t1, GMP_RNDN);
	}

	if (n > 2) 
		mpfr_sin(sin + i, t2, GMP_RNDN);

	_mpfr_poly_FHT_recursive(coeffs, n, prec, cos, sin, 1);

	mpfr_clear(t1);
	mpfr_clear(t2);
    _mpfr_vec_clear(cos, length/8 + 1);
	_mpfr_vec_clear(sin, length/8 + 1);
}