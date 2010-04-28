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
#include "mpfr_poly.h"

void _mpfr_poly_FHT(mpfr * coeffs, ulong n, mp_bitcnt_t prec)
{
   ulong i;
   mpfr_t temp, temp2, temp3, temp4, cos, sin;
   mpfr_init2(temp, prec);
   mpfr_init2(temp2, prec);
   mpfr_init2(temp3, prec);
   mpfr_init2(temp4, prec);
   mpfr_init2(sin, prec);
   mpfr_init2(cos, prec);

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
      mpfr_const_pi(temp, GMP_RNDN);
	  mpfr_div_2exp(temp, temp, n - 1, GMP_RNDN);
	  mpfr_set(temp2, temp, GMP_RNDN);

	  for (i = 1; i < len2/4; i++)
	  {
         mpfr_sin_cos(sin, cos, temp2, GMP_RNDN);
		 
		 mpfr_mul(temp3, coeffs + len2 + i, cos, GMP_RNDN);
         mpfr_mul(temp4, coeffs + 2*len2 - i, sin, GMP_RNDN);
         
         mpfr_mul(coeffs + len2 + i, coeffs + len2 + i, sin, GMP_RNDN);
         mpfr_mul(coeffs + 2*len2 - i, coeffs + 2*len2 - i, cos, GMP_RNDN);
         
		 mpfr_sub(coeffs + 2*len2 - i, coeffs + len2 + i, coeffs + 2*len2 - i, GMP_RNDN);
		 mpfr_add(coeffs + len2 + i, temp3, temp4, GMP_RNDN);

         mpfr_mul(temp3, coeffs + 3*len2/2 - i, sin, GMP_RNDN);
         mpfr_mul(temp4, coeffs + 3*len2/2 + i, cos, GMP_RNDN);
         
         mpfr_mul(coeffs + 3*len2/2 + i, coeffs + 3*len2/2 + i, sin, GMP_RNDN);
         mpfr_mul(coeffs + 3*len2/2 - i, coeffs + 3*len2/2 - i, cos, GMP_RNDN);
         
         mpfr_sub(coeffs + 3*len2/2 + i, coeffs + 3*len2/2 - i, coeffs + 3*len2/2 + i, GMP_RNDN);
		 mpfr_add(coeffs + 3*len2/2 - i, temp3, temp4, GMP_RNDN);

	     mpfr_add(temp2, temp2, temp, GMP_RNDN);   
	  }

	  if (len2 != 2)
	  {
		  mpfr_sin(sin, temp2, GMP_RNDN);
	      mpfr_add(temp3, coeffs + 5*len2/4, coeffs + 7*len2/4, GMP_RNDN);
	      mpfr_sub(coeffs + 7*len2/4, coeffs + 5*len2/4, coeffs + 7*len2/4, GMP_RNDN);
          mpfr_swap(temp3, coeffs + 5*len2/4);
	      mpfr_mul(coeffs + 5*len2/4, coeffs + 5*len2/4, sin, GMP_RNDN);
	      mpfr_mul(coeffs + 7*len2/4, coeffs + 7*len2/4, sin, GMP_RNDN);
	  }
   }

   mpfr_clear(temp);
   mpfr_clear(temp2);
   mpfr_clear(temp3);
   mpfr_clear(temp4);
   mpfr_clear(cos);
   mpfr_clear(sin);

   _mpfr_poly_FHT(coeffs, n - 1, prec);
   _mpfr_poly_FHT(coeffs + len2, n - 1, prec);
}

