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
#include <mpfr.h>
#include "flint.h"
#include "mpfr_vec.h"
#include "mpfr_poly.h"

void _mpfr_poly_convolution_trans(mpfr * coeffs1, mpfr * coeffs2, ulong n, mp_bitcnt_t prec)
{
   ulong len2 = (1UL<<n)/2;

   ulong i, j;

   mpfr_t temp, temp2, temp3, temp4;
   mpfr_init2(temp, prec);
   mpfr_init2(temp2, prec);
   mpfr_init2(temp3, prec);
   mpfr_init2(temp4, prec);

   for (i = 1, j = 2*len2 - 1; i < len2; i++, j--)
   {
      mpfr_add(temp, coeffs1 + i, coeffs1 + j, GMP_RNDN);
      mpfr_sub(temp2, coeffs1 + i, coeffs1 + j, GMP_RNDN);

	  mpfr_mul(temp3, coeffs2 + i, temp, GMP_RNDN);
	  mpfr_mul(temp4, coeffs2 + j, temp2, GMP_RNDN);
      mpfr_add(temp3, temp3, temp4, GMP_RNDN);
	  mpfr_div_2exp(coeffs1 + i, temp3, 1, GMP_RNDN);

	  mpfr_mul(temp3, coeffs2 + j, temp, GMP_RNDN);
	  mpfr_mul(temp4, coeffs2 + i, temp2, GMP_RNDN);
      mpfr_sub(temp3, temp3, temp4, GMP_RNDN);
	  mpfr_div_2exp(coeffs1 + j, temp3, 1, GMP_RNDN);
   }

   mpfr_mul(coeffs1, coeffs1, coeffs2, GMP_RNDN);
   if (n != 0) 
      mpfr_mul(coeffs1 + len2, coeffs1 + len2, coeffs2 + len2, GMP_RNDN);


   mpfr_clear(temp);
   mpfr_clear(temp2);
   mpfr_clear(temp3);
   mpfr_clear(temp4);
}
