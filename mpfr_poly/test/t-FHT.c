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

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include <mpfr.h>
#include <math.h>
#include "flint.h"
#include "mpfr_poly.h"
#include "ulong_extras.h"

const ulong revtab1[2] = { 0, 1 };
const ulong revtab2[4] = { 0, 2, 1, 3 };
const ulong revtab3[8] = { 0, 4, 2, 6, 1, 5, 3, 7 };
const ulong revtab4[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };

const ulong * revtab[5] = { NULL, revtab1, revtab2, revtab3, revtab4 };

// computes the reverse binary of a binary number of the given number of bits
ulong revbin(ulong in, ulong bits)
{
   ulong out = 0, i;
   
   if (bits == 0) return 0;
   if (bits <= 4)
      return revtab[bits][in];

   for (i = 0; i < bits; i++)
   {   
      out <<= 1;
      out += (in & 1);
      in >>= 1;
   }

   return out;
}

int main(void)
{
   int result;
   printf("FHT....");
   fflush(stdout);

   mpfr_poly_randinit();
   
   for (ulong i = 0; i < 1000UL; i++) 
   {
      mpfr_poly_t a, b;
      ulong n = n_randint(10);
	  ulong length = (1L<<n);
	  ulong prec = n_randint(100) + 50*MPFR_PREC_MIN;
      
      mpfr_poly_init2(a, length, prec);
      mpfr_poly_init2(b, length, prec);
      mpfr_poly_randtest(a, length);
      
	  //_mpfr_vec_copy(b->coeffs, a->coeffs, length);
	  for (ulong j = 0; j < length; j++)
	     mpfr_set(b->coeffs + j, a->coeffs + j, GMP_RNDN);

	  _mpfr_poly_FHT(a->coeffs, n, prec);
      for (ulong j = 0; j < length; j++)
	  {
		  ulong s = revbin(j, n);
		  if (j < s)
			 mpfr_swap(a->coeffs + j, a->coeffs + s);
	  }

      _mpfr_poly_FHT(a->coeffs, n, prec);
      
	  for (ulong j = 0; j < length; j++)
	  {
		  ulong s = revbin(j, n);
		  if (j < s)
			  mpfr_swap(a->coeffs + j, a->coeffs + s);
	  }

	  for (ulong j = 0; j < length; j++)
	     mpfr_div_2exp(a->coeffs + j, a->coeffs + j, n, GMP_RNDN);

	  for (ulong j = 0; j < length; j++)
	  {
	     mpfr_sub(a->coeffs + j, a->coeffs + j, b->coeffs + j, GMP_RNDN);
         double d = mpfr_get_d(a->coeffs + j, GMP_RNDN);
		 if (fabs(d) > 0.1)
		 {
			 printf("d = %f\n", d);
			 printf("Error: length = %ld, prec = %ld\n", length, prec);
			 printf("Error in coeff %ld\n", j);
			 abort();
		 }
	  }

      mpfr_poly_clear(a);
      mpfr_poly_clear(b);
   }
   
   mpfr_poly_randclear();
      
   printf("PASS\n");
   return 0;
}
