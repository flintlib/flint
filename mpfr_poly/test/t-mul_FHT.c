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

int main(void)
{
   int result;
   printf("mul_FHT....");
   fflush(stdout);

   mpfr_poly_randinit();
   
   for (ulong i = 0; i < 1000UL; i++) 
   {
      mpfr_poly_t a, b, c, d;
      ulong len1 = n_randint(200) + 1;
	  ulong len2 = n_randint(200) + 1;
	  ulong len_out = len1 + len2 - 1;
	  ulong prec = n_randint(100) + 50*MPFR_PREC_MIN;
      
      mpfr_poly_init2(a, len1, prec);
      mpfr_poly_init2(b, len2, prec);
      mpfr_poly_init2(c, len_out, prec);
      mpfr_poly_init2(d, len_out, prec);
      
	  mpfr_poly_randtest(a, len1);
      mpfr_poly_randtest(b, len2);
      
	  mpfr_poly_mul_FHT(c, a, b);
      mpfr_poly_mul_classical(d, a, b);
      
	  for (ulong j = 0; j < len_out; j++)
	  {
	     mpfr_sub(c->coeffs + j, c->coeffs + j, d->coeffs + j, GMP_RNDN);
         double d = mpfr_get_d(c->coeffs + j, GMP_RNDN);
		 if (fabs(d) > 0.1)
		 {
			 printf("d = %f\n", d);
			 printf("Error: len1 = %ld, len2 = %ld, prec = %ld\n", len1, len2, prec);
			 printf("Error in coeff %ld\n", j);
			 abort();
		 }
	  }

      mpfr_poly_clear(a);
      mpfr_poly_clear(b);
      mpfr_poly_clear(c);
      mpfr_poly_clear(d);
   }
   
  /*mpfr_poly_t a, b;
  ulong prec = 53;
      
  mpfr_poly_init2(a, 1, prec);
  mpfr_poly_init2(b, 2, prec);
      
  mpfr_set_ui(a->coeffs, 1, GMP_RNDN);
  a->length = 1;

  mpfr_set_ui(b->coeffs, 1, GMP_RNDN);
  mpfr_set_ui(b->coeffs + 1, 1, GMP_RNDN);
  b->length = 2;

  for (ulong j = 0; j < 10000; j++)
	 mpfr_poly_mul_classical(a, a, b);

  for (ulong j = 0; j < 100; j++)
     //mpfr_poly_mul_classical(b, a, a);
     mpfr_poly_mul_FHT(b, a, a);

   mpfr_poly_clear(a);
   mpfr_poly_clear(b);*/

   mpfr_poly_randclear();
      
   printf("PASS\n");
   return 0;
}
