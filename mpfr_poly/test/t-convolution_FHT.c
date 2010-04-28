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
   printf("convolution_FHT....");
   fflush(stdout);

   mpfr_poly_randinit();
   
   for (ulong i = 0; i < 10UL; i++) 
   {
      mpfr_poly_t a, b, c1, c2, d1, d2;
      ulong n = 10;//n_randint(10);
	  ulong length = (1L<<n);
	  ulong prec = 100;//n_randint(100) + 50*MPFR_PREC_MIN;
      
      mpfr_poly_init2(a, length, prec);
      mpfr_poly_init2(b, length, prec);
      /*mpfr_poly_init2(c1, length, prec);
      mpfr_poly_init2(d1, length, prec);
      mpfr_poly_init2(c2, length, prec);
      mpfr_poly_init2(d2, length, prec);*/
      
	  mpfr_poly_randtest(a, length);
      mpfr_poly_randtest(b, length);
      /*mpfr_poly_randtest(c1, length);
      mpfr_poly_randtest(d1, length);*/
      
	  /*for (ulong j = 0; j < length; j++)
		  mpfr_set(b->coeffs + j, a->coeffs + j, GMP_RNDN);

	  for (ulong j = 0; j < length; j++)
		  mpfr_set(c2->coeffs + j, c1->coeffs + j, GMP_RNDN);

	  for (ulong j = 0; j < length; j++)
		  mpfr_set(d2->coeffs + j, d1->coeffs + j, GMP_RNDN);*/

	  for (ulong j = 0; j < 100; j++)
		  _mpfr_poly_convolution_FHT(a->coeffs, b->coeffs, n, prec);
      /*_mpfr_poly_convolution_FHT(a->coeffs, d1->coeffs, n, prec);
      _mpfr_poly_convolution_FHT(b->coeffs, d2->coeffs, n, prec);
      _mpfr_poly_convolution_FHT(b->coeffs, c2->coeffs, n, prec);
      
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
	  }*/

      mpfr_poly_clear(a);
      mpfr_poly_clear(b);
     /* mpfr_poly_clear(c1);
      mpfr_poly_clear(d1);
      mpfr_poly_clear(c2);
      mpfr_poly_clear(d2);*/
   }
   
   mpfr_poly_randclear();
      
   printf("PASS\n");
   return 0;
}
