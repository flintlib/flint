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

ulong bits_lost(mpfr_poly_t a, mpfr_poly_t b)
{
   mpfr_t t;
   ulong i;
   mpfr_init2(t, a->prec);
   ulong lost = 0, bits, bits1, bits2;
   double d1, d2;

   for (i = 0; i < b->length; i++)
   {
      d1 = mpfr_get_d_2exp(&bits1, b->coeffs + i, GMP_RNDN);
	  mpfr_sub(t, a->coeffs + i, b->coeffs + i, GMP_RNDN);
	  d2 = mpfr_get_d_2exp(&bits2, t, GMP_RNDN);
	  bits = a->prec - (bits1 - bits2);
	  if (d2 == 0.0) bits = 0;
	  if ((long) bits < 0L) bits = 0;
	  if (bits > lost)
		  lost = bits;
   }
   
   mpfr_clear(t);

   return lost;
}

int main(void)
{
   int result;
   printf("mul....");
   fflush(stdout);

   mpfr_poly_randinit();
   
   for (ulong i = 0; i < 1000UL; i++) 
   {
      mpfr_poly_t a, b, c, d;
      ulong len1 = n_randint(200) + 1;
	  ulong len2 = n_randint(200) + 1;
	  ulong len_out = len1 + len2 - 1;
	  ulong prec = n_randint(100) + 53;
      
      mpfr_poly_init2(a, len1, 2*prec);
      mpfr_poly_init2(b, len2, 2*prec);
      mpfr_poly_init2(c, len_out, 2*prec);
      mpfr_poly_init2(d, len_out, 2*prec);
      
	  mpfr_poly_randtest(a, len1);
      mpfr_poly_randtest(b, len2);
      
	  mpfr_poly_mul(c, a, b, 8);
	  mpfr_poly_mul_classical(d, a, b);
      
	  mpfr_poly_set_prec(c, prec);
      mpfr_poly_set_prec(d, prec);

	  if ((bits_lost(c, d) > 10) || (c->length != d->length))
	  {
		  printf("Error: mul_classical and mul don't agree within 10 bits\n");
		  abort();
	  }

      mpfr_poly_clear(a);
      mpfr_poly_clear(b);
      mpfr_poly_clear(c);
      mpfr_poly_clear(d);
   }
   
   mpfr_poly_randclear();
      
   printf("PASS\n");
   return 0;
}
