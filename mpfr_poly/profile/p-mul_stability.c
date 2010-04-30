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
/******************************************************************************

 (C) 2010 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include <mpfr.h>
#include <math.h>
#include "flint.h"
#include "ulong_extras.h"
#include "mpfr_vec.h"
#include "mpfr_poly.h"

#define ITERS 30

ulong bits_lost(mpfr_poly_t a, mpfr_poly_t b)
{
   mpfr_t t;
   ulong i;
   mpfr_init2(t, a->prec);
   ulong lost = 0, bits, bits1, bits2;
   double d;

   for (i = 0; i < FLINT_MIN(2, b->length); i++)
   {
      d = mpfr_get_d_2exp(&bits1, b->coeffs + i, GMP_RNDN);
	  mpfr_sub(t, a->coeffs + i, b->coeffs + i, GMP_RNDN);
	  d = mpfr_get_d_2exp(&bits2, t, GMP_RNDN);
	  bits = a->prec - (bits1 - bits2);
	  if ((long) bits < 0L) bits = 0;
	  if (bits > lost)
		  lost = bits;
   }
   
   mpfr_clear(t);

   return lost;
}

void sample(mp_bitcnt_t * min, mp_bitcnt_t * av, mp_bitcnt_t * max, ulong length, ulong prec)
{
   ulong log_len = 0;
   while ((1UL<<log_len) < length) log_len++;
   ulong prec2 = 2*prec + 2*log_len + 1;
   ulong lost;

   mpfr_poly_t a, b, c, d;
   mpfr_poly_init2(a, length, prec);
   mpfr_poly_init2(b, length, prec);
   mpfr_poly_init2(c, length, prec);
   mpfr_poly_init2(d, length, prec);

   *min = -1L;
   *max = 0;
   *av = 0;

   for (ulong i = 0; i < ITERS; i++)
   {
      mpfr_poly_randtest(a, length);
      mpfr_poly_randtest(b, length);	  
	  //_mpfr_vec_scalar_mul_2exp(a->coeffs, a->coeffs, length, prec/2);
	  //_mpfr_vec_scalar_mul_2exp(b->coeffs, b->coeffs, length, prec/2);
	  for (ulong j = 0; j < length; j++)
	  {	  
		  mpfr_mul_2exp(a->coeffs + j, a->coeffs + j, n_randint(prec), GMP_RNDN);
	      mpfr_mul_2exp(b->coeffs + j, b->coeffs + j, n_randint(prec), GMP_RNDN);
	      if (n_randint(2)) mpfr_neg(a->coeffs + j, a->coeffs + j, GMP_RNDN);
		  if (n_randint(2)) mpfr_neg(b->coeffs + j, b->coeffs + j, GMP_RNDN);
	  }
	  
	  mpfr_poly_mul_FHT(c, a, b);
	  
	  mpfr_poly_set_prec(a, prec2);
	  mpfr_poly_set_prec(b, prec2);
	  mpfr_poly_set_prec(d, prec2);
	  mpfr_poly_mul_classical(d, a, b);
      mpfr_poly_set_prec(d, prec);
	  
	  lost = bits_lost(c, d);
	  if (lost < *min) *min = lost;
	  if (lost > *max) *max = lost;
	  (*av) += lost; 
   }
   
   (*av) /= ITERS;

   mpfr_poly_clear(a);
   mpfr_poly_clear(b);
   mpfr_poly_clear(c);
   mpfr_poly_clear(d);
}

int main(void)
{
   ulong min, av, max, length, len_old = 0, prec = 106;
   
   printf("p-mul_stability:\n");
   
   mpfr_poly_randinit();

   length = 1;

   for (ulong i = 0; i <= 100; i++)
   {
	  sample(&min, &av, &max, length, prec);
	  printf("len = %ld, min = %ld, av = %ld, max = %ld, prec = %ld\n", length, min, av, max, prec);

      len_old = length;
	  length = ceil((double) length * 1.3);
	  if (length == len_old) length++;
   }

   mpfr_poly_randclear();

   return 0;
}
