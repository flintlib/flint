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
#define PREC 300
#define GB 8 // guard bits
#define HAVE_MPFRCX 0 // whether or not we are linked against mpfrcx and mpc

//pick one of these only
#define FATEMAN1 0 // A = (x^10-10^100)/(x-10), B = (x-10)
#define FATEMAN2 1 // A = (x+1)*(sum(x^i,i,1,n) , B = (x-1)*sum((-1)^i*x^i,i,1,n) 
#define EXP 0 // A = B = sum 1 + x + x^2/2 + x^3/6 + ... + x^n/n! + ... 
#define RANDOM 0 // uniform random coefficients (modified by the flags below)

// these affect RANDOM only
#define SIGN 0 // coeffs are signed
#define VARY 0 // coefficient magnitudes vary

#if RANDOM || EXP
#define PREC_SCALE 10 // perform check at prec*10
#else
#define PREC_SCALE 1000
#endif

#if EXP 
#define PREC_SCALE2 5 // work at prec*5
#elif FATEMAN1 || FATEMAN2
#define PREC_SCALE2 1
#elif RANDOM 
#define PREC_SCALE2 3
#endif

#define PRINT 0 // print output polys

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

#if HAVE_MPFRCX
void mul_mpfrcx(mpfr_poly_t c, mpfr_poly_t a, mpfr_poly_t b, mpfr_poly_t temp)
{
	ulong len_out = a->length + b->length - 1;
	mpfr_poly_fit_length(c, len_out);
	mpfrx_array_mul_fft (c->coeffs, a->coeffs, b->coeffs, a->length, b->length);
	mpfrx_array_mul_karatsuba (c->coeffs, a->coeffs, b->coeffs, a->length, b->length, 1, 1, temp->coeffs);
	c->length = len_out;
}
#endif

void sample(mp_bitcnt_t * min, mp_bitcnt_t * av, mp_bitcnt_t * max, ulong length, ulong prec)
{
   ulong prec2 = PREC_SCALE*prec; // full precision computation occurs at this precision
   ulong prec3 = PREC_SCALE2*prec + GB; // prec at which we'll do our computation (including guard bits)
   ulong lost;
   
   mpfr_t scale, mult;
   mpfr_poly_t a, b, c, d, temp;

   mpfr_poly_init2(a, length + 1, prec);
   mpfr_poly_init2(b, length + 1, prec);
   mpfr_poly_init2(c, 2*length + 1, prec);
   mpfr_poly_init2(d, 2*length + 1, prec);
   mpfr_poly_init2(temp, 2*length + 12, prec3);
  
   mpfr_init2(scale, prec3);
   mpfr_init2(mult, prec3);
      
   *min = -1L;
   *max = 0;
   *av = 0;

   for (ulong i = 0; i < ITERS; i++)
   {
      mpfr_poly_set_prec(a, prec3); // compute at prec3 (includes guard bits)
	  mpfr_poly_set_prec(b, prec3);
	  mpfr_poly_set_prec(c, prec3);

#if EXP
	  mpfr_set_ui(a->coeffs, 1, GMP_RNDN);
	  mpfr_set_ui(b->coeffs, 1, GMP_RNDN);  
	  
	  if (length > 1)
	  {
		 mpfr_set_ui(a->coeffs + 1, 1, GMP_RNDN);
		 mpfr_set_ui(b->coeffs + 1, 1, GMP_RNDN);
	  }
	  
	  for (ulong j = 2; j < length; j++)
	  {
		  mpfr_div_ui(a->coeffs + j, a->coeffs + j - 1, j, GMP_RNDN);
		  mpfr_div_ui(b->coeffs + j, b->coeffs + j - 1, j, GMP_RNDN);
	  }

	  a->length = length;
	  b->length = length;
#endif

#if FATEMAN1
	  mpfr_set_ui(scale, 10, GMP_RNDN);
	  mpfr_set_ui(mult, 1, GMP_RNDN);
	  for (ulong j = 0; j < length; j++)
	  {
         mpfr_set(a->coeffs + length - j - 1, mult, GMP_RNDN);
		 mpfr_mul(mult, mult, scale, GMP_RNDN);
	  }
	  
	  mpfr_set_ui(b->coeffs + 1, 1, GMP_RNDN);
	  mpfr_set_si(b->coeffs, -10, GMP_RNDN);

	  a->length = length;
	  b->length = 2;
#endif

#if FATEMAN2
	  mpfr_set_ui(a->coeffs, 1, GMP_RNDN);
	  for (ulong j = 1; j < length; j++)
         mpfr_set_ui(a->coeffs + j, 2, GMP_RNDN);
      mpfr_set_ui(a->coeffs + length, 1, GMP_RNDN);

	  mpfr_set_ui(b->coeffs, 1, GMP_RNDN);
	  for (ulong j = 1; j < length; j++)
	  {
		  if (j & 1)
			  mpfr_set_si(b->coeffs + j, -2, GMP_RNDN);
		  else
			  mpfr_set_si(b->coeffs + j, 2, GMP_RNDN);
	  }
	  if (length & 1)
	     mpfr_set_si(b->coeffs + length, -2, GMP_RNDN);
	  else
	     mpfr_set_si(b->coeffs + length, 2, GMP_RNDN);

	  a->length = length + 1;
	  b->length = length + 1;
#endif

#if RANDOM
	  mpfr_poly_randtest(a, length);
      mpfr_poly_randtest(b, length);
	  for (ulong j = 0; j < length; j++)
	  {	  
#if VARY
		  mpfr_mul_2exp(a->coeffs + j, a->coeffs + j, n_randint(prec), GMP_RNDN);
	      mpfr_mul_2exp(b->coeffs + j, b->coeffs + j, n_randint(prec), GMP_RNDN);
#endif
#if SIGN
		  if (j&1) mpfr_neg(a->coeffs + j, a->coeffs + j, GMP_RNDN);
		  if (j&1) mpfr_neg(b->coeffs + j, b->coeffs + j, GMP_RNDN);
#endif
	  }
#endif

	  mpfr_poly_mul(c, a, b, GB);
	  //mpfr_poly_mul_classical(c, a, b);
	  //mpfr_poly_mul_FHT(c, a, b);
      //mpfrcx_mul(c, a, b, temp); // requires flint to be linked with mpfrcx and mpc

	  mpfr_poly_set_prec(c, prec);
	  
#if PRINT
	  ulong bits;
	  double d1;
	  for (ulong k = 0; k < 2*length + 1; k++)
	  {
		  d1 = mpfr_get_d_2exp(&bits, c->coeffs + k, GMP_RNDN);
		  printf("d = %lf, b = %ld, ", d1, bits);
	  }
	  printf("\n\n");
#endif

	  mpfr_poly_set_prec(a, prec2);
	  mpfr_poly_set_prec(b, prec2);
	  mpfr_poly_set_prec(d, prec2);
	  mpfr_poly_mul_classical(d, a, b);
      mpfr_poly_set_prec(a, prec);
	  mpfr_poly_set_prec(b, prec);
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
   mpfr_poly_clear(temp);

   mpfr_clear(scale);
   mpfr_clear(mult);
}

int main(void)
{
   ulong min, av, max, length, len_old = 0, prec = PREC;
   
   printf("p-mul_stability:\n");
   
   mpfr_poly_randinit();

   length = 1;

   for (ulong i = 0; i <= 100; i++)
   {
	  sample(&min, &av, &max, length, prec);
	  printf("len = %ld, min = %ld, av = %ld, max = %ld, prec = %ld\n", length, min, av, max, prec);

      len_old = length;
#if RANDOM
	  length = ceil((double) length * 1.3);
#endif
	  if (length == len_old) length++;
   }

   mpfr_poly_randclear();

   return 0;
}
