/*
    Copyright 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "profiler.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_poly.h"

#define BITS 100

typedef struct
{
   int algo;
   slong length;
} info_t;

void random_fmpz_poly(fmpz_poly_t pol, flint_rand_t state, slong length)
{
   fmpz * arr;
   slong i;

   fmpz_poly_fit_length(pol, length);

   arr = pol->coeffs;

   for (i = 0; i < length; i++)
      fmpz_randbits(arr + i, state, BITS);

   _fmpz_poly_set_length(pol, length);
   _fmpz_poly_normalise(pol);
}

void sample(void * arg, ulong count)
{
   info_t * info = (info_t *) arg;
   slong length = info->length, i, j;
   int algo = info->algo;
   int scale;

   fmpz_poly_t p1, p2, b, q, r;
   fmpz_poly_powers_precomp_t pinv;

   FLINT_TEST_INIT(state);

   scale = 100;
   if (length >= 50) scale = 10;
   if (length >= 500) scale = 4;

   fmpz_poly_init(p1);
   fmpz_poly_init(p2);
   fmpz_poly_init(b);
   fmpz_poly_init(q);
   fmpz_poly_init(r);

   for (i = 0; i < count; i++)
   {
      random_fmpz_poly(p1, state, 2*length - 1);
      random_fmpz_poly(p2, state, length);
      fmpz_set_ui(p2->coeffs + p2->length - 1, 1); /* p2 must be monic */

      fmpz_poly_powers_precompute(pinv, p2);

      if (algo)
      {
         prof_start();
         for (j = 0; j < scale; j++)
         {
            fmpz_poly_rem_powers_precomp(r, p1, p2, pinv);
         }
	      prof_stop();
      } else
      {
         prof_start();
         for (j = 0; j < scale; j++)
         {
            fmpz_poly_rem(r, p1, p2);
         }
	      prof_stop();
      }

      fmpz_poly_powers_clear(pinv);
   }

   fmpz_poly_clear(p1);
   fmpz_poly_clear(p2);
   fmpz_poly_clear(b);
   fmpz_poly_clear(q);
   fmpz_poly_clear(r);

   FLINT_TEST_CLEAR(state);
}

int main(void)
{
   double min, max;
   info_t info;
   slong k, scale;

   printf("fmpz_poly_rem with precomputed powers\n");
   flint_printf("bits = %wd\n", BITS);

   for (k = 1; k <= 1000; k = (slong) ceil(1.1*k))
   {
      info.length = k;
      info.algo = 0;

      scale = 100;
      if (k >= 50) scale = 10;
      if (k >= 500) scale = 4;

      prof_repeat(&min, &max, sample, (void *) &info);

      flint_printf("Standard: length %wd, min %.3g ms, max %.3g ms\n",
           info.length,
		   ((min/(double)FLINT_CLOCK_SCALE_FACTOR)/scale)/2400000.0,
           ((max/(double)FLINT_CLOCK_SCALE_FACTOR)/scale)/2400000.0
	     );

     info.algo = 1;

     prof_repeat(&min, &max, sample, (void *) &info);

      flint_printf("With powers: length %wd, min %.3g ms, max %.3g ms\n\n",
           info.length,
		   ((min/(double)FLINT_CLOCK_SCALE_FACTOR)/scale)/2400000.0,
           ((max/(double)FLINT_CLOCK_SCALE_FACTOR)/scale)/2400000.0
	     );
   }

   return 0;
}
