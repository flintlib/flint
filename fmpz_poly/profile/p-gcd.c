/*
    Copyright 2015 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "profiler.h"
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_poly.h"

#define LENGTH 100

typedef struct
{
   slong bits;
} info_t;

void random_fmpz_poly(fmpz_poly_t pol, flint_rand_t state, slong bits)
{
   fmpz * arr;
   slong i;

   fmpz_poly_fit_length(pol, LENGTH);

   arr = pol->coeffs;

   for (i = 0; i < LENGTH; i++)
      fmpz_randbits(arr + i, state, bits);

   _fmpz_poly_set_length(pol, LENGTH);
   _fmpz_poly_normalise(pol);
}

void sample(void * arg, ulong count)
{
   info_t * info = (info_t *) arg;
   slong bits = info->bits, i, j;
   int scale;
   
   scale = 100;
   if (bits >= 50) scale = 10;
   if (bits >= 500) scale = 4;
   
   FLINT_TEST_INIT(state);
   
   fmpz_poly_t p1, p2, a, b, c, g;
           
   fmpz_poly_init(p1);
   fmpz_poly_init(p2);
   fmpz_poly_init(a);
   fmpz_poly_init(b);
   fmpz_poly_init(c);
   fmpz_poly_init(g);

   for (i = 0; i < count; i++)
   {
      random_fmpz_poly(a, state, bits);
      random_fmpz_poly(b, state, bits);
      random_fmpz_poly(c, state, bits);
      
      fmpz_poly_mul(p1, a, c);
      fmpz_poly_mul(p2, b, c);

      prof_start();
      for (j = 0; j < scale; j++)
      {
         fmpz_poly_gcd(g, p1, p2);
      }
      prof_stop();
   }
  
   fmpz_poly_clear(p1);
   fmpz_poly_clear(p2);
   fmpz_poly_clear(a);
   fmpz_poly_clear(b);
   fmpz_poly_clear(c);
   fmpz_poly_clear(g);

   flint_randclear(state);
}

int main(void)
{
   double min, max;
   info_t info;
   slong k, scale;

   printf("fmpz_poly_gcd\n");
   flint_printf("length = %wd\n", LENGTH);

   for (k = 31; k <= 1000; k += 31)
   {
      info.bits = k;
     
      scale = 100;
      if (k >= 50) scale = 10;
      if (k >= 500) scale = 4;
      
      prof_repeat(&min, &max, sample, (void *) &info);
      
      flint_printf("bits %wd, min %.3g ms, max %.3g ms\n", 
           info.bits,
		   (min/scale)/1000.0,
                   (max/scale)/1000.0
	     );
   }

   return 0;
}
