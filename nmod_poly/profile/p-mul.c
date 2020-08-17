/*
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "profiler.h"
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

typedef struct
{
   slong n;
} info_t;

void sample(void * arg, ulong count)
{
   info_t * info = (info_t *) arg;
   slong n = info->n, i, j;
   slong scale;

   FLINT_TEST_INIT(state);
   
   nmod_poly_t a, b, c;
   mp_limb_t m = n_randint(state, 1<<((48-FLINT_BIT_COUNT(n))/2));
   if (m == 0) m = 2;
   
   nmod_poly_init2(a, m, n);
   nmod_poly_init2(b, m, n);
   nmod_poly_init2(c, m, 2*n - 1);
   
   for (i = 0; i < n; i++)
   {
       a->coeffs[i] = n_randint(state, m);
       b->coeffs[i] = n_randint(state, m);
   }
   a->length = n;
   b->length = n;

   scale = 1;
   if (n < 100000) scale = 10;
   if (n < 10000) scale = 100;
   if (n < 100) scale = 1000;
      
   for (i = 0; i < count; i++)
   {
	  prof_start();
      for (j = 0; j < scale; j++)
      {
          nmod_poly_mul(c, a, b);
      }
	  prof_stop();
      if (c->coeffs[n - 2] == 123) abort();
   }
  
   nmod_poly_clear(a);
   nmod_poly_clear(b);
   nmod_poly_clear(c);
   flint_randclear(state);
}

int main(void)
{
   double min, max;
   info_t info;
   slong k, scale;

   for (k = 2; k <= 30; k++)
   {
      info.n = 1<<k;

      scale = 1;
      if (info.n < 100000) scale = 10;
      if (info.n < 10000) scale = 100;
      if (info.n < 100) scale = 1000;
   
      prof_repeat(&min, &max, sample, (void *) &info);
         
      flint_printf("length %wd, min %.3g ms, max %.3g ms, norm %.3g\n", 
           info.n,
		   ((min/(double)FLINT_CLOCK_SCALE_FACTOR)/scale)/2400000.0,
           ((max/(double)FLINT_CLOCK_SCALE_FACTOR)/scale)/2400000.0,
           (((min/(double)FLINT_CLOCK_SCALE_FACTOR)/scale)/2400000.0)
              *500000.0/info.n/FLINT_BIT_COUNT(info.n)
	     );
   }

   return 0;
}
