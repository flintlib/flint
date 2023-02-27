/*
    Copyright 2010 William Hart
    Copyright 2013 Martin Lee

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
   slong s;
   slong alg;
} info_t;

void sample(void * arg, ulong count)
{
   info_t * info = (info_t *) arg;
   slong n = info->n, i, j, s= info->s, alg=info->alg;
   slong scale;

   FLINT_TEST_INIT(state);
   
   nmod_poly_t a, b, c, d, dinv;
   mp_limb_t m = n_randint(state, 1<<((48-FLINT_BIT_COUNT(n))/2)); /* modulus */
   if (m == 0) m = 2;
   
   nmod_poly_init2(a, m, n);
   nmod_poly_init2(b, m, n);
   nmod_poly_init2(c, m, 2*n - 1);
   nmod_poly_init2(d, m, s);
   nmod_poly_init2(dinv, m, s);

   for (i = 0; i < n; i++)
   {
       a->coeffs[i] = n_randint(state, m);
       b->coeffs[i] = n_randint(state, m);
   }
   for (i = 0; i < s; i++)
       d->coeffs[i] = n_randint(state, m);
   a->length = n;
   b->length = n;
   d->length = s;

   nmod_poly_reverse(dinv, d, s);
   nmod_poly_inv_series(dinv, dinv, s);

   scale = 1;
   if (n < 100000) scale = 10;
   if (n < 10000) scale = 100;
   if (n < 100) scale = 1000;
      
   for (i = 0; i < count; i++)
   {
      if (alg == 1)
      {
	  prof_start();
          for (j = 0; j < scale; j++)
          {
              nmod_poly_mulmod_preinv(c, a, b, d, dinv);
          }
	  prof_stop();
      }
      else
      {
	  prof_start();
          for (j = 0; j < scale; j++)
          {
              nmod_poly_mulmod(c, a, b, d);
          }
	  prof_stop();
      }
      if (c->coeffs[n - 2] == 123) flint_abort();
   }
  
   nmod_poly_clear(a);
   nmod_poly_clear(b);
   nmod_poly_clear(c);
   nmod_poly_clear(d);
   nmod_poly_clear(dinv);
   flint_randclear(state);
}

int main(void)
{
   double min, max;
   info_t info;
   slong i, k, scale;


   for (k = 2; k <= 30; k++)
   {
      info.n = 1<<k;

      for (i= 0; i < 3; i++)
      {
        if (i == 0)
          info.s= info.n+1;
        else if (i== 1)
          info.s= info.n+ ((1<<(k+1))-(1<<k))/2*i;
        else if (i== 2)
          info.s= info.n*2-1;
      scale = 1;
      if (info.n < 100000) scale = 10;
      if (info.n < 10000) scale = 100;
      if (info.n < 100) scale = 1000;

      info.alg= 1;
      prof_repeat(&min, &max, sample, (void *) &info);
         
      flint_printf("length %wd, modulus degree %wd, min %.3g ms, max %.3g ms, norm %.3g\n", 
           info.n, info.s, 
		   ((min/(double)FLINT_CLOCK_SCALE_FACTOR)/scale)/2400000.0,
           ((max/(double)FLINT_CLOCK_SCALE_FACTOR)/scale)/2400000.0,
           (((min/(double)FLINT_CLOCK_SCALE_FACTOR)/scale)/2400000.0)
              *500000.0/info.n/FLINT_BIT_COUNT(info.n)
	     );

      info.alg= 2;
      prof_repeat(&min, &max, sample, (void *) &info);
         
      flint_printf("length %wd, modulus degree %wd, min %.3g ms, max %.3g ms, norm %.3g\n", 
           info.n, info.s, 
		   ((min/(double)FLINT_CLOCK_SCALE_FACTOR)/scale)/2400000.0,
           ((max/(double)FLINT_CLOCK_SCALE_FACTOR)/scale)/2400000.0,
           (((min/(double)FLINT_CLOCK_SCALE_FACTOR)/scale)/2400000.0)
              *500000.0/info.n/FLINT_BIT_COUNT(info.n)
	     );
      }
   }

   return 0;
}
