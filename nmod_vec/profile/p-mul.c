/*
    Copyright 2010 William Hart

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
#include "nmod_vec.h"

typedef struct
{
   flint_bitcnt_t bits;
   int fullword;
} info_t;

void sample(void * arg, ulong count)
{
   mp_limb_t n;
   nmod_t mod;
   info_t * info = (info_t *) arg;
   flint_bitcnt_t bits = info->bits;
   int fullword = info->fullword;
   mp_size_t j;
   slong i;
   FLINT_TEST_INIT(state);
      
   n = n_randbits(state, bits);
   if (n == UWORD(0)) n++;
      
   nmod_init(&mod, n);

   mp_ptr vec1 = _nmod_vec_init(1000);
   mp_ptr vec2 = _nmod_vec_init(1000);
   mp_ptr res = _nmod_vec_init(1000);
     
   for (j = 0; j < 1000; j++)
      vec1[j] = n_randint(state, n);

   for (j = 0; j < 1000; j++)
      vec2[j] = n_randint(state, n);

   switch (fullword)
   {
   case 0:
      prof_start();
      for (i = 0; i < count; i++)
         for (j = 0; j < 1000; j++)
           res[j] = nmod_mul(vec1[j], vec2[j], mod);
      prof_stop();
      break;

   case 1:
      prof_start();
      for (i = 0; i < count; i++)
         for (j = 0; j < 1000; j++)
           res[j] = _nmod_mul_fullword(vec1[j], vec2[j], mod);
      prof_stop();
      break;
   }

   flint_randclear(state);
   _nmod_vec_clear(vec1);
   _nmod_vec_clear(vec2);
}

int main(void)
{
   double min1, min2, max;
   info_t info;
   flint_bitcnt_t i;

   for (i = 2; i <= FLINT_BITS; i++)
   {
      info.bits = i;

      info.fullword = 0;
      prof_repeat(&min1, &max, sample, (void *) &info);

      if (i != FLINT_BITS)
      {
        flint_printf("bits %wd, mul = %.1lf c/l\n", 
           i, (min1/(double)FLINT_CLOCK_SCALE_FACTOR)/1000);
      }
      else
      {
        info.fullword = 1;
        prof_repeat(&min2, &max, sample, (void *) &info);

        flint_printf("bits %wd, mul = %.1lf c/l, mul_fullword = %.1lf c/l\n", 
           i, (min1/(double)FLINT_CLOCK_SCALE_FACTOR)/1000,
              (min2/(double)FLINT_CLOCK_SCALE_FACTOR)/1000);
      }
   }

   return 0;
}
