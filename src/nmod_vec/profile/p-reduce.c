/*
    Copyright 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "profiler.h"
#include "ulong_extras.h"
#include "nmod.h"
#include "nmod_vec.h"

typedef struct
{
   flint_bitcnt_t bits;
   slong len;
   int lazy;
} info_t;

void sample(void * arg, ulong count)
{
   ulong n;
   nmod_t mod;
   info_t * info = (info_t *) arg;
   flint_bitcnt_t bits = info->bits;
   slong len = info->len;
   nn_ptr vec = _nmod_vec_init(len);
   nn_ptr vec2 = _nmod_vec_init(len);
   slong j;
   slong i;
   FLINT_TEST_INIT(state);

   for (j = 0; j < len; j++)
      vec[j] = n_randlimb(state);

   n = n_randbits(state, bits);
   if (n == UWORD(0)) n++;
   nmod_init(&mod, n);

   prof_start();
   if (info->lazy)
      for (i = 0; i < count; i++)
         _nmod_vec_reduce_lazy(vec2, vec, len, mod);
   else
      for (i = 0; i < count; i++)
         _nmod_vec_reduce(vec2, vec, len, mod);
   prof_stop();

   _nmod_vec_clear(vec);
   _nmod_vec_clear(vec2);
   FLINT_TEST_CLEAR(state);
}

int main(void)
{
   double min, minl, max;
   info_t info;
   flint_bitcnt_t i;
   slong lens[4] = {8, 32, 256, 1000};
   int k;

   flint_printf("cycles per entry (c/l), random words, by bit size of n and length;\n"
                "up to 62 bits, then _nmod_vec_reduce_lazy\n");
   for (i = 2; i <= FLINT_BITS; i++)
   {
      info.bits = i;
      flint_printf("bits %2wd:", i);
      for (k = 0; k < 4; k++)
      {
         info.len = lens[k];
         info.lazy = 0;
         prof_repeat(&min, &max, sample, (void *) &info);
         flint_printf("  len %4wd %5.2lf", lens[k],
               (min / (double) FLINT_CLOCK_SCALE_FACTOR) / lens[k]);
         if (i <= 62)
         {
            info.lazy = 1;
            prof_repeat(&minl, &max, sample, (void *) &info);
            flint_printf(" %5.2lf", (minl / (double) FLINT_CLOCK_SCALE_FACTOR) / lens[k]);
         }
      }
      flint_printf("\n");
   }

   return 0;
}
