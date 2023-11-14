/*
    Copyright 2015 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "profiler.h"
#include "flint.h"
#include "ulong_extras.h"

typedef struct
{
   int algo;
} algo_s;

#define invert_limb_naive(ninv, n)                    \
   do {                                               \
      mp_limb_t dummy;                                \
      udiv_qrnnd (ninv, dummy, ~(n), ~(WORD(0)), n);  \
   } while (0)

void sample(void * arg, ulong count)
{
   ulong * array = (ulong *) flint_malloc(200 * sizeof(ulong));
   ulong i, ninv, sum = 0;
   slong j;
   algo_s * alg = (algo_s *) arg;

   FLINT_TEST_INIT(state);

   for (i = 0; i < count; i++)
   {
      for (j = 0; j < 200; j++)
      {
         array[j] = n_randlimb(state) | (UWORD(1) << (FLINT_BITS - 1));
      }

      prof_start();
      if (alg->algo == 0)
      {
         for (j = 0; j < 200; j++)
         {
            ninv = n_preinvert_limb_prenorm(array[j]);
            sum += ninv;
         }
      } else
      {
         for (j = 0; j < 200; j++)
         {
            invert_limb_naive(ninv, array[j]);
            sum += ninv;
         }
      }
      prof_stop();

      if (sum == 0) flint_printf("\r");
   }

   flint_randclear(state);
   flint_free(array);
}

int main(void)
{
   double min, max;
   algo_s alg;

   alg.algo = 0;

   prof_repeat(&min, &max, sample, &alg);

   flint_printf("invert_limb min time is %.3f cycles, max time is %.3f cycles\n",
           min/FLINT_CLOCK_SCALE_FACTOR/200, max/FLINT_CLOCK_SCALE_FACTOR/200);

   alg.algo = 1;

   prof_repeat(&min, &max, sample, &alg);

   flint_printf("invert_limb_naive min time is %.3f cycles, max time is %.3f cycles\n",
           min/FLINT_CLOCK_SCALE_FACTOR/200, max/FLINT_CLOCK_SCALE_FACTOR/200);

   return 0;
}
