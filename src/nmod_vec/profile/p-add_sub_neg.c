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
   int type;
} info_t;

void sample(void * arg, ulong count)
{
   ulong n;
   nmod_t mod;
   info_t * info = (info_t *) arg;
   flint_bitcnt_t bits = info->bits;
   int type = info->type;
   slong j;
   slong i;
   nn_ptr vec1, vec2, res;
   FLINT_TEST_INIT(state);


   n = n_randbits(state, bits);
   if (n == UWORD(0)) n++;

   nmod_init(&mod, n);

   vec1 = _nmod_vec_init(1000);
   vec2 = _nmod_vec_init(1000);
   res = _nmod_vec_init(1000);

   for (j = 0; j < 1000; j++)
      vec1[j] = n_randint(state, n);

   for (j = 0; j < 1000; j++)
      vec2[j] = n_randint(state, n);

   switch (type)
   {
   case 1:
	  prof_start();
      for (i = 0; i < count; i++)
      {
         _nmod_vec_add(res, vec1, vec2, 1000, mod);
      }
      prof_stop();
      break;

   case 2:
	  prof_start();
      for (i = 0; i < count; i++)
      {
         _nmod_vec_sub(res, vec1, vec2, 1000, mod);
      }
      prof_stop();
      break;

   case 3:
	  prof_start();
      for (i = 0; i < count; i++)
      {
         _nmod_vec_neg(res, vec1, 1000, mod);
      }
      prof_stop();
      break;
   }

   _nmod_vec_clear(vec1);
   _nmod_vec_clear(vec2);
   FLINT_TEST_CLEAR(state);
}

int main(void)
{
   double min1, min2, min3, max;
   info_t info;
   flint_bitcnt_t i;

   for (i = 2; i <= FLINT_BITS; i++)
   {
      info.bits = i;

	  info.type = 1;
	  prof_repeat(&min1, &max, sample, (void *) &info);

	  info.type = 2;
	  prof_repeat(&min2, &max, sample, (void *) &info);

	  info.type = 3;
	  prof_repeat(&min3, &max, sample, (void *) &info);

      flint_printf("bits %wd, add = %.1lf c/l, sub = %.1lf c/l, neg = %.1lf c/l\n",
         i, (min1/(double)FLINT_CLOCK_SCALE_FACTOR)/1000,
            (min2/(double)FLINT_CLOCK_SCALE_FACTOR)/1000,
            (min3/(double)FLINT_CLOCK_SCALE_FACTOR)/1000
	  );
   }

   return 0;
}
