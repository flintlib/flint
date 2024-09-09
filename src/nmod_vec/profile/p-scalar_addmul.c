/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2024 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "profiler.h"
#include "ulong_extras.h"
#include "nmod_vec.h"

typedef struct
{
   flint_bitcnt_t bits;
   slong length;
} info_t;

void sample(void * arg, ulong count)
{
   ulong n, c;
   nmod_t mod;
   info_t * info = (info_t *) arg;
   flint_bitcnt_t bits = info->bits;
   slong length = info->length;
   slong i, j;
   nn_ptr vec = _nmod_vec_init(length);
   nn_ptr vec2 = _nmod_vec_init(length);
   FLINT_TEST_INIT(state);


   for (i = 0; i < count; i++)
   {
      n = n_randbits(state, bits);
      if (n == UWORD(0)) n++;
      c = n_randint(state, n);
      for (j = 0; j < length; j++)
         vec[j] = n_randint(state, n);

      nmod_init(&mod, n);

      prof_start();
      for (j = 0; j < 100; j++)
         _nmod_vec_scalar_addmul_nmod(vec2, vec, length, c, mod);
      prof_stop();
   }

   _nmod_vec_clear(vec);
   _nmod_vec_clear(vec2);
   FLINT_TEST_CLEAR(state);
}

void sample_shoup(void * arg, ulong count)
{
   ulong n, c;
   nmod_t mod;
   info_t * info = (info_t *) arg;
   flint_bitcnt_t bits = info->bits;
   slong length = info->length;
   slong i, j;
   nn_ptr vec = _nmod_vec_init(length);
   nn_ptr vec2 = _nmod_vec_init(length);
   FLINT_TEST_INIT(state);


   for (i = 0; i < count; i++)
   {
      n = n_randbits(state, bits);
      if (n == UWORD(0)) n++;
      c = n_randint(state, n);
      for (j = 0; j < length; j++)
         vec[j] = n_randint(state, n);

      nmod_init(&mod, n);

      prof_start();
      for (j = 0; j < 100; j++)
         _nmod_vec_scalar_addmul_nmod_shoup(vec2, vec, length, c, mod);
      prof_stop();
   }

   _nmod_vec_clear(vec);
   _nmod_vec_clear(vec2);
   FLINT_TEST_CLEAR(state);
}

void sample_generic(void * arg, ulong count)
{
   ulong n, c;
   nmod_t mod;
   info_t * info = (info_t *) arg;
   flint_bitcnt_t bits = info->bits;
   slong length = info->length;
   slong i, j;
   nn_ptr vec = _nmod_vec_init(length);
   nn_ptr vec2 = _nmod_vec_init(length);
   FLINT_TEST_INIT(state);


   for (i = 0; i < count; i++)
   {
      n = n_randbits(state, bits);
      if (n == UWORD(0)) n++;
      c = n_randint(state, n);
      for (j = 0; j < length; j++)
         vec[j] = n_randint(state, n);

      nmod_init(&mod, n);

      prof_start();
      for (j = 0; j < 100; j++)
         _nmod_vec_scalar_addmul_nmod_generic(vec2, vec, length, c, mod);
      prof_stop();
   }

   _nmod_vec_clear(vec);
   _nmod_vec_clear(vec2);
   FLINT_TEST_CLEAR(state);
}

int main(void)
{
   double min, max;
   double mins[18]; // note: max seems to be consistently identical or extremely close to min
   double mins_shoup[18];
   double mins_generic[18];
   info_t info;
   flint_bitcnt_t i;

   flint_printf("unit: all measurements in c/l\n");
   flint_printf("profiled: general branching | precomp shoup | generic\n");
   flint_printf("bit/len\t");
   for (int len = 1; len <= 16; ++len)
       flint_printf("%d\t\t", len);
   flint_printf("1024\t\t");
   flint_printf("65536\n");

   for (i = 2; i <= FLINT_BITS; i++)
   {
      info.bits = i;

      for (int len = 1; len <= 16; ++len)
      {
          info.length = len;

          prof_repeat(&min, &max, sample, (void *) &info);
          mins[len-1] = min;
          prof_repeat(&min, &max, sample_shoup, (void *) &info);
          mins_shoup[len-1] = min;
          prof_repeat(&min, &max, sample_generic, (void *) &info);
          mins_generic[len-1] = min;
      }

      info.length = 1024;
      prof_repeat(&min, &max, sample, (void *) &info);
      mins[16] = min;
      prof_repeat(&min, &max, sample_shoup, (void *) &info);
      mins_shoup[16] = min;
      prof_repeat(&min, &max, sample_generic, (void *) &info);
      mins_generic[16] = min;

      info.length = 65536;
      prof_repeat(&min, &max, sample, (void *) &info);
      mins[17] = min;
      prof_repeat(&min, &max, sample_shoup, (void *) &info);
      mins_shoup[17] = min;
      prof_repeat(&min, &max, sample_generic, (void *) &info);
      mins_generic[17] = min;

      if (i < FLINT_BITS)
      {
          flint_printf("%wd", i);
          for (int len = 1; len <= 16; ++len)
              flint_printf("\t%.1lf|%.1lf|%.1lf",
                      (mins[len-1]/(double)FLINT_CLOCK_SCALE_FACTOR)/(len*100),
                      (mins_shoup[len-1]/(double)FLINT_CLOCK_SCALE_FACTOR)/(len*100),
                      (mins_generic[len-1]/(double)FLINT_CLOCK_SCALE_FACTOR)/(len*100));
          flint_printf("\t%.1lf|%.1lf|%.1lf",
                  (mins[16]/(double)FLINT_CLOCK_SCALE_FACTOR)/(1024*100),
                  (mins_shoup[16]/(double)FLINT_CLOCK_SCALE_FACTOR)/(1024*100),
                  (mins_generic[16]/(double)FLINT_CLOCK_SCALE_FACTOR)/(1024*100));
          flint_printf("\t%.1lf|%.1lf|%.1lf",
                  (mins[17]/(double)FLINT_CLOCK_SCALE_FACTOR)/(65536*100),
                  (mins_shoup[17]/(double)FLINT_CLOCK_SCALE_FACTOR)/(65536*100),
                  (mins_generic[17]/(double)FLINT_CLOCK_SCALE_FACTOR)/(65536*100));
          flint_printf("\n");
      }
      else
      {
          flint_printf("%wd", i);
          for (int len = 1; len <= 16; ++len)
              flint_printf("\t%.1lf| na |%.1lf",
                      (mins[len-1]/(double)FLINT_CLOCK_SCALE_FACTOR)/(len*100),
                      (mins_generic[len-1]/(double)FLINT_CLOCK_SCALE_FACTOR)/(len*100));
          flint_printf("\t%.1lf| na |%.1lf",
                  (mins[16]/(double)FLINT_CLOCK_SCALE_FACTOR)/(1024*100),
                  (mins_generic[16]/(double)FLINT_CLOCK_SCALE_FACTOR)/(1024*100));
          flint_printf("\t%.1lf| na |%.1lf",
                  (mins[17]/(double)FLINT_CLOCK_SCALE_FACTOR)/(65536*100),
                  (mins_generic[17]/(double)FLINT_CLOCK_SCALE_FACTOR)/(65536*100));
          flint_printf("\n");
      }
   }

   return 0;
}
