/*
    Copyright (C) 2013 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <math.h>
#include "profiler.h"
#include "flint.h"
#include "longlong.h"
#include "mpn_extras.h"
#include "ulong_extras.h"

#define mock_mulmod_preinvn(rxx, axx, bxx, nnn, nxx, ninv, norm)    \
   do {                                                             \
      mp_ptr __t;                                                   \
      TMP_INIT;                                                     \
                                                                    \
      TMP_START;                                                    \
      __t = TMP_ALLOC(3*(nnn)*sizeof(mp_limb_t));                   \
                                                                    \
      mpn_mul_n(__t, axx, bxx, nnn);                                \
      if (norm)                                                     \
         mpn_rshift(__t, __t, 2*(nnn), norm);                       \
                                                                    \
      mpn_tdiv_qr(__t + 2*(nnn), rxx, 0, __t, 2*(nnn), nxx, nnn);   \
      TMP_END;                                                      \
   } while (0)

typedef struct
{
   slong limbs;
   int algo;
} info_t;

void sample(void * arg, ulong count)
{
    info_t * info = (info_t *) arg;
    slong size = info->limbs, i, j;
    int algo = info->algo;
    int scale = 200;
   
    mpz_t a, b, d, r2;

    gmp_randstate_t st;
    FLINT_TEST_INIT(state);

    mp_ptr dinv;
    flint_bitcnt_t norm;
    
    mpz_init(a);
    mpz_init(b);
    mpz_init(d);
    /* don't init r2 */

    gmp_randinit_default(st);
    

    for (i = 0; i < count; i++)
    {
       mpz_rrandomb(a, st, size*FLINT_BITS);
       mpz_rrandomb(b, st, size*FLINT_BITS);
       do {
          mpz_rrandomb(d, st, size*FLINT_BITS);
       } while (mpz_sgn(d) == 0);
       
       /* reduce a, b mod d */
       mpz_fdiv_r(a, a, d);
       mpz_fdiv_r(b, b, d);

       /* normalise */
       count_leading_zeros(norm, d->_mp_d[d->_mp_size - 1]);
       mpz_mul_2exp(a, a, norm);
       mpz_mul_2exp(b, b, norm);
       mpz_mul_2exp(d, d, norm);

       dinv = flint_malloc(size*sizeof(mp_limb_t));
       flint_mpn_preinvn(dinv, d->_mp_d, size);

       r2->_mp_d = flint_malloc(size*sizeof(mp_limb_t));
       
       prof_start();
       if (algo == 1)
       {
          for (j = 0; j < scale; j++)
          {
             flint_mpn_mulmod_preinvn(r2->_mp_d, a->_mp_d, b->_mp_d, size, d->_mp_d, dinv, norm); 
          }
       } else
       {
          for (j = 0; j < scale; j++)
          {
             mock_mulmod_preinvn(r2->_mp_d, a->_mp_d, b->_mp_d, size, d->_mp_d, dinv, norm); 
          }
       }
       prof_stop();

       flint_free(r2->_mp_d);
       flint_free(dinv);
    }

    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(d);
    /* don't init r2 */

    gmp_randclear(st);
    flint_randclear(state);
}

int main(void)
{
   double min, max;
   info_t info;
   slong k, scale;

   printf("1: With precomputed inverse\n");
   printf("2: Without precomputed inverse\n\n");

   for (k = 1; k <= 10000; k = (slong) ceil(1.1*k))
   {
      info.limbs = k;
      info.algo = 1;

      scale = 200;
   
      prof_repeat(&min, &max, sample, (void *) &info);
         
      flint_printf("1: limbs %wd, min %.3g ms, max %.3g ms\n", 
           info.limbs,
		   ((min/(double)FLINT_CLOCK_SCALE_FACTOR)/scale)/2400000.0,
           ((max/(double)FLINT_CLOCK_SCALE_FACTOR)/scale)/2400000.0
	     );

     info.algo = 2;
     
     prof_repeat(&min, &max, sample, (void *) &info);
         
      flint_printf("2: limbs %wd, min %.3g ms, max %.3g ms\n\n", 
           info.limbs,
		   ((min/(double)FLINT_CLOCK_SCALE_FACTOR)/scale)/2400000.0,
           ((max/(double)FLINT_CLOCK_SCALE_FACTOR)/scale)/2400000.0
	     );
   }

   return 0;
}
