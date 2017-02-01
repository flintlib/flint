/* 
    Copyright (C) 2009, 2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fft.h"

int
main(void)
{
    mp_bitcnt_t depth, w;
    
    FLINT_TEST_INIT(state);

    flint_printf("mul_fft_main....");
    fflush(stdout);

    
    _flint_rand_init_gmp(state);

    for (depth = 11; depth <= 11; depth++)
    {
        w = 1;
        {
            mp_size_t n = (UWORD(1)<<depth);
            mp_bitcnt_t bits1 = (n*w - (depth + 1))/2; 
            mp_size_t len1 = 2*n + n_randint(state, 2*n) + 1;
            mp_size_t len2 = len1;

            int iter, i;
            
            for (i = 0; i < 1; i++)
            {
               mp_bitcnt_t b1 = len1*bits1, b2;
               mp_size_t n1, n2;
               mp_size_t j;
               mp_limb_t * i1, *i2, *r1, *r2;

               if (len2 <= 0)
                  len2 = 2*n + n_randint(state, 2*n) + 1;
               
               b2 = len2*bits1;
               
               n1 = (b1 - 1)/FLINT_BITS + 1;
               n2 = (b2 - 1)/FLINT_BITS + 1;
                    
               if (n1 < n2) /* ensure b1 >= b2 */
               {
                  mp_size_t t = n1;
                  mp_bitcnt_t tb = b1;
                  n1 = n2;
                  b1 = b2;
                  n2 = t;
                  b2 = tb;
               }

               i1 = flint_malloc(3*(n1 + n2)*sizeof(mp_limb_t));
               i2 = i1 + n1;
               r1 = i2 + n2;
               r2 = r1 + n1 + n2;
   
               flint_mpn_urandomb(i1, state->gmp_state, b1);
               flint_mpn_urandomb(i2, state->gmp_state, b2);
  
               /*mpn_mul(r2, i1, n1, i2, n2);*/
               if (i == 0) printf("n1 = %ld, n2 = %ld\n", n1, n2);

               flint_mpn_mul_fft_main(r1, i1, n1, i2, n2);
           
               /*for (j = 0; j < n1 + n2; j++)
               {
                   if (r1[j] != r2[j]) 
                   {
                       flint_printf("error in limb %wd, %wx != %wx\n", j, r1[j], r2[j]);
                       abort();
                   }
               }*/

               flint_free(i1);
            }
        }
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
