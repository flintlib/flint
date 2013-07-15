/* 

Copyright 2009, 2011 William Hart. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of
      conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list
      of conditions and the following disclaimer in the documentation and/or other materials
      provided with the distribution.

THIS SOFTWARE IS PROVIDED BY William Hart ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL William Hart OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those of the
authors and should not be interpreted as representing official policies, either expressed
or implied, of William Hart.

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
    
    flint_rand_t state;

    printf("mul_fft_main....");
    fflush(stdout);

    flint_randinit(state);
    _flint_rand_init_gmp(state);

    for (depth = 6; depth <= 13; depth++)
    {
        for (w = 1; w <= 3 - (depth >= 12); w++)
        {
            int iter = 1 + 200*(depth <= 8) + 80*(depth <= 9) + 10*(depth <= 10), i;
            
            for (i = 0; i < iter; i++)
            {
               mp_size_t n = (1UL<<depth);
               mp_bitcnt_t bits1 = (n*w - (depth + 1))/2; 
               mp_size_t len1 = 2*n + n_randint(state, 2*n) + 1;
               mp_size_t len2 = 2*n + 2 - len1 + n_randint(state, 2*n);

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
  
               mpn_mul(r2, i1, n1, i2, n2);
               flint_mpn_mul_fft_main(r1, i1, n1, i2, n2);
           
               for (j = 0; j < n1 + n2; j++)
               {
                   if (r1[j] != r2[j]) 
                   {
                       printf("error in limb %ld, %lx != %lx\n", j, r1[j], r2[j]);
                       abort();
                   }
               }

               flint_free(i1);
            }
        }
    }

    flint_randclear(state);
    
    printf("PASS\n");
    return 0;
}
