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

    printf("fft/ifft_radix2....");
    fflush(stdout);

    flint_randinit(state);
    _flint_rand_init_gmp(state);

    for (depth = 6; depth <= 12; depth++)
    {
        for (w = 1; w <= 5; w++)
        {
            mp_size_t n = (1UL<<depth);
            mp_size_t limbs = (n*w)/GMP_LIMB_BITS;
            mp_size_t size = limbs + 1;
            mp_size_t i;
            mp_limb_t * ptr;
            mp_limb_t ** ii, ** jj, *t1, *t2;
        
            ii = flint_malloc((2*(n + n*size) + 2*size)*sizeof(mp_limb_t));
            for (i = 0, ptr = (mp_limb_t *) ii + 2*n; i < 2*n; i++, ptr += size) 
            {
                ii[i] = ptr;
                random_fermat(ii[i], state, limbs);
            }
            t1 = ptr;
            t2 = t1 + size;
   
            for (i = 0; i < 2*n; i++)
               mpn_normmod_2expp1(ii[i], limbs);
    
            jj = flint_malloc(2*(n + n*size)*sizeof(mp_limb_t));
            for (i = 0, ptr = (mp_limb_t *) jj + 2*n; i < 2*n; i++, ptr += size) 
            {
                jj[i] = ptr;
                flint_mpn_copyi(jj[i], ii[i], size);
            }
   
            fft_radix2(ii, n, w, &t1, &t2);
            ifft_radix2(ii, n, w, &t1, &t2);
            for (i = 0; i < 2*n; i++)
            {
                mpn_div_2expmod_2expp1(ii[i], ii[i], limbs, depth + 1);
                mpn_normmod_2expp1(ii[i], limbs);
            }

            for (i = 0; i < 2*n; i++)
            {
                if (mpn_cmp(ii[i], jj[i], size) != 0)
                {
                    printf("FAIL:\n");
                    printf("Error in entry %ld\n", i);
                    abort();
                }
            }

            flint_free(ii);
            flint_free(jj);
        }
    }

    flint_randclear(state);
    
    printf("PASS\n");
    return 0;
}
