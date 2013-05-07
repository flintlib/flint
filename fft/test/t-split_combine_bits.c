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
    int i;
    mp_size_t j;

    flint_rand_t state;

    printf("split/combine_bits....");
    fflush(stdout);

    flint_randinit(state);
    _flint_rand_init_gmp(state);

    for (i = 0; i < 10000; i++)
    {
        mp_size_t total_limbs = n_randint(state, 1000) + 1;
        mp_limb_t * in = flint_malloc(total_limbs*sizeof(mp_limb_t));
        mp_limb_t * out = flint_calloc(total_limbs, sizeof(mp_limb_t));
        
        mp_bitcnt_t bits = n_randint(state, 200) + 1;
        mp_size_t limbs = (2*bits - 1)/FLINT_BITS + 1;
        long length = (total_limbs*FLINT_BITS - 1)/bits + 1;
        
        mp_limb_t ** poly;
        poly = flint_malloc(length*sizeof(mp_limb_t *));
        for (j = 0; j < length; j++)
           poly[j] = flint_malloc((limbs + 1)*sizeof(mp_limb_t));

        flint_mpn_urandomb(in, state->gmp_state, total_limbs*FLINT_BITS);

        fft_split_bits(poly, in, total_limbs, bits, limbs);
        fft_combine_bits(out, poly, length, bits, limbs, total_limbs);
        
        for (j = 0; j < total_limbs; j++)
        {
           if (in[j] != out[j])
           {
              printf("FAIL:\n");
              printf("Error in limb %ld, %lu != %lu\n", j, in[j], out[j]);
              abort();
           }
        }

        flint_free(in);
        flint_free(out);

        for (j = 0; j < length; j++)
           flint_free(poly[j]);

        flint_free(poly);
    }

    flint_randclear(state);
    
    printf("PASS\n");
    return 0;
}
