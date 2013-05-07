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
    mp_size_t iters, j;
    double truncation;

    flint_rand_t state;

    printf("mul_truncate_sqrt2....");
    fflush(stdout);

    flint_randinit(state);
    _flint_rand_init_gmp(state);

    depth = 13;
    w = 1;
    iters = 1;
    truncation = 1.0;

    {
       mp_size_t n = (1UL<<depth);
       mp_bitcnt_t bits1 = (n*w - (depth + 1))/2; 
       mp_bitcnt_t bits = 2*n*bits1;
       mp_size_t int_limbs = ((mp_size_t)(truncation*bits))/FLINT_BITS;
       mp_size_t j;
       mp_limb_t * i1, *i2, *r1, *r2;
        
       printf("limbs = %ld\n", int_limbs);
       
       i1 = flint_malloc(6*int_limbs*sizeof(mp_limb_t));
       i2 = i1 + int_limbs;
       r1 = i2 + int_limbs;
       r2 = r1 + 2*int_limbs;
   
       flint_mpn_urandomb(i1, state->gmp_state, int_limbs*FLINT_BITS);
       flint_mpn_urandomb(i2, state->gmp_state, int_limbs*FLINT_BITS);
  
       //mpn_mul(r2, i1, int_limbs, i2, int_limbs);
       for (j = 0; j < iters; j++)
          mul_truncate_sqrt2(r1, i1, int_limbs, i2, int_limbs, depth, w);

       flint_free(i1);
    }

    flint_randclear(state);
    
    printf("done\n");
    return 0;
}
