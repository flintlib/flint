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
    mp_size_t iters;
    
    flint_rand_t state;

    printf("mul_fft_main....");
    fflush(stdout);

    flint_randinit(state);
    _flint_rand_init_gmp(state);

    iters = 1;
    
    {
       mp_size_t int_limbs = 1000000;
       mp_size_t j;
       mp_limb_t * i1, *i2, *r1, *r2;
        
       printf("bits = %ld\n", int_limbs*FLINT_BITS);
       
       i1 = flint_malloc(6*int_limbs*sizeof(mp_limb_t));
       i2 = i1 + int_limbs;
       r1 = i2 + int_limbs;
       r2 = r1 + 2*int_limbs;
   
       flint_mpn_urandomb(i1, state->gmp_state, int_limbs*FLINT_BITS);
       flint_mpn_urandomb(i2, state->gmp_state, int_limbs*FLINT_BITS);
  
       for (j = 0; j < iters; j++)
          //mpn_mul(r2, i1, int_limbs, i2, int_limbs);
          flint_mpn_mul_fft_main(r1, i1, int_limbs, i2, int_limbs);

       flint_free(i1);
    }

    flint_randclear(state);
    
    printf("done\n");
    return 0;
}
