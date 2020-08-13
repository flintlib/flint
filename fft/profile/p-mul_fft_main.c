/* 
    Copyright (C) 2009, 2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
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
    
    FLINT_TEST_INIT(state);

    flint_printf("mul_fft_main....");
    fflush(stdout);

    
    _flint_rand_init_gmp(state);

    iters = 1;
    
    {
       mp_size_t int_limbs = 1000000;
       mp_size_t j;
       mp_limb_t * i1, *i2, *r1, *r2;
        
       flint_printf("bits = %wd\n", int_limbs*FLINT_BITS);
       
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
    
    flint_printf("done\n");
    return 0;
}
