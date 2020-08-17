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
    int i;
    mp_size_t j;

    FLINT_TEST_INIT(state);

    flint_printf("split/combine_bits....");
    fflush(stdout);

    
    _flint_rand_init_gmp(state);

    for (i = 0; i < 10000; i++)
    {
        mp_size_t total_limbs = n_randint(state, 1000) + 1;
        mp_limb_t * in = flint_malloc(total_limbs*sizeof(mp_limb_t));
        mp_limb_t * out = flint_calloc(total_limbs, sizeof(mp_limb_t));
        
        flint_bitcnt_t bits = n_randint(state, 200) + 1;
        mp_size_t limbs = (2*bits - 1)/FLINT_BITS + 1;
        slong length = (total_limbs*FLINT_BITS - 1)/bits + 1;
        
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
              flint_printf("FAIL:\n");
              flint_printf("Error in limb %wd, %wu != %wu\n", j, in[j], out[j]);
              abort();
           }
        }

        flint_free(in);
        flint_free(out);

        for (j = 0; j < length; j++)
           flint_free(poly[j]);

        flint_free(poly);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
