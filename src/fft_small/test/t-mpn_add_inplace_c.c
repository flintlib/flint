/* 
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "mpn_extras.h"
#include "fft_small.h"

int main(void)
{
    FLINT_TEST_INIT(state);
    slong iter;

    flint_printf("mpn_add_inplace_c....");
    fflush(stdout);

    _flint_rand_init_gmp(state);

    for (iter = 0; iter < 1000; iter++)
    {
        mp_limb_t a[10], b[10], c[10];
        mp_size_t an, bn;
        unsigned char cf, c1, c2;

        bn = 1 + n_randint(state, 4);
        an = bn + n_randint(state, 4);

        flint_mpn_rrandom(a, state->gmp_state, an);
        flint_mpn_rrandom(b, state->gmp_state, bn);
        flint_mpn_copyi(c, a, an);
        cf = n_randint(state, 2);

        c1 = flint_mpn_add_inplace_c(a, an, b, bn, cf);

        c2 = mpn_add(c, c, an, b, bn);
        c2 += mpn_add_1(c, c, an, cf);

        if (c1 != c2 || mpn_cmp(a, c, an) != 0)
        {
            flint_printf("FAIL\n");
            flint_abort();
        }
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
