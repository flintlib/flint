/*
    Copyright (C) 2011 Fredrik Johansson

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
#include "nmod_vec.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);
    

    flint_printf("dot_bound_limbs....");
    fflush(stdout);

    for (i = 0; i < 100000 * flint_test_multiplier(); i++)
    {
        slong len;
        nmod_t mod;
        mp_limb_t m;
        int limbs1, limbs2;
        mpz_t t;

        len = n_randint(state, 10000) + 1;
        m = n_randtest_not_zero(state);

        nmod_init(&mod, m);

        limbs1 = _nmod_vec_dot_bound_limbs(len, mod);

        mpz_init2(t, 4*FLINT_BITS);
        flint_mpz_set_ui(t, m-1);
        mpz_mul(t, t, t);
        flint_mpz_mul_ui(t, t, len);
        limbs2 = mpz_size(t);

        if (limbs1 != limbs2)
        {
            flint_printf("FAIL:\n");
            flint_printf("m = %wu\n", m);
            flint_printf("len = %wd\n", len);
            flint_printf("limbs1 = %d\n", limbs1);
            flint_printf("limbs2 = %d\n", limbs2);
            gmp_printf("bound: %Zd\n", t);
            abort();
        }

        mpz_clear(t);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
