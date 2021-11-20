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
    

    flint_printf("dot....");
    fflush(stdout);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        slong len;
        nmod_t mod;
        mp_limb_t m, res;
        mp_ptr x, y;
        int limbs1;
        mpz_t s, t;
        slong j;

        len = n_randint(state, 1000) + 1;
        m = n_randtest_not_zero(state);

        nmod_init(&mod, m);

        x = _nmod_vec_init(len);
        y = _nmod_vec_init(len);

        _nmod_vec_randtest(x, state, len, mod);
        _nmod_vec_randtest(y, state, len, mod);

        limbs1 = _nmod_vec_dot_bound_limbs(len, mod);

        res = _nmod_vec_dot(x, y, len, mod, limbs1);

        mpz_init(s);
        mpz_init(t);

        for (j = 0; j < len; j++)
        {
            flint_mpz_set_ui(t, x[j]);
            flint_mpz_addmul_ui(s, t, y[j]);
        }

        flint_mpz_mod_ui(s, s, m);

        if (flint_mpz_get_ui(s) != res)
        {
            flint_printf("FAIL:\n");
            flint_printf("m = %wu\n", m);
            flint_printf("len = %wd\n", len);
            flint_printf("limbs1 = %d\n", limbs1);
            fflush(stdout);
            flint_abort();
        }

        mpz_clear(s);
        mpz_clear(t);

        _nmod_vec_clear(x);
        _nmod_vec_clear(y);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
