/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod.h"
#include "nmod_vec.h"

TEST_FUNCTION_START(nmod_vec_discrete_log_pohlig_hellman, state)
{
    slong i, j, k;

    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        nmod_discrete_log_pohlig_hellman_t L;
        nmod_discrete_log_pohlig_hellman_init(L);
        for (j = 0; j < 10; j++)
        {
            double score;
            nmod_t fpctx;
            mp_limb_t p;

            p = n_randtest_prime(state, 1);
            nmod_init(&fpctx, p);
            score = nmod_discrete_log_pohlig_hellman_precompute_prime(L, p);
            if (score > 10000)
            {
                continue;
            }

            for (k = 0; k < 10; k++)
            {
                ulong x;
                mp_limb_t y, alpha = nmod_discrete_log_pohlig_hellman_primitive_root(L);

                x = n_urandint(state, p - 1);
                y = nmod_pow_ui(alpha, x, fpctx);
                if (x != nmod_discrete_log_pohlig_hellman_run(L, y))
                {

                    printf("FAIL\n");
                    flint_printf("modulo %wu log base %wu of %wu"
                                           " should be %wu\n", p, alpha, y, x);
                    fflush(stdout);
                    flint_abort();
                }
            }
        }
        nmod_discrete_log_pohlig_hellman_clear(L);
    }

    TEST_FUNCTION_END(state);
}
