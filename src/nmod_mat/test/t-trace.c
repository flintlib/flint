/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_mat.h"

TEST_FUNCTION_START(nmod_mat_trace, state)
{
    slong i;

    /* Test trace(AB) = trace(BA) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_mat_t A, B, AB, BA;
        mp_limb_t mod, trab, trba;
        slong m, n;

        mod = n_randtest_prime(state, 0);
        m = n_randint(state, 10);
        n = n_randint(state, 10);

        nmod_mat_init(A, m, n, mod);
        nmod_mat_init(B, n, m, mod);
        nmod_mat_init(AB, m, m, mod);
        nmod_mat_init(BA, n, n, mod);

        nmod_mat_randtest(A, state);
        nmod_mat_randtest(B, state);

        nmod_mat_mul(AB, A, B);
        nmod_mat_mul(BA, B, A);

        trab = nmod_mat_trace(AB);
        trba = nmod_mat_trace(BA);

        if (trab != trba)
            TEST_FUNCTION_FAIL(
                    "A = %{nmod_mat}\n"
                    "B = %{nmod_mat}\n"
                    "AB = %{nmod_mat}\n"
                    "BA = %{nmod_mat}\n"
                    "tr(AB) = %wu\n"
                    "tr(BA) = %wu\n",
                    A, B, AB, BA, trab, trba);

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(AB);
        nmod_mat_clear(BA);
    }

    TEST_FUNCTION_END(state);
}
