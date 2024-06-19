/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2024 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gmpcompat.h"
#include "nmod.h"
#include "nmod_vec.h"

TEST_FUNCTION_START(_nmod_vec_dot_params, state)
{
    int i;

    for (i = 0; i < 100000 * flint_test_multiplier(); i++)
    {
        slong len;
        nmod_t mod;
        ulong m;
        int nlimbs1, nlimbs2;
        dot_params_t params;
        mpz_t t;

        len = n_randint(state, 10000) + 1;
        m = n_randtest_not_zero(state);

        nmod_init(&mod, m);

        params = _nmod_vec_dot_params(len, mod);
        if (params.method == _DOT0)
            nlimbs1 = 0;
        else if (params.method == _DOT1)
            nlimbs1 = 1;
        else if (params.method > _DOT2)
            nlimbs1 = 3;
        else
            nlimbs1 = 2;

        mpz_init2(t, 4*FLINT_BITS);
        flint_mpz_set_ui(t, m-1);
        mpz_mul(t, t, t);
        flint_mpz_mul_ui(t, t, len);
        nlimbs2 = mpz_size(t);

        if (params.method != _DOT_POW2 && nlimbs1 != nlimbs2)
            TEST_FUNCTION_FAIL(
                    "m = %wu\n"
                    "len = %wd\n"
                    "nlimbs1 = %d\n"
                    "nlimbs2 = %d\n"
                    "bound: %{mpz}\n",
                    m, len, nlimbs1, nlimbs2, t);

// DOT_SPLIT only defined for 64-bit config
#if (FLINT_BITS == 64)
        ulong pow2;
        NMOD_RED(pow2, UWORD(1)<<DOT_SPLIT_BITS, mod);
        if (params.method == _DOT2_SPLIT && params.pow2_precomp != pow2)
            TEST_FUNCTION_FAIL(
                    "m = %wu\n"
                    "len = %wd\n"
                    "pow2_precomp = %d\n"
                    "actual pow2 = %d\n"
                    "(method _DOT2_SPLIT)\n",
                    m, len, params.pow2_precomp, pow2, t);
#endif  // FLINT_BITS == 64

        mpz_clear(t);
    }

    TEST_FUNCTION_END(state);
}
