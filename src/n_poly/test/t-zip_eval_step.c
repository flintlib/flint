/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod.h"
#include "n_poly.h"

/* a residue, often n - 1 */
static ulong
_zip_residue(flint_rand_t state, nmod_t mod)
{
    return n_randint(state, 4) == 0 ? mod.n - 1 : n_randint(state, mod.n);
}

/* _nmod_zip_eval_step and _n_fqp_zip_eval_step against nmod_mul and
   nmod_add, for lengths on both sides of the split into vector functions
   (from 16 terms for the first, up to 50 bits for the second) */
TEST_FUNCTION_START(zip_eval_step, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        nmod_t mod;
        slong i, j, length = n_randint(state, 50), d = 1 + n_randint(state, 5);
        nn_ptr cur, cur0, inc, coeffs, res;
        ulong t, s;

        nmod_init(&mod, n_randtest_not_zero(state));

        cur = flint_malloc((length + 1) * sizeof(ulong));
        cur0 = flint_malloc((length + 1) * sizeof(ulong));
        inc = flint_malloc((length + 1) * sizeof(ulong));
        coeffs = flint_malloc((d * length + 1) * sizeof(ulong));
        res = flint_malloc(d * sizeof(ulong));

        for (i = 0; i < length; i++)
        {
            cur0[i] = cur[i] = _zip_residue(state, mod);
            inc[i] = _zip_residue(state, mod);
        }
        for (i = 0; i < d * length; i++)
            coeffs[i] = _zip_residue(state, mod);

        t = _nmod_zip_eval_step(cur, inc, coeffs, length, mod);
        for (s = 0, i = 0; i < length; i++)
        {
            s = nmod_add(s, nmod_mul(cur0[i], coeffs[i], mod), mod);
            if (cur[i] != nmod_mul(cur0[i], inc[i], mod))
                TEST_FUNCTION_FAIL("nmod: n = %wu, length = %wd, cur[%wd]\n",
                                   mod.n, length, i);
        }
        if (t != s)
            TEST_FUNCTION_FAIL("nmod: n = %wu, length = %wd, sum\n", mod.n, length);

        for (i = 0; i < length; i++)
            cur[i] = cur0[i];
        _n_fqp_zip_eval_step(res, cur, inc, coeffs, length, d, mod);
        for (i = 0; i < length; i++)
            if (cur[i] != nmod_mul(cur0[i], inc[i], mod))
                TEST_FUNCTION_FAIL("fqp: n = %wu, length = %wd, d = %wd, cur[%wd]\n",
                                   mod.n, length, d, i);
        for (j = 0; j < d; j++)
        {
            for (s = 0, i = 0; i < length; i++)
                s = nmod_add(s, nmod_mul(cur0[i], coeffs[d * i + j], mod), mod);
            if (res[j] != s)
                TEST_FUNCTION_FAIL("fqp: n = %wu, length = %wd, d = %wd, res[%wd]\n",
                                   mod.n, length, d, j);
        }

        flint_free(cur);
        flint_free(cur0);
        flint_free(inc);
        flint_free(coeffs);
        flint_free(res);
    }

    TEST_FUNCTION_END(state);
}
