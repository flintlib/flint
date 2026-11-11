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

/* a polynomial of length len (or less, after normalisation) */
static void
_addmul_poly(n_poly_t A, flint_rand_t state, slong len, nmod_t mod)
{
    slong i;

    n_poly_fit_length(A, len);
    for (i = 0; i < len; i++)
        A->coeffs[i] = n_randint(state, 4) == 0 ? mod.n - 1 : n_randint(state, mod.n);
    A->length = len;
    _n_poly_normalise(A);
}

/* checks A = B + d0 C coefficient by coefficient */
static int
_addmul_check(const n_poly_t A, const n_poly_t B, const n_poly_t C,
              ulong d0, nmod_t mod)
{
    slong i, len = FLINT_MAX(B->length, C->length);
    ulong t;

    if (A->length > len || (A->length > 0 && A->coeffs[A->length - 1] == 0))
        return 0;
    for (i = 0; i < len; i++)
    {
        t = nmod_mul(d0, i < C->length ? C->coeffs[i] : 0, mod);
        t = nmod_add(t, i < B->length ? B->coeffs[i] : 0, mod);
        if (t != (i < A->length ? A->coeffs[i] : 0))
            return 0;
    }
    return 1;
}

/* lengths on both sides of the switch to the vector functions (16) */
TEST_FUNCTION_START(n_poly_mod_scalar_addmul_nmod, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        n_poly_t A, B, C, B0, C0;
        nmod_t mod;
        ulong d0;
        int alias = n_randint(state, 4);

        nmod_init(&mod, n_randtest_not_zero(state));
        switch (n_randint(state, 4))
        {
            case 0: d0 = 0; break;
            case 1: d0 = 1 % mod.n; break;
            case 2: d0 = mod.n - 1; break;
            default: d0 = n_randint(state, mod.n);
        }

        n_poly_init(A);
        n_poly_init(B);
        n_poly_init(C);
        n_poly_init(B0);
        n_poly_init(C0);

        _addmul_poly(A, state, n_randint(state, 40), mod);
        _addmul_poly(B, state, n_randint(state, 40), mod);
        _addmul_poly(C, state, n_randint(state, 40), mod);
        n_poly_set(B0, B);
        n_poly_set(C0, C);

        if (alias == 0)
        {
            n_poly_mod_scalar_addmul_nmod(A, B, C, d0, mod);
            if (!_addmul_check(A, B0, C0, d0, mod))
                TEST_FUNCTION_FAIL("n = %wu, d0 = %wu, lengths %wd, %wd\n",
                                   mod.n, d0, B0->length, C0->length);
        }
        else if (alias == 1)
        {
            n_poly_mod_scalar_addmul_nmod(B, B, C, d0, mod);
            if (!_addmul_check(B, B0, C0, d0, mod))
                TEST_FUNCTION_FAIL("aliasing A = B: n = %wu, d0 = %wu, lengths %wd, %wd\n",
                                   mod.n, d0, B0->length, C0->length);
        }
        else if (alias == 2)
        {
            n_poly_mod_scalar_addmul_nmod(C, B, C, d0, mod);
            if (!_addmul_check(C, B0, C0, d0, mod))
                TEST_FUNCTION_FAIL("aliasing A = C: n = %wu, d0 = %wu, lengths %wd, %wd\n",
                                   mod.n, d0, B0->length, C0->length);
        }
        else
        {
            n_poly_mod_scalar_addmul_nmod(B, B, B, d0, mod);
            if (!_addmul_check(B, B0, B0, d0, mod))
                TEST_FUNCTION_FAIL("aliasing A = B = C: n = %wu, d0 = %wu, length %wd\n",
                                   mod.n, d0, B0->length);
        }

        n_poly_clear(A);
        n_poly_clear(B);
        n_poly_clear(C);
        n_poly_clear(B0);
        n_poly_clear(C0);
    }

    TEST_FUNCTION_END(state);
}
