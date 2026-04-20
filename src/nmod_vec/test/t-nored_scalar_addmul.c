/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_vec.h"
#include "test_helpers.h"

static ulong
n_randtest_upto_bits(flint_rand_t state, ulong bits)
{
    ulong x = n_randtest(state);
    if (bits != FLINT_BITS)
        x &= ((UWORD(1) << bits) - 1);
    return x;
}

TEST_FUNCTION_START(nmod_vec_nored_scalar_addmul, state)
{
    slong iter;

    for (iter = 0; iter < 1000000 * flint_test_multiplier(); iter++)
    {
        nn_ptr Ain, A, B, Aref;
        nmod_t mod;
        ulong c;
        ulong n = n_randtest_not_zero(state);
        slong len = n_randint(state, 30);
        slong i;

        nmod_init(&mod, n);

        Ain = flint_malloc(sizeof(ulong) * 3 * len);
        A = flint_malloc(sizeof(ulong) * 3 * len);
        B = flint_malloc(sizeof(ulong) * len);
        Aref = flint_malloc(sizeof(ulong) * 3 * len);

        /* halflimb */

        for (i = 0; i < len; i++)
        {
            Ain[i] = A[i] = n_randtest_upto_bits(state, FLINT_BITS - 1);
            B[i] = n_randtest_upto_bits(state, FLINT_BITS / 2);
        }

        c = n_randtest_upto_bits(state, FLINT_BITS / 2);

        _nmod_vec_nored_scalar_addmul_halflimb(A, B, len, c);

        for (i = 0; i < len; i++)
            Aref[i] = Ain[i] + B[i] * c;

        if (!_nmod_vec_equal(A, Aref, len))
        {
            flint_printf("FAIL: _nmod_vec_nored_scalar_addmul_halflimb\n");
            flint_printf("Ain = %{ulong*}\n", Ain, len);
            flint_printf("A = %{ulong*}\n", A, len);
            flint_printf("B = %{ulong*}\n", B, len);
            flint_printf("c = %wu\n", c);
            flint_printf("Aref = %{ulong*}\n", Aref, len);
            flint_abort();
        }

        /* ll_halflimb */

        for (i = 0; i < len; i++)
        {
            Ain[2 * i + 1] = A[2 * i + 1] = n_randtest_upto_bits(state, FLINT_BITS - 1);
            Ain[2 * i] = A[2 * i] = n_randtest(state);
            B[i] = n_randtest_upto_bits(state, FLINT_BITS / 2);
        }

        c = n_randtest_upto_bits(state, FLINT_BITS / 2 - 1);

        _nmod_vec_nored_ll_scalar_addmul_halflimb(A, B, len, c);

        for (i = 0; i < len; i++)
        {
            ulong hi, lo;
            umul_ppmm(hi, lo, B[i], c);
            add_ssaaaa(Aref[2 * i + 1], Aref[2 * i], Ain[2 * i + 1], Ain[2 * i], hi, lo);
        }

        if (!_nmod_vec_equal(A, Aref, 2 * len))
        {
            flint_printf("FAIL: _nmod_vec_nored_ll_scalar_addmul_halflimb\n");
            flint_printf("Ain = %{ulong*}\n", Ain, 2 * len);
            flint_printf("A = %{ulong*}\n", A, 2 * len);
            flint_printf("B = %{ulong*}\n", B, len);
            flint_printf("c = %wu\n", c);
            flint_printf("Aref = %{ulong*}\n", Aref, 2 * len);
            flint_abort();
        }

        /* ll */

        for (i = 0; i < len; i++)
        {
            Ain[2 * i + 1] = A[2 * i + 1] = n_randtest_upto_bits(state, FLINT_BITS - 1);
            Ain[2 * i] = A[2 * i] = n_randtest(state);
            B[i] = n_randtest(state);
        }

        c = n_randtest_upto_bits(state, FLINT_BITS - 1);

        _nmod_vec_nored_ll_scalar_addmul(A, B, len, c);

        for (i = 0; i < len; i++)
        {
            ulong hi, lo;
            umul_ppmm(hi, lo, B[i], c);
            add_ssaaaa(Aref[2 * i + 1], Aref[2 * i], Ain[2 * i + 1], Ain[2 * i], hi, lo);
        }

        if (!_nmod_vec_equal(A, Aref, 2 * len))
        {
            flint_printf("FAIL: _nmod_vec_nored_ll_scalar_addmul\n");
            flint_printf("Ain = %{ulong*}\n", Ain, 2 * len);
            flint_printf("A = %{ulong*}\n", A, 2 * len);
            flint_printf("B = %{ulong*}\n", B, len);
            flint_printf("c = %wu\n", c);
            flint_printf("Aref = %{ulong*}\n", Aref, 2 * len);
            flint_abort();
        }

        /* lll */

        for (i = 0; i < len; i++)
        {
            Ain[3 * i + 2] = A[3 * i + 2] = n_randtest_upto_bits(state, FLINT_BITS - 1);
            Ain[3 * i + 1] = A[3 * i + 1] = n_randtest(state);
            Ain[3 * i] = A[3 * i] = n_randtest(state);
            B[i] = n_randtest(state);
        }

        c = n_randtest(state);

        _nmod_vec_nored_lll_scalar_addmul(A, B, len, c);

        for (i = 0; i < len; i++)
        {
            ulong hi, lo;
            umul_ppmm(hi, lo, B[i], c);
            add_sssaaaaaa(Aref[3 * i + 2], Aref[3 * i + 1], Aref[3 * i],
                          Ain[3 * i + 2], Ain[3 * i + 1], Ain[3 * i], 0, hi, lo);
        }

        if (!_nmod_vec_equal(A, Aref, 3 * len))
        {
            flint_printf("FAIL: _nmod_vec_nored_lll_scalar_addmul\n");
            flint_printf("Ain = %{ulong*}\n", Ain, 3 * len);
            flint_printf("A = %{ulong*}\n", A, 3 * len);
            flint_printf("B = %{ulong*}\n", B, len);
            flint_printf("c = %wu\n", c);
            flint_printf("Aref = %{ulong*}\n", Aref, 3 * len);
            flint_abort();
        }

        flint_free(Ain);
        flint_free(A);
        flint_free(B);
        flint_free(Aref);
    }

    TEST_FUNCTION_END(state);
}
