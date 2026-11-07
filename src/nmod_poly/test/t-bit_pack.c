/*
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod.h"
#include "nmod_poly.h"
#include "fmpz.h"

TEST_FUNCTION_START(nmod_poly_bit_pack, state)
{
    int i, result;

    /* Check aliasing of a and c */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b;
        ulong n;
        ulong bits;
        nn_ptr mpn;
        slong nlo;

        do
            n = n_randtest_not_zero(state);
        while (n == 1);
        bits = 2 * FLINT_BIT_COUNT(n) + n_randint(state, FLINT_BITS);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        do
            nmod_poly_randtest(a, state, n_randint(state, 100));
        while (a->length == 0);

        mpn =
            flint_malloc(sizeof(ulong) *
                   ((bits * a->length - 1) / FLINT_BITS + 1));

        _nmod_poly_bit_pack(mpn, a->coeffs, a->length, bits);
        nmod_poly_fit_length(b, a->length);

        nlo = n_randint(state, a->length);
        _nmod_poly_bit_unpack(b->coeffs, 0, nlo, mpn, bits, a->mod);
        _nmod_poly_bit_unpack(b->coeffs + nlo, nlo, a->length, mpn, bits, a->mod);

        b->length = a->length;

        result = (nmod_poly_equal(a, b));
        if (!result)
            TEST_FUNCTION_FAIL(
                    "a = %{nmod_poly}\n"
                    "b = %{nmod_poly}\n",
                    a, b);

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        flint_free(mpn);
    }

    for (i = 0; i < 2000 * flint_test_multiplier(); i++)
    {
        fmpz_t f;
        nmod_poly_t A, B;
        slong b;
        ulong n;

        do
            n = n_randtest_not_zero(state);
        while (n == 1);

        fmpz_init(f);
        nmod_poly_init(A, n);
        nmod_poly_init(B, n);

        nmod_poly_randtest(A, state, 1+n_randint(state,100));

        b = FLINT_BIT_COUNT(n) + n_randint(state, FLINT_BITS);

        nmod_poly_bit_pack(f, A, b);
        nmod_poly_bit_unpack(B, f, b);

        if (!nmod_poly_equal(A, B))
            TEST_FUNCTION_FAIL(
                    "Input polynomial: %{nmod_poly}\n"
                    "Bit-packed integer: %{fmpz}\n"
                    "Output polynomial: %{nmod_poly}\n",
                    A, f, B);

        fmpz_clear(f);
        nmod_poly_clear(A);
        nmod_poly_clear(B);
    }

    /*
        _nmod_poly_bit_unpack with FLINT_BITS < bits < 2 FLINT_BITS on
        arbitrary fields (not packed residues), sometimes all ones; against
        n_ll_mod_preinv on the fields read with shifts
    */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_t mod;
        ulong n, ninv, bits;
        slong j, len, nlo, limbs;
        nn_ptr mpn, res;
        flint_bitcnt_t nu = 1 + n_randint(state, FLINT_BITS);

        n = n_randbits(state, nu);
        if (n == 1)
            n = 2;
        nmod_init(&mod, n);
        ninv = n_preinvert_limb(n);

        bits = FLINT_BITS + 1 + n_randint(state, FLINT_BITS - 1);
        len = 1 + n_randint(state, 40);
        nlo = n_randint(state, len);
        limbs = (len * bits - 1) / FLINT_BITS + 3;

        mpn = flint_malloc(sizeof(ulong) * limbs);
        res = flint_malloc(sizeof(ulong) * len);

        for (j = 0; j < limbs; j++)
            mpn[j] = n_randint(state, 4) ? n_randtest(state) : UWORD_MAX;

        _nmod_poly_bit_unpack(res, nlo, len, mpn, bits, mod);

        for (j = nlo; j < len; j++)
        {
            ulong lo, hi, r;
            ulong start = j * bits;
            slong w = start / FLINT_BITS;
            unsigned int sh = start % FLINT_BITS;

            /* the field of bits bits starting at bit start */
            lo = mpn[w] >> sh;
            hi = mpn[w + 1] >> sh;
            if (sh != 0)
            {
                lo |= mpn[w + 1] << (FLINT_BITS - sh);
                hi |= mpn[w + 2] << (FLINT_BITS - sh);
            }
            hi &= (UWORD(1) << (bits - FLINT_BITS)) - 1;
            r = n_ll_mod_preinv(hi, lo, n, ninv);

            if (res[j - nlo] != r)
                TEST_FUNCTION_FAIL("arbitrary fields\nn = %wu, bits = %wu, len = %wd, nlo = %wd, j = %wd\n"
                                   "hi = %wu, lo = %wu\ngot %wu, expected %wu\n",
                                   n, bits, len, nlo, j, hi, lo, res[j - nlo], r);
        }

        flint_free(mpn);
        flint_free(res);
    }

    TEST_FUNCTION_END(state);
}
