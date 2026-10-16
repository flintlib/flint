/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

TEST_FUNCTION_START(fmpz_poly_mulmid_classical_fft_small, state)
{
    slong iter;

    for (iter = 0; iter < 30 * flint_test_multiplier(); iter++)
    {
        slong len1 = 1 + n_randint(state, 12);
        slong len2 = 1 + n_randint(state, 12);
        slong big = FLINT_MAX(len1, len2);
        slong bits = 12000 + n_randint(state, 12000);
        slong nlo, nhi, i;
        int sq = (n_randint(state, 4) == 0);
        int ok;
        fmpz * a, * b, * r1, * r2;

        if (sq)
            len2 = len1;

        nlo = n_randint(state, len1 + len2 - 1);
        nhi = nlo + 1 + n_randint(state, len1 + len2 - 1 - nlo);

        a = _fmpz_vec_init(len1);
        b = sq ? a : _fmpz_vec_init(len2);
        r1 = _fmpz_vec_init(nhi - nlo);
        r2 = _fmpz_vec_init(nhi - nlo);

        for (i = 0; i < len1; i++)
        {
            fmpz_randbits(a + i, state, bits);
            if (n_randint(state, 2))
                fmpz_neg(a + i, a + i);
        }
        if (!sq)
            for (i = 0; i < len2; i++)
            {
                fmpz_randbits(b + i, state, bits);
                if (n_randint(state, 2))
                    fmpz_neg(b + i, b + i);
            }
        if (n_randint(state, 3) == 0)
            fmpz_zero(a + n_randint(state, len1));

        if (n_randint(state, 2))
            ok = _fmpz_poly_mulmid_classical_fft_small(r1, a, len1,
                                                       b, len2, nlo, nhi);
        else
            ok = _fmpz_poly_mulmid_classical_fft_small(r1, b, len2,
                                                       a, len1, nlo, nhi);

        if (!ok)
        {
            /* the density gate may refuse (e.g. a planted zero at
               small k): the fallback contract, not a failure */
            _fmpz_vec_clear(a, len1);
            if (!sq)
                _fmpz_vec_clear(b, len2);
            _fmpz_vec_clear(r1, nhi - nlo);
            _fmpz_vec_clear(r2, nhi - nlo);
            continue;
        }

        /* classical requires the longer operand first */
        if (len1 >= len2)
            _fmpz_poly_mulmid_classical(r2, a, len1, b, len2, nlo, nhi);
        else
            _fmpz_poly_mulmid_classical(r2, b, len2, a, len1, nlo, nhi);

        if (!_fmpz_vec_equal(r1, r2, nhi - nlo))
            TEST_FUNCTION_FAIL("len1 = %wd, len2 = %wd, bits = %wd, "
                "window [%wd, %wd), sq = %d\n",
                len1, len2, bits, nlo, nhi, sq);

        _fmpz_vec_clear(a, len1);
        if (!sq)
            _fmpz_vec_clear(b, len2);
        _fmpz_vec_clear(r1, nhi - nlo);
        _fmpz_vec_clear(r2, nhi - nlo);
        (void) big;
    }

    /* sparse and structured operands: coefficients drawn from
       {0, +1, -1, big} exercise the trivial-coefficient fast paths
       of the transformed ring (free conversion, pointwise ops
       degenerating to additions or integer bookkeeping) */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong len1 = 2 + n_randint(state, 10);
        slong len2 = 2 + n_randint(state, 10);
        slong bits = 12000 + n_randint(state, 8000);
        slong nlo, nhi, i;
        int ok;
        fmpz * a, * b, * r1, * r2;

        nlo = n_randint(state, len1 + len2 - 1);
        nhi = nlo + 1 + n_randint(state, len1 + len2 - 1 - nlo);

        a = _fmpz_vec_init(len1);
        b = _fmpz_vec_init(len2);
        r1 = _fmpz_vec_init(nhi - nlo);
        r2 = _fmpz_vec_init(nhi - nlo);

        for (i = 0; i < len1; i++)
            switch (n_randint(state, 4))
            {
                case 0: break;
                case 1: fmpz_set_si(a + i,
                            n_randint(state, 2) ? 1 : -1); break;
                default: fmpz_randtest(a + i, state, bits);
            }
        for (i = 0; i < len2; i++)
            switch (n_randint(state, 4))
            {
                case 0: break;
                case 1: fmpz_set_si(b + i,
                            n_randint(state, 2) ? 1 : -1); break;
                default: fmpz_randtest(b + i, state, bits);
            }
        /* one genuinely large coefficient per operand, so the
           transformed representation stays engaged (an all-trivial
           operand may legitimately refuse and fall back) */
        fmpz_randbits(a + n_randint(state, len1), state, bits);
        fmpz_randbits(b + n_randint(state, len2), state, bits);

        ok = _fmpz_poly_mulmid_classical_fft_small(r1, a, len1,
                                                   b, len2, nlo, nhi);
        if (!ok)
            goto sparse_next;   /* the density gate may refuse
                                   trivial-heavy inputs: the fallback
                                   contract, not a failure */
        if (len1 >= len2)
            _fmpz_poly_mulmid_classical(r2, a, len1, b, len2, nlo, nhi);
        else
            _fmpz_poly_mulmid_classical(r2, b, len2, a, len1, nlo, nhi);
        if (!_fmpz_vec_equal(r1, r2, nhi - nlo))
            TEST_FUNCTION_FAIL("sparse mismatch: len1 = %wd, "
                "len2 = %wd, bits = %wd, window [%wd, %wd)\n",
                len1, len2, bits, nlo, nhi);

sparse_next:
        _fmpz_vec_clear(a, len1);
        _fmpz_vec_clear(b, len2);
        _fmpz_vec_clear(r1, nhi - nlo);
        _fmpz_vec_clear(r2, nhi - nlo);
    }

    /* wrapper: aliasing, zero operands, empty windows; small
       coefficients exercise the refusal-and-fallback contract */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        fmpz_poly_t a, b, c, d;
        slong bits = n_randint(state, 2) ? 20000 : 100;
        slong len1, len2, nlo, nhi;
        int ok;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_init(d);
        fmpz_poly_randtest(a, state, 1 + n_randint(state, 10), bits);
        fmpz_poly_randtest(b, state, 1 + n_randint(state, 10), bits);
        len1 = a->length;
        len2 = b->length;

        if (len1 == 0 || len2 == 0)
        {
            nlo = 0;
            nhi = 1;
        }
        else
        {
            nlo = n_randint(state, len1 + len2 - 1);
            nhi = nlo + 1 + n_randint(state, len1 + len2 - 1 - nlo);
        }

        ok = fmpz_poly_mulmid_classical_fft_small(c, a, b, nlo, nhi);
        if (ok)
        {
            fmpz_poly_mulmid_classical(d,
                len1 >= len2 ? a : b, len1 >= len2 ? b : a, nlo, nhi);
            if (!fmpz_poly_equal(c, d))
                TEST_FUNCTION_FAIL("wrapper: bits = %wd, window "
                    "[%wd, %wd)\n", bits, nlo, nhi);

            /* aliased call must agree with the unaliased result */
            fmpz_poly_set(d, a);
            ok = fmpz_poly_mulmid_classical_fft_small(d, d, b, nlo, nhi);
            if (!ok || !fmpz_poly_equal(c, d))
                TEST_FUNCTION_FAIL("wrapper aliasing: bits = %wd\n",
                                   bits);
        }
        /* wrapper refusals fall back by contract (small bits or the
           density gate); nothing further to check */

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
        fmpz_poly_clear(d);
    }

    /* a zero operand yields a zero window and succeeds (deterministic
       coverage of the zero-operand branch) */
#if FLINT_HAVE_FFT_SMALL
    {
        fmpz * a = _fmpz_vec_init(5);
        fmpz * b = _fmpz_vec_init(3);
        fmpz * r = _fmpz_vec_init(4);
        slong i;
        int ok;
        for (i = 0; i < 5; i++)
            fmpz_randbits(a + i, state, 15000);
        fmpz_one(r + 0);   /* stale content must be overwritten */
        ok = _fmpz_poly_mulmid_classical_fft_small(r, a, 5, b, 3,
                                                   2, 6);
        if (!ok || !_fmpz_vec_is_zero(r, 4))
            TEST_FUNCTION_FAIL("zero operand window\n");
        _fmpz_vec_clear(a, 5);
        _fmpz_vec_clear(b, 3);
        _fmpz_vec_clear(r, 4);
    }
#endif

    TEST_FUNCTION_END(state);
}
