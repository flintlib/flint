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

TEST_FUNCTION_START(fmpz_poly_mul_toom_karatsuba, state)
{
    slong i;

    for (i = 0; i < 300 * flint_test_multiplier(); i++)
    {
        fmpz *a, *b, *res, *ref;
        slong len1, len2, n;
        flint_bitcnt_t bits = 1 + n_randint(state, 500);
        int squaring = n_randint(state, 2);

        /* both balanced and very unbalanced shapes, beyond the base case */
        len1 = 1 + n_randint(state, 90);
        if (n_randint(state, 2))
            len2 = squaring ? len1 : 1 + n_randint(state, 90);
        else
            len2 = len1;
        if (squaring)
            len2 = len1;
        n = len1 + len2 - 1;

        a = _fmpz_vec_init(len1);
        b = squaring ? a : _fmpz_vec_init(len2);
        res = _fmpz_vec_init(n);
        ref = _fmpz_vec_init(n);

        _fmpz_vec_randtest(a, state, len1, bits);
        if (!squaring)
            _fmpz_vec_randtest(b, state, len2, bits);

        _fmpz_poly_mul_toom_karatsuba(res, a, len1, b, len2);
        if (len1 >= len2)
            _fmpz_poly_mul_classical(ref, a, len1, b, len2);
        else
            _fmpz_poly_mul_classical(ref, b, len2, a, len1);

        if (!_fmpz_vec_equal(res, ref, n))
            TEST_FUNCTION_FAIL(
                "len1 = %wd, len2 = %wd, squaring = %d\n"
                "a = %{fmpz*}\n"
                "b = %{fmpz*}\n"
                "got = %{fmpz*}\n"
                "expected = %{fmpz*}\n",
                len1, len2, squaring, a, len1, b, len2,
                res, n, ref, n);

        _fmpz_vec_clear(a, len1);
        if (!squaring)
            _fmpz_vec_clear(b, len2);
        _fmpz_vec_clear(res, n);
        _fmpz_vec_clear(ref, n);
    }

    /* the fmpz_poly wrapper, including aliasing */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, res, ref;
        flint_bitcnt_t bits = 1 + n_randint(state, 300);
        int alias = n_randint(state, 3);

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(res);
        fmpz_poly_init(ref);

        fmpz_poly_randtest(a, state, n_randint(state, 60), bits);
        fmpz_poly_randtest(b, state, n_randint(state, 60), bits);

        fmpz_poly_mul_classical(ref, a, b);

        if (alias == 1)
        {
            fmpz_poly_set(res, a);
            fmpz_poly_mul_toom_karatsuba(res, res, b);
        }
        else if (alias == 2)
        {
            fmpz_poly_set(res, b);
            fmpz_poly_mul_toom_karatsuba(res, a, res);
        }
        else
        {
            fmpz_poly_mul_toom_karatsuba(res, a, b);
        }

        if (!fmpz_poly_equal(res, ref))
            TEST_FUNCTION_FAIL(
                "alias = %d\n"
                "a = %{fmpz_poly}\n"
                "b = %{fmpz_poly}\n"
                "got = %{fmpz_poly}\n"
                "expected = %{fmpz_poly}\n",
                alias, a, b, res, ref);

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(res);
        fmpz_poly_clear(ref);
    }

    TEST_FUNCTION_END(state);
}
