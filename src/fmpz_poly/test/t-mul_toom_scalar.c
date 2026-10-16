/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Opus 5

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

TEST_FUNCTION_START(fmpz_poly_mul_toom_scalar, state)
{
    slong len1, len2, i;
    int result;

    /* every supported shape, including zero-padded inputs */
    for (len1 = 1; len1 <= 19; len1++)
    {
        for (len2 = 1; len2 <= len1; len2++)
        {
            slong n = len1 + len2 - 1;

            if (n > (FMPZ_POLY_TOOM_SCALAR_N_MAX))
                continue;

            for (i = 0; i < 10 * flint_test_multiplier(); i++)
            {
                fmpz *a, *b, *res, *ref;
                flint_bitcnt_t bits = 1 + n_randint(state, 2000);

                a = _fmpz_vec_init(len1);
                b = _fmpz_vec_init(len2);
                res = _fmpz_vec_init(n);
                ref = _fmpz_vec_init(n);

                _fmpz_vec_randtest(a, state, len1, bits);
                _fmpz_vec_randtest(b, state, len2, bits);

                result = _fmpz_poly_mul_toom_scalar(res, a, len1, b, len2);
                _fmpz_poly_mul_classical(ref, a, len1, b, len2);

                if (!result || !_fmpz_vec_equal(res, ref, n))
                    TEST_FUNCTION_FAIL(
                        "len1 = %wd, len2 = %wd, returned %d\n"
                        "a = %{fmpz*}\n"
                        "b = %{fmpz*}\n"
                        "got = %{fmpz*}\n"
                        "expected = %{fmpz*}\n",
                        len1, len2, result, a, len1, b, len2,
                        res, n, ref, n);

                _fmpz_vec_clear(a, len1);
                _fmpz_vec_clear(b, len2);
                _fmpz_vec_clear(res, n);
                _fmpz_vec_clear(ref, n);
            }
        }
    }

    /* squaring */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz *a, *res, *ref;
        slong len = 1 + n_randint(state, 10);
        slong n = 2 * len - 1;

        if (n > (FMPZ_POLY_TOOM_SCALAR_N_MAX))
            continue;

        a = _fmpz_vec_init(len);
        res = _fmpz_vec_init(n);
        ref = _fmpz_vec_init(n);

        _fmpz_vec_randtest(a, state, len, 1 + n_randint(state, 1000));

        result = _fmpz_poly_mul_toom_scalar(res, a, len, a, len);
        _fmpz_poly_mul_classical(ref, a, len, a, len);

        if (!result || !_fmpz_vec_equal(res, ref, n))
            TEST_FUNCTION_FAIL("squaring, len = %wd\n"
                               "a = %{fmpz*}\n"
                               "got = %{fmpz*}\n"
                               "expected = %{fmpz*}\n",
                               len, a, len, res, n, ref, n);

        _fmpz_vec_clear(a, len);
        _fmpz_vec_clear(res, n);
        _fmpz_vec_clear(ref, n);
    }

    /* fmpz_poly_t interface, including aliasing of the output */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c, d;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_init(d);

        fmpz_poly_randtest(a, state, 1 + n_randint(state, 10),
                           1 + n_randint(state, 500));
        fmpz_poly_randtest(b, state, 1 + n_randint(state, 10),
                           1 + n_randint(state, 500));

        result = fmpz_poly_mul_toom_scalar(c, a, b);
        fmpz_poly_mul(d, a, b);

        if (result && !fmpz_poly_equal(c, d))
            TEST_FUNCTION_FAIL("a = %{fmpz_poly}\nb = %{fmpz_poly}\n"
                               "got = %{fmpz_poly}\n"
                               "expected = %{fmpz_poly}\n", a, b, c, d);

        if (result)
        {
            fmpz_poly_set(c, a);
            fmpz_poly_mul_toom_scalar(c, c, b);

            if (!fmpz_poly_equal(c, d))
                TEST_FUNCTION_FAIL("aliasing res = poly1\n"
                                   "a = %{fmpz_poly}\nb = %{fmpz_poly}\n"
                                   "got = %{fmpz_poly}\n"
                                   "expected = %{fmpz_poly}\n", a, b, c, d);
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
        fmpz_poly_clear(d);
    }

    /* shapes above the length bound must be refused */
    for (len1 = 1; len1 <= 25; len1++)
    {
        for (len2 = 1; len2 <= len1; len2++)
        {
            slong n = len1 + len2 - 1;
            fmpz *a, *b, *res;

            if (n <= (FMPZ_POLY_TOOM_SCALAR_N_MAX))
                continue;

            a = _fmpz_vec_init(len1);
            b = _fmpz_vec_init(len2);
            res = _fmpz_vec_init(n);

            _fmpz_vec_randtest(a, state, len1, 100);
            _fmpz_vec_randtest(b, state, len2, 100);

            if (_fmpz_poly_mul_toom_scalar(res, a, len1, b, len2) != 0)
                TEST_FUNCTION_FAIL("expected refusal for len1 = %wd, "
                                   "len2 = %wd\n", len1, len2);

            _fmpz_vec_clear(a, len1);
            _fmpz_vec_clear(b, len2);
            _fmpz_vec_clear(res, n);
        }
    }

    TEST_FUNCTION_END(state);
}
