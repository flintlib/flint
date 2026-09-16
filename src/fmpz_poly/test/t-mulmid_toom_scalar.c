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

static int
mulmid_supported(slong len1, slong len2, slong nlo, slong nhi)
{
    slong sa = FLINT_MAX(0, nlo - (len2 - 1));
    slong ea = FLINT_MIN(len1 - 1, nhi - 1);
    slong sb = FLINT_MAX(0, nlo - (len1 - 1));
    slong eb = FLINT_MIN(len2 - 1, nhi - 1);
    slong m1 = ea - sa + 1;
    slong n1 = eb - sb + 1;
    slong w = nhi - nlo;
    slong nmax = FMPZ_POLY_TOOM_SCALAR_N_MAX;

    if (FLINT_MIN(m1, n1) == 1)
        return 1;

    return (m1 + n1 - 1 <= nmax) || (w + FLINT_MIN(m1, n1) - 1 <= nmax);
}

TEST_FUNCTION_START(fmpz_poly_mulmid_toom_scalar, state)
{
    slong len1, len2, nlo, nhi, i;
    int result;

    /* exhaustive sweep over (len1, len2, nlo, nhi) with small coefficients */
    for (len1 = 1; len1 <= 12; len1++)
    {
        for (len2 = 1; len2 <= 12; len2++)
        {
            slong prodlen = len1 + len2 - 1;

            for (nlo = 0; nlo < prodlen; nlo++)
            {
                for (nhi = nlo + 1; nhi <= prodlen; nhi++)
                {
                    fmpz *a, *b, *res, *ref;
                    slong w = nhi - nlo;

                    a = _fmpz_vec_init(len1);
                    b = _fmpz_vec_init(len2);
                    res = _fmpz_vec_init(w);
                    ref = _fmpz_vec_init(prodlen);

                    _fmpz_vec_randtest(a, state, len1, 8);
                    _fmpz_vec_randtest(b, state, len2, 8);

                    result = _fmpz_poly_mulmid_toom_scalar(res, a, len1,
                                 b, len2, nlo, nhi);

                    if (result != mulmid_supported(len1, len2, nlo, nhi))
                        TEST_FUNCTION_FAIL(
                            "unexpected return value %d\n"
                            "len1 = %wd, len2 = %wd, nlo = %wd, nhi = %wd\n",
                            result, len1, len2, nlo, nhi);

                    if (result)
                    {
                        _fmpz_poly_mul_classical(ref, a, len1, b, len2);

                        if (!_fmpz_vec_equal(res, ref + nlo, w))
                            TEST_FUNCTION_FAIL(
                                "len1 = %wd, len2 = %wd, "
                                "nlo = %wd, nhi = %wd\n"
                                "a = %{fmpz*}\n"
                                "b = %{fmpz*}\n"
                                "got = %{fmpz*}\n"
                                "expected = %{fmpz*}\n",
                                len1, len2, nlo, nhi, a, len1, b, len2,
                                res, w, ref + nlo, w);
                    }

                    _fmpz_vec_clear(a, len1);
                    _fmpz_vec_clear(b, len2);
                    _fmpz_vec_clear(res, w);
                    _fmpz_vec_clear(ref, prodlen);
                }
            }
        }
    }

    /* large coefficients, both argument orders, against fmpz_poly_mul */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz *a, *b, *r1, *r2, *ref;
        slong prodlen, w;
        int ok1, ok2;

        len1 = 1 + n_randint(state, 20);
        len2 = 1 + n_randint(state, 20);
        prodlen = len1 + len2 - 1;
        nlo = n_randint(state, prodlen);
        nhi = nlo + 1 + n_randint(state, prodlen - nlo);
        w = nhi - nlo;

        a = _fmpz_vec_init(len1);
        b = _fmpz_vec_init(len2);
        r1 = _fmpz_vec_init(w);
        r2 = _fmpz_vec_init(w);
        ref = _fmpz_vec_init(prodlen);

        _fmpz_vec_randtest(a, state, len1, 1 + n_randint(state, 2000));
        _fmpz_vec_randtest(b, state, len2, 1 + n_randint(state, 2000));

        ok1 = _fmpz_poly_mulmid_toom_scalar(r1, a, len1, b, len2, nlo, nhi);
        ok2 = _fmpz_poly_mulmid_toom_scalar(r2, b, len2, a, len1, nlo, nhi);

        if (ok1 != ok2 || (ok1 && !_fmpz_vec_equal(r1, r2, w)))
            TEST_FUNCTION_FAIL("argument order changed the result\n"
                               "len1 = %wd, len2 = %wd, nlo = %wd, nhi = %wd\n"
                               "a = %{fmpz*}\nb = %{fmpz*}\n",
                               len1, len2, nlo, nhi, a, len1, b, len2);

        if (ok1)
        {
            _fmpz_poly_mul_classical(ref, a, len1, b, len2);

            if (!_fmpz_vec_equal(r1, ref + nlo, w))
                TEST_FUNCTION_FAIL("len1 = %wd, len2 = %wd, nlo = %wd, "
                                   "nhi = %wd\n"
                                   "a = %{fmpz*}\nb = %{fmpz*}\n"
                                   "got = %{fmpz*}\nexpected = %{fmpz*}\n",
                                   len1, len2, nlo, nhi, a, len1, b, len2,
                                   r1, w, ref + nlo, w);
        }

        _fmpz_vec_clear(a, len1);
        _fmpz_vec_clear(b, len2);
        _fmpz_vec_clear(r1, w);
        _fmpz_vec_clear(r2, w);
        _fmpz_vec_clear(ref, prodlen);
    }

    /* fmpz_poly_t interface, including aliasing of the output */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c, d;
        slong prodlen;

        len1 = 1 + n_randint(state, 10);
        len2 = 1 + n_randint(state, 10);

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_init(d);

        fmpz_poly_randtest_not_zero(a, state, len1, 1 + n_randint(state, 500));
        fmpz_poly_randtest_not_zero(b, state, len2, 1 + n_randint(state, 500));

        prodlen = a->length + b->length - 1;
        nlo = n_randint(state, prodlen);
        nhi = nlo + 1 + n_randint(state, prodlen - nlo);

        result = fmpz_poly_mulmid_toom_scalar(c, a, b, nlo, nhi);

        if (result)
        {
            fmpz_poly_mul(d, a, b);
            fmpz_poly_shift_right(d, d, nlo);
            fmpz_poly_truncate(d, nhi - nlo);

            if (!fmpz_poly_equal(c, d))
                TEST_FUNCTION_FAIL("nlo = %wd, nhi = %wd\n"
                                   "a = %{fmpz_poly}\nb = %{fmpz_poly}\n"
                                   "got = %{fmpz_poly}\n"
                                   "expected = %{fmpz_poly}\n",
                                   nlo, nhi, a, b, c, d);

            fmpz_poly_set(c, a);
            fmpz_poly_mulmid_toom_scalar(c, c, b, nlo, nhi);

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

    TEST_FUNCTION_END(state);
}
