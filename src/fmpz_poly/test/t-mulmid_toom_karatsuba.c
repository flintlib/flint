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

TEST_FUNCTION_START(fmpz_poly_mulmid_toom_karatsuba, state)
{
    slong i;

    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        fmpz *a, *b, *res, *ref;
        slong len1, len2, prodlen, nlo, nhi, w;
        flint_bitcnt_t bits = 1 + n_randint(state, 400);
        int shape;

        if (n_randint(state, 2))
        {
            /* unbalanced shapes, including the middle product shape */
            len2 = 1 + n_randint(state, 40);
            len1 = len2 + n_randint(state, 3 * len2 + 1);
            if (n_randint(state, 3) == 0)
                len1 = 2 * len2 - 1;
        }
        else
        {
            len1 = 1 + n_randint(state, 80);
            len2 = 1 + n_randint(state, 80);
        }
        prodlen = len1 + len2 - 1;

        shape = n_randint(state, 5);
        if (shape == 0)         /* low */
        {
            nlo = 0;
            nhi = 1 + n_randint(state, prodlen);
        }
        else if (shape == 1)    /* high */
        {
            nlo = n_randint(state, prodlen);
            nhi = prodlen;
        }
        else if (shape == 2 && len1 >= len2)    /* middle */
        {
            nlo = len2 - 1;
            nhi = len1;
        }
        else
        {
            nlo = n_randint(state, prodlen);
            nhi = nlo + 1 + n_randint(state, prodlen - nlo);
        }
        w = nhi - nlo;

        a = _fmpz_vec_init(len1);
        b = _fmpz_vec_init(len2);
        res = _fmpz_vec_init(w);
        ref = _fmpz_vec_init(w);

        _fmpz_vec_randtest(a, state, len1, bits);
        _fmpz_vec_randtest(b, state, len2, bits);

        _fmpz_poly_mulmid_toom_karatsuba(res, a, len1, b, len2, nlo, nhi);
        _fmpz_poly_mulmid_classical(ref, a, len1, b, len2, nlo, nhi);

        if (!_fmpz_vec_equal(res, ref, w))
            TEST_FUNCTION_FAIL(
                "len1 = %wd, len2 = %wd, nlo = %wd, nhi = %wd\n"
                "a = %{fmpz*}\n"
                "b = %{fmpz*}\n"
                "got = %{fmpz*}\n"
                "expected = %{fmpz*}\n",
                len1, len2, nlo, nhi, a, len1, b, len2,
                res, w, ref, w);

        _fmpz_vec_clear(a, len1);
        _fmpz_vec_clear(b, len2);
        _fmpz_vec_clear(res, w);
        _fmpz_vec_clear(ref, w);
    }

    /* the fmpz_poly wrapper, including aliasing and out-of-range windows */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, res, ref;
        slong nlo, nhi;
        flint_bitcnt_t bits = 1 + n_randint(state, 200);
        int alias = n_randint(state, 3);

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(res);
        fmpz_poly_init(ref);

        fmpz_poly_randtest(a, state, n_randint(state, 50), bits);
        fmpz_poly_randtest(b, state, n_randint(state, 50), bits);

        nlo = n_randint(state, 100);
        nhi = n_randint(state, 100);

        fmpz_poly_mulmid(ref, a, b, nlo, nhi);

        if (alias == 1)
        {
            fmpz_poly_set(res, a);
            fmpz_poly_mulmid_toom_karatsuba(res, res, b, nlo, nhi);
        }
        else if (alias == 2)
        {
            fmpz_poly_set(res, b);
            fmpz_poly_mulmid_toom_karatsuba(res, a, res, nlo, nhi);
        }
        else
        {
            fmpz_poly_mulmid_toom_karatsuba(res, a, b, nlo, nhi);
        }

        if (!fmpz_poly_equal(res, ref))
            TEST_FUNCTION_FAIL(
                "alias = %d, nlo = %wd, nhi = %wd\n"
                "a = %{fmpz_poly}\n"
                "b = %{fmpz_poly}\n"
                "got = %{fmpz_poly}\n"
                "expected = %{fmpz_poly}\n",
                alias, nlo, nhi, a, b, res, ref);

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(res);
        fmpz_poly_clear(ref);
    }

    TEST_FUNCTION_END(state);
}
