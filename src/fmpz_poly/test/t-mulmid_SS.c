/*
    Copyright (C) 2026 Fredrik Johansson

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
#include "ulong_extras.h"

TEST_FUNCTION_START(fmpz_poly_mulmid_SS, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        fmpz_poly_t a, b, c, d;
        slong nlo, nhi, bits;
        int aliasing, result;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_init(d);

        bits = 1 + n_randint(state, 200);
        fmpz_poly_randtest(a, state, n_randint(state, 10), bits);
        fmpz_poly_randtest(b, state, n_randint(state, 10), bits);
        fmpz_poly_randtest(c, state, n_randint(state, 10), bits);
        fmpz_poly_randtest(d, state, n_randint(state, 10), bits);

        nlo = n_randint(state, 10);
        nhi = n_randint(state, 10);
        aliasing = n_randint(state, 5);

        if (aliasing == 3 || aliasing == 4)
            fmpz_poly_set(b, a);

        fmpz_poly_mul_classical(c, a, b);
        fmpz_poly_shift_right(c, c, nlo);
        fmpz_poly_truncate(c, FLINT_MAX(0, nhi - nlo));

        if (aliasing == 0)
        {
            fmpz_poly_mulmid_SS(d, a, b, nlo, nhi);
        }
        else if (aliasing == 1)
        {
            fmpz_poly_set(d, a);
            fmpz_poly_mulmid_SS(d, d, b, nlo, nhi);
        }
        else if (aliasing == 2)
        {
            fmpz_poly_set(d, b);
            fmpz_poly_mulmid_SS(d, a, d, nlo, nhi);
        }
        else if (aliasing == 3)
        {
            fmpz_poly_mulmid_SS(d, a, a, nlo, nhi);
        }
        else if (aliasing == 4)
        {
            fmpz_poly_set(d, a);
            fmpz_poly_mulmid_SS(d, d, d, nlo, nhi);
        }

        result = (fmpz_poly_equal(c, d));
        if (!result)
        {
            flint_printf("FAIL: fmpz_poly_mulmid_SS\n");
            flint_printf("aliasing = %d, nlo = %wd, nhi = %wd\n", aliasing, nlo, nhi);
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(b), flint_printf("\n\n");
            fmpz_poly_print(c), flint_printf("\n\n");
            fmpz_poly_print(d), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
        fmpz_poly_clear(d);
    }

    /* longer operands and middle-product windows, which exercise the
       wraparound (cyclic) transforms at all depths including the
       matrix Fourier algorithm */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        fmpz_poly_t a, b, c, d;
        slong len1, len2, len, nlo, nhi, bits;
        int aliasing, result;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_init(d);

        bits = 1 + n_randint(state, 300);
        len1 = 1 + n_randint(state, n_randint(state, 2) ? 100 : 1200);
        len2 = 1 + n_randint(state, len1);

        fmpz_poly_randtest(a, state, len1, bits);
        fmpz_poly_randtest(b, state, len2, bits);
        aliasing = n_randint(state, 4);
        if (aliasing == 3)
            fmpz_poly_set(b, a);

        len = a->length + b->length - 1;
        if (a->length == 0 || b->length == 0)
            len = 0;

        switch (n_randint(state, 4))
        {
            case 0:
                /* the classical middle product (len1 - len2 + 1 terms) */
                nlo = FLINT_MIN(a->length, b->length) - 1;
                nhi = FLINT_MAX(a->length, b->length);
                break;
            case 1:
                /* a window around the middle */
                nlo = len / 3 + n_randint(state, len / 6 + 1);
                nhi = nlo + 1 + n_randint(state, len / 3 + 1);
                break;
            default:
                nlo = n_randint(state, len + 2);
                nhi = n_randint(state, len + 2);
                break;
        }

        nlo = FLINT_MAX(nlo, 0);
        nhi = FLINT_MAX(nhi, 0);

        fmpz_poly_mul(c, a, b);
        fmpz_poly_shift_right(c, c, nlo);
        fmpz_poly_truncate(c, FLINT_MAX(0, nhi - nlo));

        if (aliasing == 0 || aliasing == 3)
        {
            fmpz_poly_mulmid_SS(d, a, b, nlo, nhi);
        }
        else if (aliasing == 1)
        {
            fmpz_poly_set(d, a);
            fmpz_poly_mulmid_SS(d, d, b, nlo, nhi);
        }
        else
        {
            fmpz_poly_set(d, b);
            fmpz_poly_mulmid_SS(d, a, d, nlo, nhi);
        }

        result = (fmpz_poly_equal(c, d));
        if (!result)
        {
            flint_printf("FAIL: fmpz_poly_mulmid_SS (long)\n");
            flint_printf("aliasing = %d, len1 = %wd, len2 = %wd, nlo = %wd, nhi = %wd, bits = %wd\n",
                aliasing, a->length, b->length, nlo, nhi, bits);
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
        fmpz_poly_clear(d);
    }

    /* the same array as both operands at different lengths (a product
       with a truncation of itself, not a square) */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        fmpz * a, * r1, * r2;
        slong len1, len2, nlo, nhi, bits, i;

        len1 = 2 + n_randint(state, 200);
        len2 = 1 + n_randint(state, len1 - 1);
        bits = 1 + n_randint(state, 300);
        nlo = n_randint(state, len1 + len2 - 1);
        nhi = nlo + 1 + n_randint(state, len1 + len2 - 1 - nlo);

        a = _fmpz_vec_init(len1);
        r1 = _fmpz_vec_init(nhi - nlo);
        r2 = _fmpz_vec_init(nhi - nlo);
        for (i = 0; i < len1; i++)
            fmpz_randtest(a + i, state, bits);

        _fmpz_poly_mulmid_classical(r1, a, len1, a, len2, nlo, nhi);
        if (n_randint(state, 2))
            _fmpz_poly_mulmid_SS(r2, a, len1, a, len2, nlo, nhi);
        else
            _fmpz_poly_mulmid_SS(r2, a, len2, a, len1, nlo, nhi);

        if (!_fmpz_vec_equal(r1, r2, nhi - nlo))
        {
            flint_printf("FAIL: fmpz_poly_mulmid_SS (same array, different lengths)\n");
            flint_printf("len1 = %wd, len2 = %wd, nlo = %wd, nhi = %wd\n", len1, len2, nlo, nhi);
            fflush(stdout);
            flint_abort();
        }

        _fmpz_vec_clear(a, len1);
        _fmpz_vec_clear(r1, nhi - nlo);
        _fmpz_vec_clear(r2, nhi - nlo);
    }

    TEST_FUNCTION_END(state);
}
