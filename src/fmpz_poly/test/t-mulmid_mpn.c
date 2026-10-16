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
#include "fmpz_poly.h"
#include "mpn_extras.h"
#include "fmpz_poly/impl.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(fmpz_poly_mulmid_mpn, state)
{
    slong iter;

    for (iter = 0; iter < 3000 * flint_test_multiplier(); iter++)
    {
        fmpz_poly_t a, b, c, d;
        slong nlo, nhi, bits, len1, len2, i;
        int aliasing, result, method, boundary;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_init(d);

        method = n_randint(state, 5) - 1;   /* -1 = automatic */
        _flint_mpn_poly_mulmid_force_method = method;
        _flint_mpn_poly_mulmid_force_cutoff = n_randint(state, 2) ? -1 : (slong) n_randint(state, 8);

        len1 = n_randint(state, n_randint(state, 4) ? 20 : 150);
        len2 = n_randint(state, n_randint(state, 4) ? 20 : 150);

        /* bit sizes close to limb boundaries are the most delicate */
        boundary = n_randint(state, 3);
        if (boundary == 0)
            bits = 1 + n_randint(state, 400);
        else
            bits = FLINT_BITS * (1 + n_randint(state, 5)) - 2 + n_randint(state, 4);

        fmpz_poly_randtest(a, state, len1, bits);
        fmpz_poly_randtest(b, state, len2, bits);
        fmpz_poly_randtest(c, state, len1 + len2, bits);
        fmpz_poly_randtest(d, state, len1 + len2, bits);

        /* random sign pattern: all nonnegative, all nonpositive, mixed */
        if (n_randint(state, 3) == 0)
        {
            for (i = 0; i < a->length; i++)
                fmpz_abs(a->coeffs + i, a->coeffs + i);
            if (n_randint(state, 2))
                for (i = 0; i < b->length; i++)
                    fmpz_abs(b->coeffs + i, b->coeffs + i);
        }
        else if (n_randint(state, 3) == 0)
        {
            /* extreme values */
            for (i = 0; i < a->length; i++)
                if (n_randint(state, 2))
                {
                    fmpz_one(a->coeffs + i);
                    fmpz_mul_2exp(a->coeffs + i, a->coeffs + i, bits - 1);
                    if (n_randint(state, 2)) fmpz_neg(a->coeffs + i, a->coeffs + i);
                    if (n_randint(state, 2)) fmpz_sub_ui(a->coeffs + i, a->coeffs + i, 1);
                }
            for (i = 0; i < b->length; i++)
                if (n_randint(state, 2))
                {
                    fmpz_one(b->coeffs + i);
                    fmpz_mul_2exp(b->coeffs + i, b->coeffs + i, bits - 1);
                    if (n_randint(state, 2)) fmpz_neg(b->coeffs + i, b->coeffs + i);
                    if (n_randint(state, 2)) fmpz_sub_ui(b->coeffs + i, b->coeffs + i, 1);
                }
        }
        _fmpz_poly_normalise(a);
        _fmpz_poly_normalise(b);

        nlo = n_randint(state, 40);
        nhi = n_randint(state, 200);
        aliasing = n_randint(state, 5);

        if (aliasing == 3 || aliasing == 4)
            fmpz_poly_set(b, a);

        fmpz_poly_mul_classical(c, a, b);
        fmpz_poly_shift_right(c, c, nlo);
        fmpz_poly_truncate(c, FLINT_MAX(0, nhi - nlo));

        if (aliasing == 0)
        {
            fmpz_poly_mulmid_mpn(d, a, b, nlo, nhi);
        }
        else if (aliasing == 1)
        {
            fmpz_poly_set(d, a);
            fmpz_poly_mulmid_mpn(d, d, b, nlo, nhi);
        }
        else if (aliasing == 2)
        {
            fmpz_poly_set(d, b);
            fmpz_poly_mulmid_mpn(d, a, d, nlo, nhi);
        }
        else if (aliasing == 3)
        {
            fmpz_poly_mulmid_mpn(d, a, a, nlo, nhi);
        }
        else if (aliasing == 4)
        {
            fmpz_poly_set(d, a);
            fmpz_poly_mulmid_mpn(d, d, d, nlo, nhi);
        }

        result = (fmpz_poly_equal(c, d));
        if (!result)
        {
            flint_printf("FAIL: fmpz_poly_mulmid_mpn\n");
            flint_printf("aliasing = %d, nlo = %wd, nhi = %wd, method = %d, bits = %wd\n", aliasing, nlo, nhi, method, bits);
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

    _flint_mpn_poly_mulmid_force_method = -1;
    _flint_mpn_poly_mulmid_force_cutoff = -1;

    TEST_FUNCTION_END(state);
}
