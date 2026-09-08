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

TEST_FUNCTION_START(fmpz_div, state)
{
    slong i;

    for (i = 0; i < 5000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, q, q2, r;
        int exact, exact2, which;
        slong bits;

        fmpz_init(a); fmpz_init(b); fmpz_init(q); fmpz_init(q2); fmpz_init(r);

        bits = n_randint(state, 20) == 0 ? 100000 : 300;
        fmpz_randtest(a, state, bits);
        fmpz_randtest(b, state, n_randint(state, 3) == 0 ? 64 : bits);
        if (n_randint(state, 3) == 0)
            fmpz_mul(a, a, b);      /* exact case */

        /* reference: exact iff the truncated remainder vanishes */
        if (fmpz_is_zero(b))
        {
            exact2 = fmpz_is_zero(a);
            fmpz_zero(q2);
        }
        else
        {
            fmpz_tdiv_qr(q2, r, a, b);
            exact2 = fmpz_is_zero(r);
            if (!exact2)
                fmpz_zero(q2);
        }

        which = n_randint(state, 4);
        if (which == 0)
            exact = fmpz_div(q, a, b);
        else if (which == 1)
        {
            fmpz_set(q, a);
            exact = fmpz_div(q, q, b);
        }
        else if (which == 2 && fmpz_fits_si(b))
            exact = fmpz_div_si(q, a, fmpz_get_si(b));
        else if (which == 3 && fmpz_abs_fits_ui(b) && fmpz_sgn(b) >= 0)
            exact = fmpz_div_ui(q, a, fmpz_get_ui(b));
        else
        {
            fmpz_set(q, b);
            exact = fmpz_div(q, a, q);
        }

        if (exact != exact2 || !fmpz_equal(q, q2))
            TEST_FUNCTION_FAIL("which = %d, exact = %d, exact2 = %d\na = %{fmpz}\nb = %{fmpz}\nq = %{fmpz}\nq2 = %{fmpz}\n", which, exact, exact2, a, b, q, q2);

        fmpz_clear(a); fmpz_clear(b); fmpz_clear(q); fmpz_clear(q2); fmpz_clear(r);
    }

    TEST_FUNCTION_END(state);
}
