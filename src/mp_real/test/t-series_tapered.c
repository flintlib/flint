/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"
#include "fmpz.h"
#include "arb.h"
#include "mp_real.h"
#include "mp_real/impl.h"

/* _mp_real_series_tapered: |res - f(t)| < 6 ulps for t < 2^-r and
   f = tan, atan, atanh, sin, 1 - cos over the range of the static
   tables (n <= nmax, r >= rmin), including t near 2^-r and tiny t */
TEST_FUNCTION_START(mp_real_series_tapered, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        int func = n_randint(state, 5);
        slong nmax = _mp_real_series_tapered_nmax(func);
        slong rmin = _mp_real_series_tapered_rmin(func);
        slong n = 1 + n_randint(state, (iter % 10 == 0) ? nmax : FLINT_MIN(nmax, 16));
        slong r = rmin + n_randint(state, (iter % 3 == 0) ? 400 : 64);
        slong b, j;
        nn_ptr x, y;
        arb_t xa, ya, t;
        fmpz_t f;
        double u;

        if (r >= FLINT_BITS * n)
            r = FLINT_BITS * n - 1;

        x = flint_malloc(n * sizeof(ulong));
        y = flint_malloc(n * sizeof(ulong));
        arb_init(xa);
        arb_init(ya);
        arb_init(t);
        fmpz_init(f);

        /* t < 2^-r: bits below position 64 n - r */
        flint_mpn_rrandom(x, state, n);
        b = FLINT_BITS * n - r;
        if (n_randint(state, 4) == 0)
            flint_mpn_store(x, n, ~UWORD(0));
        for (j = b / FLINT_BITS; j < n; j++)
            x[j] = (j == b / FLINT_BITS)
                ? (x[j] & ((UWORD(1) << (b % FLINT_BITS)) - 1)) : 0;
        if (n_randint(state, 8) == 0)
            mpn_rshift(x, x, n, 1 + n_randint(state, FLINT_BITS - 1));

        _mp_real_series_tapered(y, x, n, (flint_bitcnt_t) r, func);

        fmpz_set_ui_array(f, x, n);
        arb_set_fmpz(xa, f);
        arb_mul_2exp_si(xa, xa, -FLINT_BITS * n);
        if (func == MP_REAL_SERIES_TAN)
            arb_tan(t, xa, FLINT_BITS * n + 64);
        else if (func == MP_REAL_SERIES_ATAN)
            arb_atan(t, xa, FLINT_BITS * n + 64);
        else if (func == MP_REAL_SERIES_ATANH)
            arb_atanh(t, xa, FLINT_BITS * n + 64);
        else if (func == MP_REAL_SERIES_SIN)
            arb_sin(t, xa, FLINT_BITS * n + 64);
        else
        {
            arb_cos(t, xa, FLINT_BITS * n + 64);
            arb_sub_ui(t, t, 1, FLINT_BITS * n + 64);
            arb_neg(t, t);
        }
        fmpz_set_ui_array(f, y, n);
        arb_set_fmpz(ya, f);
        arb_mul_2exp_si(ya, ya, -FLINT_BITS * n);
        arb_sub(t, ya, t, FLINT_BITS * n + 64);
        arb_mul_2exp_si(t, t, FLINT_BITS * n);
        arb_abs(t, t);
        u = arf_get_d(arb_midref(t), ARF_RND_UP)
            + mag_get_d(arb_radref(t));

        if (u >= 6.0)
            TEST_FUNCTION_FAIL("func = %d, n = %wd, r = %wd, error = %f ulp\n",
                func, n, r, u);

        flint_free(x);
        flint_free(y);
        arb_clear(xa);
        arb_clear(ya);
        arb_clear(t);
        fmpz_clear(f);
    }

    TEST_FUNCTION_END(state);
}
