/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "longlong.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "radix.h"
#include "padic.h"
#include "padic_radix.h"
#include "gr.h"

TEST_FUNCTION_START(padic_radix_exp_rectangular, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        padic_radix_t x, y, xy, ex, ey, exey, exy;
        int status = GR_SUCCESS;

        gr_ctx_init_padic_radix_randtest(ctx, state, n_randint(state, 4) ? 10 : 100);

        padic_radix_init(x, ctx);
        padic_radix_init(y, ctx);
        padic_radix_init(xy, ctx);
        padic_radix_init(ex, ctx);
        padic_radix_init(ey, ctx);
        padic_radix_init(exy, ctx);
        padic_radix_init(exey, ctx);

        GR_IGNORE(gr_randtest(x, state, ctx));
        GR_IGNORE(gr_randtest(y, state, ctx));
        status |= padic_radix_add(xy, x, y, ctx);
        GR_IGNORE(gr_randtest(ex, state, ctx));

        status |= padic_radix_exp_rectangular(ex, x, ctx);
        status |= padic_radix_set(ey, y, ctx);  /* aliasing */
        status |= padic_radix_exp_rectangular(ey, ey, ctx);
        status |= padic_radix_exp_rectangular(exy, xy, ctx);
        status |= padic_radix_mul(exey, ex, ey, ctx);

        if (status == GR_SUCCESS && padic_radix_equal(exy, exey, ctx) == T_FALSE)
        {
            flint_printf("FAIL: padic_radix_exp_rectangular\n");
            flint_printf("x = "); gr_println(x, ctx);
            flint_printf("y = "); gr_println(y, ctx);
            flint_printf("xy = "); gr_println(xy, ctx);
            flint_printf("ex = "); gr_println(ex, ctx);
            flint_printf("ey = "); gr_println(ey, ctx);
            flint_printf("exy = "); gr_println(exy, ctx);
            flint_printf("exey = "); gr_println(exey, ctx);
            flint_abort();
        }

        padic_radix_clear(x, ctx);
        padic_radix_clear(y, ctx);
        padic_radix_clear(xy, ctx);
        padic_radix_clear(ex, ctx);
        padic_radix_clear(ey, ctx);
        padic_radix_clear(exy, ctx);
        padic_radix_clear(exey, ctx);

        gr_ctx_clear(ctx);
    }

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        ulong p;
        radix_t radix;
        slong thr, v, N;
        fmpz_t uf, pf, ref, got;
        radix_integer_t u, g;
        int status;

        radix_init_randtest_prime(radix, state);
        p = DIGIT_RADIX(radix);

        thr = (p == 2) ? 2 : 1;
        v = thr + (slong) n_randint(state, 6);
        N = 1 + (slong) n_randint(state, 80);

        fmpz_init(uf);
        fmpz_init(pf);
        fmpz_init(ref);
        fmpz_init(got);
        fmpz_set_ui(pf, p);

        fmpz_randtest_unsigned(uf, state, 1 + n_randint(state, 400));
        if (fmpz_is_zero(uf))
            fmpz_one(uf);
        if (fmpz_divisible(uf, pf))
            fmpz_add_ui(uf, uf, 1);

        radix_integer_init(u, radix);
        radix_integer_init(g, radix);
        radix_integer_set_fmpz(u, uf, radix);

        status = _padic_radix_exp_rectangular(g, u, v, N, radix);

        if (status == GR_SUCCESS)
        {
            _padic_exp(ref, uf, v, pf, N);
            radix_integer_get_fmpz(got, g, radix);

            if (!fmpz_equal(ref, got))
            {
                flint_printf("FAIL: _padic_radix_exp_balanced vs _padic_exp\n");
                flint_printf("p = %wu, e = %wu, v = %wd, N = %wd\n", p, radix->exp, v, N);
                flint_printf("u = "); fmpz_print(uf); flint_printf("\n");
                flint_printf("ref = "); fmpz_print(ref); flint_printf("\n");
                flint_printf("got = "); fmpz_print(got); flint_printf("\n");
                flint_abort();
            }
        }
        else if (status == GR_UNABLE)
        {
        }
        else
        {
            flint_printf("FAIL: unexpected status %d (p = %wu, e = %wu)\n", status, p, radix->exp);
            flint_abort();
        }

        radix_integer_clear(u, radix);
        radix_integer_clear(g, radix);
        fmpz_clear(uf);
        fmpz_clear(pf);
        fmpz_clear(ref);
        fmpz_clear(got);
        radix_clear(radix);
    }

    TEST_FUNCTION_END(state);
}
