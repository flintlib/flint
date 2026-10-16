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


TEST_FUNCTION_START(padic_radix_log_rectangular, state)
{
    slong iter;

    /* log(x y) == log(x) + log(y) on 1-units (built via exp for coverage) */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        padic_radix_t a, b, x, y, xy, lx, ly, lxly, lxy;
        int status = GR_SUCCESS;

        gr_ctx_init_padic_radix_randtest(ctx, state, 10);

        padic_radix_init(a, ctx);
        padic_radix_init(b, ctx);
        padic_radix_init(x, ctx);
        padic_radix_init(y, ctx);
        padic_radix_init(xy, ctx);
        padic_radix_init(lx, ctx);
        padic_radix_init(ly, ctx);
        padic_radix_init(lxly, ctx);
        padic_radix_init(lxy, ctx);

        GR_IGNORE(gr_randtest(a, state, ctx));
        GR_IGNORE(gr_randtest(b, state, ctx));

        status |= padic_radix_exp(x, a, ctx);
        status |= padic_radix_exp(y, b, ctx);
        status |= padic_radix_mul(xy, x, y, ctx);

        status |= padic_radix_log_rectangular(lx, x, ctx);
        status |= padic_radix_set(ly, y, ctx);      /* aliasing */
        status |= padic_radix_log_rectangular(ly, ly, ctx);
        status |= padic_radix_log_rectangular(lxy, xy, ctx);
        status |= padic_radix_add(lxly, lx, ly, ctx);

        if (status == GR_SUCCESS && padic_radix_equal(lxy, lxly, ctx) == T_FALSE)
        {
            flint_printf("FAIL: padic_radix_log_rectangular (homomorphism)\n");
            flint_printf("x = "); gr_println(x, ctx);
            flint_printf("y = "); gr_println(y, ctx);
            flint_printf("lxy = "); gr_println(lxy, ctx);
            flint_printf("lxly = "); gr_println(lxly, ctx);
            flint_abort();
        }

        padic_radix_clear(a, ctx);
        padic_radix_clear(b, ctx);
        padic_radix_clear(x, ctx);
        padic_radix_clear(y, ctx);
        padic_radix_clear(xy, ctx);
        padic_radix_clear(lx, ctx);
        padic_radix_clear(ly, ctx);
        padic_radix_clear(lxly, ctx);
        padic_radix_clear(lxy, ctx);

        gr_ctx_clear(ctx);
    }

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        ulong p;
        radix_t radix;
        slong thr, v, N;
        fmpz_t uf, pf, pv, yf, ref, got, pN;
        radix_integer_t y, g;

        radix_init_randtest_prime(radix, state);
        p = DIGIT_RADIX(radix);

        thr = (p == 2) ? 2 : 1;
        v = thr + (slong) n_randint(state, 6);
        N = 1 + (slong) n_randint(state, n_randint(state, 2) ? 60 : 200);

        if (v >= N)
        {
            radix_clear(radix);
            continue;
        }

        fmpz_init(uf);
        fmpz_init(pf);
        fmpz_init(pv);
        fmpz_init(yf);
        fmpz_init(ref);
        fmpz_init(got);
        fmpz_init(pN);
        fmpz_set_ui(pf, p);

        fmpz_randtest_unsigned(uf, state, 1 + n_randint(state, 300));
        if (fmpz_is_zero(uf))
            fmpz_one(uf);
        if (fmpz_divisible(uf, pf))
            fmpz_add_ui(uf, uf, 1);

        fmpz_pow_ui(pv, pf, v);
        fmpz_mul(yf, uf, pv);

        radix_integer_init(y, radix);
        radix_integer_init(g, radix);
        radix_integer_set_fmpz(y, yf, radix);

        _padic_radix_log_rectangular(g, y, N, radix);

        _padic_log(ref, yf, v, pf, N);
        fmpz_pow_ui(pN, pf, N);
        fmpz_mod(ref, ref, pN);

        radix_integer_get_fmpz(got, g, radix);

        if (!fmpz_equal(ref, got))
        {
            flint_printf("FAIL: _padic_radix_log_rectangular vs _padic_log\n");
            flint_printf("p = %wu, e = %wu, v = %wd, N = %wd\n", p, radix->exp, v, N);
            flint_printf("y = "); fmpz_print(yf); flint_printf("\n");
            flint_printf("ref = "); fmpz_print(ref); flint_printf("\n");
            flint_printf("got = "); fmpz_print(got); flint_printf("\n");
            flint_abort();
        }

        radix_integer_clear(y, radix);
        radix_integer_clear(g, radix);
        fmpz_clear(uf);
        fmpz_clear(pf);
        fmpz_clear(pv);
        fmpz_clear(yf);
        fmpz_clear(ref);
        fmpz_clear(got);
        fmpz_clear(pN);
        radix_clear(radix);
    }

    TEST_FUNCTION_END(state);
}
