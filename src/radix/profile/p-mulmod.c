/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdint.h>
#include "fmpz.h"
#include "gr.h"
#include "profiler.h"
#include "gmp.h"
#include "radix.h"
#include "fmpz_mod.h"
#include "padic.h"
#include "nmod.h"

int main()
{
    flint_rand_t state;
    ulong p;
    slong n, exp;
    int have_mpn_mod;

    flint_rand_init(state);
    radix_t radix;

    p = 7;
    radix_init(radix, p, 0);
    exp = radix->exp;

    flint_printf("radix = %wu^%wd\n", p, radix->exp);

    for (n = 1; n <= 10000000; n *= 2)
    {
        fmpz_mod_ctx_t mod_ctx;
        gr_ctx_t gr_ctx;
        radix_integer_t x, y, z;
        gr_ptr gx = NULL, gy = NULL, gz = NULL;
        fmpz_t fx, fy, fz, fp, pn;
        padic_ctx_t pctx;
        padic_t px, py, pz;

        fmpz_init_set_ui(fp, p);
        fmpz_init(pn);
        fmpz_ui_pow_ui(pn, p, exp * n);

        fmpz_mod_ctx_init(mod_ctx, pn);
        padic_ctx_init(pctx, fp, exp * n, exp * n, PADIC_SERIES);
        have_mpn_mod = (gr_ctx_init_mpn_mod(gr_ctx, pn) == GR_SUCCESS);

        fmpz_init(fx);
        fmpz_init(fy);
        fmpz_init(fz);

        radix_integer_init(x, radix);
        radix_integer_init(y, radix);
        radix_integer_init(z, radix);

        padic_init2(px, exp * n);
        padic_init2(py, exp * n);
        padic_init2(pz, exp * n);

        do
        {
            fmpz_randm(fx, state, pn);
        } while (fmpz_divisible_ui(fx, p));
        do
        {
            fmpz_randm(fy, state, pn);
        } while (fmpz_divisible_ui(fy, p));

        padic_set_fmpz(px, fx, pctx);
        padic_set_fmpz(py, fy, pctx);

        if (have_mpn_mod)
        {
            gx = gr_heap_init(gr_ctx);
            gy = gr_heap_init(gr_ctx);
            gz = gr_heap_init(gr_ctx);

            GR_MUST_SUCCEED(gr_set_fmpz(gx, fx, gr_ctx));
            GR_MUST_SUCCEED(gr_set_fmpz(gy, fy, gr_ctx));
        }

        do
        {
            radix_integer_rand_limbs(x, state, n, radix);
        } while (x->d[0] % p == 0);
        do
        {
            radix_integer_rand_limbs(y, state, n, radix);
        } while (y->d[0] % p == 0);

        double t1, t2, t3, t4, FLINT_SET_BUT_UNUSED(__);

        TIMEIT_START;
        padic_mul(pz, px, py, pctx);
        TIMEIT_STOP_VALUES(__, t1);

        TIMEIT_START;
        fmpz_mod_mul(fz, fx, fy, mod_ctx);
        TIMEIT_STOP_VALUES(__, t2);

        if (have_mpn_mod)
        {
            TIMEIT_START;
            GR_IGNORE(gr_mul(gz, gx, gy, gr_ctx));
            TIMEIT_STOP_VALUES(__, t3);
        }
        else
            t3 = 0.0;

        TIMEIT_START;
        radix_integer_mullow_limbs(z, x, y, n, radix);
        TIMEIT_STOP_VALUES(__, t4);

        char st1[20];
        char st2[20];
        char st3[20];
        char st4[20];

        sprintf(st1, "%.2e", t1);
        sprintf(st2, "%.2e", t2);
        if (t3 == 0.0)
            sprintf(st3, "-");
        else
            sprintf(st3, "%.2e", t3);
        sprintf(st4, "%.2e", t4);

        if (n == 1)
            flint_printf("       n      padic   fmpz_mod    mpn_mod   radix_integer    speedup/padic  speedup/fmpz_mod\n");

        flint_printf("%8wd   %8s   %8s   %8s   %8s         %.2fx          %.2fx\n", n, st1, st2, st3, st4, t1 / t4, t2 / t4);

        padic_clear(px);
        padic_clear(py);
        padic_clear(pz);

        radix_integer_clear(x, radix);
        radix_integer_clear(y, radix);
        radix_integer_clear(z, radix);

        fmpz_clear(fx);
        fmpz_clear(fy);
        fmpz_clear(fz);

        if (have_mpn_mod)
        {
            gr_heap_clear(gx, gr_ctx);
            gr_heap_clear(gy, gr_ctx);
            gr_heap_clear(gz, gr_ctx);
            gr_ctx_clear(gr_ctx);
        }

        fmpz_mod_ctx_clear(mod_ctx);
        padic_ctx_clear(pctx);

        fmpz_clear(pn);
        fmpz_clear(fp);
    }

    radix_clear(radix);
    flint_rand_clear(state);
    flint_cleanup_master();
    return 0;
}

