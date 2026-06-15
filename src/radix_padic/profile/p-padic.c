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
#include "radix.h"
#include "radix_padic.h"
#include "padic.h"
#include "nmod.h"

int main()
{
    flint_rand_t state;
    ulong p;
    slong n;
    int op;

    flint_rand_init(state);
    radix_t radix;

    p = 7;
    radix_init(radix, p, 0);

    flint_printf("precision %wu^n, limb radix = %wu^%wd\n\n", p, p, radix->exp);

    for (op = 0; op < 11; op++)
    {
        if (op == 0)
            flint_printf("x + y  (same valuation)\n");
        else if (op == 1)
            flint_printf("x + y  (delta val = 3)\n");
        else if (op == 2)
            flint_printf("x + y  (delta val = n / 2)\n");
        else if (op == 3)
            flint_printf("x * y\n");
        else if (op == 4)
            flint_printf("x * 100\n");
        else if (op == 5)
            flint_printf("inv(y)\n");
        else if (op == 6)
            flint_printf("inv(100)\n");
        else if (op == 7)
            flint_printf("x / y\n");
        else if (op == 8)
            flint_printf("x / 100\n");
        else if (op == 9)
            flint_printf("sqrt(y)\n");
        else if (op == 10)
            flint_printf("sqrt(100)\n");

        for (n = 10; n <= 11000000; n *= 2)
        {
            gr_ctx_t ctx;
            radix_integer_t x, y, z;
            radix_padic_t rx, ry, rz;
            padic_ctx_t pctx;
            padic_t px, py, pz;
            fmpz_t fp, pn, fx, fy, fz;

            fmpz_init_set_ui(fp, p);
            fmpz_init(pn);
            fmpz_ui_pow_ui(pn, p, n);

            padic_ctx_init(pctx, fp, n - 1, n + 1, PADIC_SERIES);
            gr_ctx_init_radix_padic(ctx, p, n, RADIX_PADIC_PREC_INF, 0);

            radix_integer_init(x, radix);
            radix_integer_init(y, radix);
            radix_integer_init(z, radix);

            padic_init2(px, n);
            padic_init2(py, n);
            padic_init2(pz, n);

            fmpz_init(fx);
            fmpz_init(fy);
            fmpz_init(fz);

            radix_padic_init(rx, ctx);
            radix_padic_init(ry, ctx);
            radix_padic_init(rz, ctx);

            do
            {
                fmpz_randm(fx, state, pn);
            } while (fmpz_divisible_ui(fx, p));
            do
            {
                fmpz_randm(fy, state, pn);
            } while (fmpz_divisible_ui(fy, p) || (op == 9 && fmpz_fdiv_ui(fy, p) != 1));


            if (op == 4 || op == 6 || op == 8 || op == 10)
                fmpz_set_ui(fy, 100);

            padic_set_fmpz(px, fx, pctx);
            padic_set_fmpz(py, fy, pctx);

            radix_padic_set_fmpz(rx, fx, ctx);
            radix_padic_set_fmpz(ry, fy, ctx);

            rx->N = n;
            ry->N = n;

            if (op == 1)
            {
                padic_shift(py, py, 3, pctx);
                ry->v = 3;
            }
            else if (op == 2)
            {
                padic_shift(py, py, n / 2, pctx);
                ry->v = n / 2;
            }

            double t1, t2, FLINT_SET_BUT_UNUSED(__);

            if (op == 0 || op == 1 || op == 2)
            {
                TIMEIT_START;
                padic_add(pz, px, py, pctx);
                TIMEIT_STOP_VALUES(__, t1);
                TIMEIT_START;
                radix_padic_add(rz, rx, ry, ctx);
                TIMEIT_STOP_VALUES(__, t2);
            }
            else if (op == 3 || op == 4)
            {
                TIMEIT_START;
                padic_mul(pz, px, py, pctx);
                TIMEIT_STOP_VALUES(__, t1);
                TIMEIT_START;
                radix_padic_mul(rz, rx, ry, ctx);
                TIMEIT_STOP_VALUES(__, t2);
            }
            else if (op == 5 || op == 6)
            {
                TIMEIT_START;
                padic_inv(pz, py, pctx);
                TIMEIT_STOP_VALUES(__, t1);
                TIMEIT_START;
                radix_padic_inv(rz, ry, ctx);
                TIMEIT_STOP_VALUES(__, t2);
            }
            else if (op == 7 || op == 8)
            {
                TIMEIT_START;
                padic_div(pz, px, py, pctx);
                TIMEIT_STOP_VALUES(__, t1);
                TIMEIT_START;
                radix_padic_div(rz, rx, ry, ctx);
                TIMEIT_STOP_VALUES(__, t2);
            }
            else if (op == 9 || op == 10)
            {
                TIMEIT_START;
                padic_sqrt(pz, py, pctx);
                TIMEIT_STOP_VALUES(__, t1);
                TIMEIT_START;
                radix_padic_sqrt(rz, ry, ctx);
                TIMEIT_STOP_VALUES(__, t2);
            }

            if (n == 10)
                flint_printf("       n      padic     radix_padic    speedup\n");

            flint_printf("%8wd   %10.2e    %10.2e     %5.2fx\n", n, t1, t2, t1 / t2);

            radix_padic_clear(rx, ctx);
            radix_padic_clear(ry, ctx);
            radix_padic_clear(rz, ctx);

            padic_clear(px);
            padic_clear(py);
            padic_clear(pz);

            radix_integer_clear(x, radix);
            radix_integer_clear(y, radix);
            radix_integer_clear(z, radix);

            fmpz_clear(fx);
            fmpz_clear(fy);
            fmpz_clear(fz);

            gr_ctx_clear(ctx);
            padic_ctx_clear(pctx);

            fmpz_clear(pn);
            fmpz_clear(fp);
        }

        flint_printf("\n");
    }

    radix_clear(radix);
    flint_rand_clear(state);
    flint_cleanup_master();
    return 0;
}

