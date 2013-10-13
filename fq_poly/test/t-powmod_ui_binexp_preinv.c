/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2013 Mike Hansen

******************************************************************************/

#include "fq_poly.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;
    flint_randinit(state);

    flint_printf("powmod_ui_binexp_preinv....");
    fflush(stdout);

    /* Aliasing of res and a */
    for (i = 0; i < 50; i++)
    {
        fq_ctx_t ctx;
        fq_poly_t a, res1, t, f, finv;
        ulong exp;

        fq_ctx_randtest(ctx, state);
        
        exp = n_randint(state, 50);

        fq_poly_init(a);
        fq_poly_init(f);
        fq_poly_init(finv);
        fq_poly_init(res1);
        fq_poly_init(t);

        fq_poly_randtest(a, state, n_randint(state, 50), ctx);
        fq_poly_randtest_not_zero(f, state, n_randint(state, 50) + 1, ctx);

        fq_poly_reverse(finv, f, f->length, ctx);
        fq_poly_inv_series_newton(finv, finv, f->length, ctx);

        fq_poly_powmod_ui_binexp_preinv(res1, a, exp, f, finv, ctx);
        fq_poly_powmod_ui_binexp_preinv(a, a, exp, f, finv, ctx);

        result = (fq_poly_equal(res1, a));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("exp: %lu\n\n", exp);
            flint_printf("a:\n"); fq_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("f:\n"); fq_poly_print(f, ctx), flint_printf("\n\n");
            flint_printf("res:\n"); fq_poly_print(res1, ctx), flint_printf("\n\n");
            abort();
        }

        fq_poly_clear(a);
        fq_poly_clear(f);
        fq_poly_clear(finv);
        fq_poly_clear(res1);
        fq_poly_clear(t);

        fq_ctx_clear(ctx);
    }

    /* Aliasing of res and f */
    for (i = 0; i < 50; i++)
    {
        fq_ctx_t ctx;
        fq_poly_t a, res1, t, f, finv;
        ulong exp;

        fq_ctx_randtest(ctx, state);

        exp = n_randint(state, 50);

        fq_poly_init(a);
        fq_poly_init(f);
        fq_poly_init(finv);
        fq_poly_init(res1);
        fq_poly_init(t);

        fq_poly_randtest(a, state, n_randint(state, 50), ctx);
        fq_poly_randtest_not_zero(f, state, n_randint(state, 50) + 1, ctx);

        fq_poly_reverse(finv, f, f->length, ctx);
        fq_poly_inv_series_newton(finv, finv, f->length, ctx);

        fq_poly_powmod_ui_binexp_preinv(res1, a, exp, f, finv, ctx);
        fq_poly_powmod_ui_binexp_preinv(f, a, exp, f, finv, ctx);

        result = (fq_poly_equal(res1, f));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("exp: %lu\n\n", exp);
            flint_printf("a:\n"); fq_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("f:\n"); fq_poly_print(f, ctx), flint_printf("\n\n");
            flint_printf("res1:\n"); fq_poly_print(res1, ctx), flint_printf("\n\n");
            abort();
        }

        fq_poly_clear(a);
        fq_poly_clear(f);
        fq_poly_clear(finv);
        fq_poly_clear(res1);
        fq_poly_clear(t);

        fq_ctx_clear(ctx);
    }

    /* No aliasing */
    for (i = 0; i < 1000; i++)
    {
        fq_ctx_t ctx;
        fq_poly_t a, res1, res2, t, f, finv;
        ulong exp;
        int j;

        fq_ctx_randtest(ctx, state);

        exp = n_randint(state, 50);

        fq_poly_init(a);
        fq_poly_init(f);
        fq_poly_init(finv);
        fq_poly_init(res1);
        fq_poly_init(res2);
        fq_poly_init(t);

        fq_poly_randtest(a, state, n_randint(state, 50), ctx);
        fq_poly_randtest_not_zero(f, state, n_randint(state, 50) + 1, ctx);

        fq_poly_reverse(finv, f, f->length, ctx);
        fq_poly_inv_series_newton(finv, finv, f->length, ctx);

        fq_poly_powmod_ui_binexp_preinv(res1, a, exp, f, finv, ctx);

        fq_poly_zero(res2, ctx);
        fq_poly_one(res2, ctx);

        for (j = 1; j <= exp; j++)
            fq_poly_mulmod(res2, res2, a, f, ctx);

        result = (fq_poly_equal(res1, res2));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("exp: %lu\n\n", exp);
            flint_printf("a:\n"); fq_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("f:\n"); fq_poly_print(f, ctx), flint_printf("\n\n");
            flint_printf("res1:\n"); fq_poly_print(res1, ctx), flint_printf("\n\n");
            flint_printf("res2:\n"); fq_poly_print(res2, ctx), flint_printf("\n\n");
            abort();
        }

        fq_poly_clear(a);
        fq_poly_clear(f);
        fq_poly_clear(finv);
        fq_poly_clear(res1);
        fq_poly_clear(res2);
        fq_poly_clear(t);

        fq_ctx_clear(ctx);
    }

    /* Check that a^(b+c) = a^b * a^c */
    for (i = 0; i < 50; i++)
    {
        fq_ctx_t ctx;
        fq_poly_t a, res1, res2, res3, res4, t, f, finv;

        ulong exp1, exp2, exp3;

        fq_ctx_randtest(ctx, state);

        exp1 = n_randint(state, 50);
        exp2 = n_randint(state, 50);

        fq_poly_init(a);
        fq_poly_init(f);
        fq_poly_init(finv);
        fq_poly_init(res1);
        fq_poly_init(res2);
        fq_poly_init(res3);
        fq_poly_init(res4);
        fq_poly_init(t);

        fq_poly_randtest(a, state, n_randint(state, 50), ctx);
        fq_poly_randtest_not_zero(f, state, n_randint(state, 50) + 1, ctx);

        fq_poly_reverse(finv, f, f->length, ctx);
        fq_poly_inv_series_newton(finv, finv, f->length, ctx);


        fq_poly_powmod_ui_binexp_preinv(res1, a, exp1, f, finv, ctx);
        fq_poly_powmod_ui_binexp_preinv(res2, a, exp2, f, finv, ctx);
        
        fq_poly_mulmod_preinv(res4, res1, res2, f, finv, ctx);
        exp3 = exp1 + exp2;
        fq_poly_powmod_ui_binexp_preinv(res3, a, exp3, f, finv, ctx);

        result = (fq_poly_equal(res4, res3));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a:\n"); fq_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("f:\n"); fq_poly_print(f, ctx), flint_printf("\n\n");
            flint_printf("res3:\n"); fq_poly_print(res3, ctx), flint_printf("\n\n");
            flint_printf("res4:\n"); fq_poly_print(res4, ctx), flint_printf("\n\n");
            abort();
        }

        fq_poly_clear(a);
        fq_poly_clear(f);
        fq_poly_clear(finv);
        fq_poly_clear(res1);
        fq_poly_clear(res2);
        fq_poly_clear(res3);
        fq_poly_clear(res4);
        fq_poly_clear(t);

        fq_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    flint_printf("PASS\n");
    return 0;
}
