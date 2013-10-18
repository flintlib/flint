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

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fq_poly.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;
    flint_randinit(state);

    flint_printf("powmod_fmpz_binexp....");
    fflush(stdout);

    /* Aliasing of res and a */
    for (i = 0; i < 5 * flint_test_multiplier(); i++)
    {
        fq_ctx_t ctx;
        fq_poly_t a, res, t, f;
        ulong exp;
        fmpz_t expz;

        fq_ctx_randtest(ctx, state);

        exp = n_randint(state, 50);
        fmpz_init_set_ui(expz, exp);

        fq_poly_init(a, ctx);
        fq_poly_init(f, ctx);
        fq_poly_init(res, ctx);
        fq_poly_init(t, ctx);

        fq_poly_randtest(a, state, n_randint(state, 20), ctx);
        fq_poly_randtest_not_zero(f, state, n_randint(state, 20) + 1, ctx);

        fq_poly_powmod_fmpz_binexp(res, a, expz, f, ctx);
        fq_poly_powmod_fmpz_binexp(a, a, expz, f, ctx);

        result = (fq_poly_equal(res, a, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a:\n"); fq_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("f:\n"); fq_poly_print(f, ctx), flint_printf("\n\n");
            flint_printf("res:\n"); fq_poly_print(res, ctx), flint_printf("\n\n");
            abort();
        }

        fq_poly_clear(a, ctx);
        fq_poly_clear(f, ctx);
        fq_poly_clear(res, ctx);
        fq_poly_clear(t, ctx);
        fmpz_clear(expz);

        fq_ctx_clear(ctx);
    }

    /* Aliasing of res and f */
    for (i = 0; i < 5 * flint_test_multiplier(); i++)
    {
        fq_ctx_t ctx;
        fq_poly_t a, res, t, f;
        ulong exp;
        fmpz_t expz;

        fq_ctx_randtest(ctx, state);

        exp = n_randint(state, 50);
        fmpz_init_set_ui(expz, exp);

        fq_poly_init(a, ctx);
        fq_poly_init(f, ctx);
        fq_poly_init(res, ctx);
        fq_poly_init(t, ctx);

        fq_poly_randtest(a, state, n_randint(state, 20), ctx);
        fq_poly_randtest_not_zero(f, state, n_randint(state, 20) + 1, ctx);

        fq_poly_powmod_fmpz_binexp(res, a, expz, f, ctx);
        fq_poly_powmod_fmpz_binexp(f, a, expz, f, ctx);

        result = (fq_poly_equal(res, f, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a:\n"); fq_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("f:\n"); fq_poly_print(f, ctx), flint_printf("\n\n");
            flint_printf("res:\n"); fq_poly_print(res, ctx), flint_printf("\n\n");
            abort();
        }

        fq_poly_clear(a, ctx);
        fq_poly_clear(f, ctx);
        fq_poly_clear(res, ctx);
        fq_poly_clear(t, ctx);
        fmpz_clear(expz);

        fq_ctx_clear(ctx);
    }

    /* No aliasing */
    for (i = 0; i < 5 * flint_test_multiplier(); i++)
    {
        fq_ctx_t ctx;
        fq_poly_t a, res1, res2, t, f;
        ulong exp;
        fmpz_t expz;

        fq_ctx_randtest(ctx, state);

        exp = n_randint(state, 50);

        fq_poly_init(a, ctx);
        fq_poly_init(f, ctx);
        fq_poly_init(res1, ctx);
        fq_poly_init(res2, ctx);
        fq_poly_init(t, ctx);

        fq_poly_randtest(a, state, n_randint(state, 20), ctx);
        fq_poly_randtest_not_zero(f, state, n_randint(state, 20) + 1, ctx);
        fmpz_init_set_ui(expz, exp);

        fq_poly_powmod_fmpz_binexp(res1, a, expz, f, ctx);
        fq_poly_powmod_ui_binexp(res2, a, exp, f, ctx);

        result = (fq_poly_equal(res1, res2, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a:\n"); fq_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("f:\n"); fq_poly_print(f, ctx), flint_printf("\n\n");
            flint_printf("res1:\n"); fq_poly_print(res1, ctx), flint_printf("\n\n");
            flint_printf("res2:\n"); fq_poly_print(res2, ctx), flint_printf("\n\n");
            abort();
        }

        fq_poly_clear(a, ctx);
        fq_poly_clear(f, ctx);
        fq_poly_clear(res1, ctx);
        fq_poly_clear(res2, ctx);
        fq_poly_clear(t, ctx);
        fmpz_clear(expz);

        fq_ctx_clear(ctx);
    }

    /* Check that a^(b+c) = a^b * a^c */
    for (i = 0; i < 5 * flint_test_multiplier(); i++)
    {
        fq_ctx_t ctx;
        fq_poly_t a, res1, res2, res3, res4, t, f;
        fmpz_t exp1, exp2, exp3;

        fq_ctx_randtest(ctx, state);

        fmpz_init(exp1);
        fmpz_init(exp2);
        fmpz_randtest(exp1, state, 200);
        if (fmpz_sgn(exp1) == -1) fmpz_neg(exp1, exp1);
        fmpz_randtest(exp2, state, 200);
        if (fmpz_sgn(exp2) == -1) fmpz_neg(exp2, exp2);

        fq_poly_init(a, ctx);
        fq_poly_init(f, ctx);
        fq_poly_init(res1, ctx);
        fq_poly_init(res2, ctx);
        fq_poly_init(res3, ctx);
        fq_poly_init(res4, ctx);
        fq_poly_init(t, ctx);

        fq_poly_randtest(a, state, n_randint(state, 20), ctx);
        fq_poly_randtest_not_zero(f, state, n_randint(state, 20) + 1, ctx);

        fq_poly_powmod_fmpz_binexp(res1, a, exp1, f, ctx);
        fq_poly_powmod_fmpz_binexp(res2, a, exp2, f, ctx);
        fq_poly_mulmod(res4, res1, res2, f, ctx);
        fmpz_init(exp3);
        fmpz_add(exp3, exp1, exp2);
        fq_poly_powmod_fmpz_binexp(res3, a, exp3, f, ctx);

        result = (fq_poly_equal(res4, res3, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a:\n"); fq_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("f:\n"); fq_poly_print(f, ctx), flint_printf("\n\n");
            flint_printf("res3:\n"); fq_poly_print(res3, ctx), flint_printf("\n\n");
            flint_printf("res4:\n"); fq_poly_print(res4, ctx), flint_printf("\n\n");
            abort();
        }

        fq_poly_clear(a, ctx);
        fq_poly_clear(f, ctx);
        fq_poly_clear(res1, ctx);
        fq_poly_clear(res2, ctx);
        fq_poly_clear(res3, ctx);
        fq_poly_clear(res4, ctx);
        fq_poly_clear(t, ctx);
        fmpz_clear(exp1);
        fmpz_clear(exp2);
        fmpz_clear(exp3);

        fq_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    flint_printf("PASS\n");
    return 0;
}
