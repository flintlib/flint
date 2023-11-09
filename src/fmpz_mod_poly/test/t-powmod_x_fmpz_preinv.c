/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Martin Lee

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"

TEST_FUNCTION_START(fmpz_mod_poly_powmod_x_fmpz_preinv, state)
{
    int i, result;
    fmpz_mod_ctx_t ctx;

    fmpz_mod_ctx_init_ui(ctx, 2);

    /* No aliasing */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t a, res1, res2, f, finv;
        fmpz_t p;
        ulong exp;
        fmpz_t expz;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, p);
        exp = n_randint(state, 50);

        fmpz_mod_poly_init(f, ctx);
        fmpz_mod_poly_init(finv, ctx);
        fmpz_mod_poly_init(res1, ctx);
        fmpz_mod_poly_init(res2, ctx);

        fmpz_mod_poly_randtest_not_zero(f, state, n_randint(state, 50) + 1, ctx);
        fmpz_init_set_ui(expz, exp);

        fmpz_mod_poly_reverse (finv, f, f->length, ctx);
        fmpz_mod_poly_inv_series (finv, finv, f->length, ctx);

        fmpz_mod_poly_init2(a, f->length-1, ctx);
        fmpz_mod_poly_set_coeff_ui (a, 1, 1, ctx);
        fmpz_mod_poly_powmod_x_fmpz_preinv(res1, expz, f, finv, ctx);
        fmpz_mod_poly_powmod_fmpz_binexp_preinv(res2, a, expz, f, finv, ctx);

        result = (fmpz_mod_poly_equal(res1, res2, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("expz:\n"); fmpz_print(expz), flint_printf("\n\n");
            flint_printf("a:\n"); fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("f:\n"); fmpz_mod_poly_print(f, ctx), flint_printf("\n\n");
            flint_printf("finv:\n"); fmpz_mod_poly_print(finv, ctx), flint_printf("\n\n");
            flint_printf("res1:\n"); fmpz_mod_poly_print(res1, ctx), flint_printf("\n\n");
            flint_printf("res2:\n"); fmpz_mod_poly_print(res2, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(p);
        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(f, ctx);
        fmpz_mod_poly_clear(res1, ctx);
        fmpz_mod_poly_clear(res2, ctx);
        fmpz_mod_poly_clear(finv, ctx);
        fmpz_clear(expz);
    }

    /* Aliasing of res and f */
    for (i = 0; i < 25 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t res, f, finv;
        fmpz_t p;
        ulong exp;
        fmpz_t expz;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, p);
        exp = n_randint(state, 50);
        fmpz_init_set_ui(expz, exp);

        fmpz_mod_poly_init(f, ctx);
        fmpz_mod_poly_init(finv, ctx);
        fmpz_mod_poly_init(res, ctx);

        fmpz_mod_poly_randtest_not_zero(f, state, n_randint(state, 50) + 1, ctx);

        fmpz_mod_poly_reverse (finv, f, f->length, ctx);
        fmpz_mod_poly_inv_series (finv, finv, f->length, ctx);

        fmpz_mod_poly_powmod_x_fmpz_preinv(res, expz, f, finv, ctx);
        fmpz_mod_poly_powmod_x_fmpz_preinv(f, expz, f, finv, ctx);

        result = (fmpz_mod_poly_equal(res, f, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("f:\n"); fmpz_mod_poly_print(f, ctx), flint_printf("\n\n");
            flint_printf("finv:\n"); fmpz_mod_poly_print(finv, ctx), flint_printf("\n\n");
            flint_printf("expz:\n"); fmpz_print(expz), flint_printf("\n\n");
            flint_printf("res:\n"); fmpz_mod_poly_print(res, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(p);
        fmpz_mod_poly_clear(f, ctx);
        fmpz_mod_poly_clear(finv, ctx);
        fmpz_mod_poly_clear(res, ctx);
        fmpz_clear(expz);
    }

    /* Aliasing of res and finv */
    for (i = 0; i < 25 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t res, f, finv;
        fmpz_t p;
        ulong exp;
        fmpz_t expz;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, p);
        exp = n_randint(state, 50);
        fmpz_init_set_ui(expz, exp);

        fmpz_mod_poly_init(f, ctx);
        fmpz_mod_poly_init(finv, ctx);
        fmpz_mod_poly_init(res, ctx);

        fmpz_mod_poly_randtest_not_zero(f, state, n_randint(state, 50) + 1, ctx);

        fmpz_mod_poly_reverse (finv, f, f->length, ctx);
        fmpz_mod_poly_inv_series (finv, finv, f->length, ctx);

        fmpz_mod_poly_powmod_x_fmpz_preinv(res, expz, f, finv, ctx);
        fmpz_mod_poly_powmod_x_fmpz_preinv(finv, expz, f, finv, ctx);

        result = (fmpz_mod_poly_equal(res, finv, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("f:\n"); fmpz_mod_poly_print(f, ctx), flint_printf("\n\n");
            flint_printf("expz:\n"); fmpz_print(expz), flint_printf("\n\n");
            flint_printf("finv:\n"); fmpz_mod_poly_print(finv, ctx), flint_printf("\n\n");
            flint_printf("res:\n"); fmpz_mod_poly_print(res, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(p);
        fmpz_mod_poly_clear(f, ctx);
        fmpz_mod_poly_clear(finv, ctx);
        fmpz_mod_poly_clear(res, ctx);
        fmpz_clear(expz);
    }

    /* Check that x^(b+c) = x^b * x^c */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t res1, res2, res3, res4, f, finv;
        fmpz_t p;
        fmpz_t exp1, exp2, exp3;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_init(exp1);
        fmpz_init(exp2);
        fmpz_randtest(exp1, state, 200);
        if (fmpz_sgn(exp1) == -1) fmpz_neg(exp1, exp1);
        fmpz_randtest(exp2, state, 200);
        if (fmpz_sgn(exp2) == -1) fmpz_neg(exp2, exp2);

        fmpz_mod_poly_init(f, ctx);
        fmpz_mod_poly_init(finv, ctx);
        fmpz_mod_poly_init(res1, ctx);
        fmpz_mod_poly_init(res2, ctx);
        fmpz_mod_poly_init(res3, ctx);
        fmpz_mod_poly_init(res4, ctx);

        fmpz_mod_poly_randtest_not_zero(f, state, n_randint(state, 50) + 1, ctx);

        fmpz_mod_poly_reverse (finv, f, f->length, ctx);
        fmpz_mod_poly_inv_series (finv, finv, f->length, ctx);

        fmpz_mod_poly_powmod_x_fmpz_preinv(res1, exp1, f, finv, ctx);
        fmpz_mod_poly_powmod_x_fmpz_preinv(res2, exp2, f, finv, ctx);
        fmpz_mod_poly_mulmod(res4, res1, res2, f, ctx);
        fmpz_init(exp3);
        fmpz_add(exp3, exp1, exp2);
        fmpz_mod_poly_powmod_x_fmpz_preinv(res3, exp3, f, finv, ctx);

        result = (fmpz_mod_poly_equal(res4, res3, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("f:\n"); fmpz_mod_poly_print(f, ctx), flint_printf("\n\n");
            flint_printf("res3:\n"); fmpz_mod_poly_print(res3, ctx), flint_printf("\n\n");
            flint_printf("res4:\n"); fmpz_mod_poly_print(res4, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(p);
        fmpz_mod_poly_clear(f, ctx);
        fmpz_mod_poly_clear(finv, ctx);
        fmpz_mod_poly_clear(res1, ctx);
        fmpz_mod_poly_clear(res2, ctx);
        fmpz_mod_poly_clear(res3, ctx);
        fmpz_mod_poly_clear(res4, ctx);
        fmpz_clear(exp1);
        fmpz_clear(exp2);
        fmpz_clear(exp3);
    }

    fmpz_mod_ctx_clear(ctx);

    TEST_FUNCTION_END(state);
}
