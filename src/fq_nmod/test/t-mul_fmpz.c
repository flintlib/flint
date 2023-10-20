/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_poly.h"
#include "fmpz.h"
#include "fq_nmod.h"

TEST_FUNCTION_START(fq_nmod_mul_fmpz, state)
{
    int i, result;

    /* Check aliasing of a, b */
    for (i = 0; i < 2000; i++)
    {
        fmpz_t p;
        slong d;
        fq_nmod_ctx_t ctx;
        fmpz_t x;
        fq_nmod_t a, b;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        fq_nmod_ctx_init_conway(ctx, p, d, "a");

        fq_nmod_init(a, ctx);
        fq_nmod_init(b, ctx);
        fmpz_init(x);

        fq_nmod_randtest(a, state, ctx);
        fmpz_randtest_mod_signed(x,state,fq_nmod_ctx_prime(ctx));
        fq_nmod_mul_fmpz(b, a, x, ctx);
        fq_nmod_mul_fmpz(a, a, x, ctx);

        result = (fq_nmod_equal(a, b, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), fq_nmod_print_pretty(a, ctx), flint_printf("\n");
            flint_printf("b = "), fq_nmod_print_pretty(b, ctx), flint_printf("\n");
            flint_printf("x = "), fmpz_print(x), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fq_nmod_clear(a, ctx);
        fq_nmod_clear(b, ctx);
        fmpz_clear(x);
        fmpz_clear(p);
        fq_nmod_ctx_clear(ctx);
    }

    /* compare with direct multiplication */
    for (i = 0; i < 2000; i++)
    {
        fmpz_t p;
        slong d;
        fq_nmod_ctx_t ctx;
        fmpz_t x;
        fq_nmod_t a, c;
        nmod_poly_t b;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        fq_nmod_ctx_init_conway(ctx, p, d, "a");

        fq_nmod_init(a, ctx);
        fq_nmod_init(c, ctx);
        fmpz_init(x);
        nmod_poly_init(b, fmpz_get_ui(p));

        fq_nmod_randtest(a, state, ctx);
        fmpz_randtest_mod_signed(x, state, fq_nmod_ctx_prime(ctx));
        fq_nmod_mul_fmpz(c, a, x, ctx);

        if (fmpz_cmp_ui(x, 0) < 0)
        {
            nmod_poly_scalar_mul_nmod(b,a,-fmpz_get_si(x));
            nmod_poly_neg(b, b);
        }
        else
        {
            nmod_poly_scalar_mul_nmod(b,a,fmpz_get_ui(x));
        }

        result = (fq_nmod_equal(c, b, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), fq_nmod_print_pretty(a, ctx), flint_printf("\n");
            flint_printf("b = "), fq_nmod_print_pretty(b, ctx), flint_printf("\n");
            flint_printf("c = "), fq_nmod_print_pretty(b, ctx), flint_printf("\n");
            flint_printf("x = "), fmpz_print(x), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fq_nmod_clear(a, ctx);
        fq_nmod_clear(c, ctx);
        nmod_poly_clear(b);
        fmpz_clear(p);
        fmpz_clear(x);
        fq_nmod_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
