/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2009 William Hart
    Copyright (C) 2013 Martin Lee

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#undef ulong
#define ulong ulongxx/* interferes with system includes */

#include <stdlib.h>
#include <stdio.h>

#undef ulong

#include <gmp.h>

#define ulong mp_limb_t

#include "flint.h"
#include "fmpz.h"
#include "fmpz_mod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    fmpz_mod_ctx_t ctx;
    FLINT_TEST_INIT(state);

    flint_printf("div_newton_n_preinv....");
    fflush(stdout);

    fmpz_mod_ctx_init_ui(ctx, 2);

    /* Compare to div_basecase */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, q, q2, binv;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(binv, ctx);
        fmpz_mod_poly_init(q, ctx);
        fmpz_mod_poly_init(q2, ctx);
        do
        fmpz_mod_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, ctx);
        while (b->length <= 2);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100), ctx);
        if (a->length > 2*(b->length)-3)
          fmpz_mod_poly_truncate (a, 2*(b->length)-3, ctx);

        {
            fmpz_t d;
            fmpz *leadB = fmpz_mod_poly_lead(b, ctx);

            fmpz_init(d);
            fmpz_gcd(d, p, leadB);
            while (!fmpz_is_one(d))
            {
                fmpz_divexact(leadB, leadB, d);
                fmpz_gcd(d, p, leadB);
            }
            fmpz_clear(d);
        }

        fmpz_mod_poly_reverse (binv, b, b->length, ctx);
        fmpz_mod_poly_inv_series_newton (binv, binv, b->length, ctx);
        fmpz_mod_poly_div_basecase(q, a, b, ctx);
        fmpz_mod_poly_div_newton_n_preinv(q2, a, b, binv, ctx);

        result = (fmpz_mod_poly_equal(q, q2, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("p = "), fmpz_print(p), flint_printf("\n\n");
            flint_printf("a = "), fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("b = "), fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            flint_printf("q = "), fmpz_mod_poly_print(q, ctx), flint_printf("\n\n");
            flint_printf("q2 = "), fmpz_mod_poly_print(q2, ctx), flint_printf("\n\n");
            flint_abort();
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(binv, ctx);
        fmpz_mod_poly_clear(q, ctx);
        fmpz_mod_poly_clear(q2, ctx);
        fmpz_clear(p);
    }

    /* Alias a and q */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, q, binv;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(binv, ctx);
        fmpz_mod_poly_init(q, ctx);
        do
        fmpz_mod_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, ctx);
        while (b->length <= 2);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100), ctx);
        if (a->length > 2*(b->length)-3)
          fmpz_mod_poly_truncate (a, 2*(b->length)-3, ctx);

        {
            fmpz_t d;
            fmpz *leadB = fmpz_mod_poly_lead(b, ctx);

            fmpz_init(d);
            fmpz_gcd(d, p, leadB);
            while (!fmpz_is_one(d))
            {
                fmpz_divexact(leadB, leadB, d);
                fmpz_gcd(d, p, leadB);
            }
            fmpz_clear(d);
        }

        fmpz_mod_poly_reverse (binv, b, b->length, ctx);
        fmpz_mod_poly_inv_series_newton (binv, binv, b->length, ctx);
        fmpz_mod_poly_div_newton_n_preinv(q, a, b, binv, ctx);
        fmpz_mod_poly_div_newton_n_preinv(a, a, b, binv, ctx);

        result = (fmpz_mod_poly_equal(q, a, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("p = "), fmpz_print(p), flint_printf("\n\n");
            flint_printf("a = "), fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("b = "), fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            flint_printf("q = "), fmpz_mod_poly_print(q, ctx), flint_printf("\n\n");
            flint_abort();
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(binv, ctx);
        fmpz_mod_poly_clear(q, ctx);
        fmpz_clear(p);
    }

    /* Alias b and q */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, q, binv;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(binv, ctx);
        fmpz_mod_poly_init(q, ctx);
        do
        fmpz_mod_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, ctx);
        while (b->length <= 2);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100), ctx);
        if (a->length > 2*(b->length)-3)
          fmpz_mod_poly_truncate (a, 2*(b->length)-3, ctx);

        {
            fmpz_t d;
            fmpz *leadB = fmpz_mod_poly_lead(b, ctx);

            fmpz_init(d);
            fmpz_gcd(d, p, leadB);
            while (!fmpz_is_one(d))
            {
                fmpz_divexact(leadB, leadB, d);
                fmpz_gcd(d, p, leadB);
            }
            fmpz_clear(d);
        }

        fmpz_mod_poly_reverse (binv, b, b->length, ctx);
        fmpz_mod_poly_inv_series_newton (binv, binv, b->length, ctx);
        fmpz_mod_poly_div_newton_n_preinv(q, a, b, binv, ctx);
        fmpz_mod_poly_div_newton_n_preinv(b, a, b, binv, ctx);

        result = (fmpz_mod_poly_equal(q, b, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("p = "), fmpz_print(p), flint_printf("\n\n");
            flint_printf("a = "), fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("b = "), fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            flint_printf("q = "), fmpz_mod_poly_print(q, ctx), flint_printf("\n\n");
            flint_abort();
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(binv, ctx);
        fmpz_mod_poly_clear(q, ctx);
        fmpz_clear(p);
    }


    /* Alias binv and q */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, q, binv;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(binv, ctx);
        fmpz_mod_poly_init(q, ctx);
        do
        fmpz_mod_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, ctx);
        while (b->length <= 2);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100), ctx);
        if (a->length > 2*(b->length)-3)
          fmpz_mod_poly_truncate (a, 2*(b->length)-3, ctx);

        {
            fmpz_t d;
            fmpz *leadB = fmpz_mod_poly_lead(b, ctx);

            fmpz_init(d);
            fmpz_gcd(d, p, leadB);
            while (!fmpz_is_one(d))
            {
                fmpz_divexact(leadB, leadB, d);
                fmpz_gcd(d, p, leadB);
            }
            fmpz_clear(d);
        }

        fmpz_mod_poly_reverse (binv, b, b->length, ctx);
        fmpz_mod_poly_inv_series_newton (binv, binv, b->length, ctx);
        fmpz_mod_poly_div_newton_n_preinv(q, a, b, binv, ctx);
        fmpz_mod_poly_div_newton_n_preinv(binv, a, b, binv, ctx);

        result = (fmpz_mod_poly_equal(q, binv, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("p = "), fmpz_print(p), flint_printf("\n\n");
            flint_printf("a = "), fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("b = "), fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            flint_printf("binv = "), fmpz_mod_poly_print(binv, ctx), flint_printf("\n\n");
            flint_printf("q = "), fmpz_mod_poly_print(q, ctx), flint_printf("\n\n");
            flint_abort();
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(binv, ctx);
        fmpz_mod_poly_clear(q, ctx);
        fmpz_clear(p);
    }

    fmpz_mod_ctx_clear(ctx);
    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
