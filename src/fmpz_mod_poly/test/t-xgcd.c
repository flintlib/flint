/*
    Copyright (C) 2012 Sebastian Pancratz

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

TEST_FUNCTION_START(fmpz_mod_poly_xgcd, state)
{
    int i, result;
    fmpz_mod_ctx_t ctx;

    fmpz_mod_ctx_init_ui(ctx, 2);

    /* Generic case, most likely co-prime arguments ******************************/

    /* Compare with result from GCD and check correctness */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, d, g, s, t, v, w;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(d, ctx);
        fmpz_mod_poly_init(g, ctx);
        fmpz_mod_poly_init(s, ctx);
        fmpz_mod_poly_init(t, ctx);
        fmpz_mod_poly_init(v, ctx);
        fmpz_mod_poly_init(w, ctx);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100), ctx);
        fmpz_mod_poly_randtest(b, state, n_randint(state, 100), ctx);

        fmpz_mod_poly_gcd(d, a, b, ctx);
        fmpz_mod_poly_xgcd(g, s, t, a, b, ctx);

        fmpz_mod_poly_mul(v, s, a, ctx);
        fmpz_mod_poly_mul(w, t, b, ctx);
        fmpz_mod_poly_add(w, v, w, ctx);

        result = (fmpz_mod_poly_equal(d, g, ctx) && fmpz_mod_poly_equal(g, w, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(d, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(g, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(s, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(t, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(v, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(w, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(d, ctx);
        fmpz_mod_poly_clear(g, ctx);
        fmpz_mod_poly_clear(s, ctx);
        fmpz_mod_poly_clear(t, ctx);
        fmpz_mod_poly_clear(v, ctx);
        fmpz_mod_poly_clear(w, ctx);
        fmpz_clear(p);
    }

    /* Special case, arguments share a factor ********************************/

    /* Compare with result from GCD and check correctness */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, d, f, g, s, t, v, w;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(d, ctx);
        fmpz_mod_poly_init(f, ctx);
        fmpz_mod_poly_init(g, ctx);
        fmpz_mod_poly_init(s, ctx);
        fmpz_mod_poly_init(t, ctx);
        fmpz_mod_poly_init(v, ctx);
        fmpz_mod_poly_init(w, ctx);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100), ctx);
        fmpz_mod_poly_randtest(b, state, n_randint(state, 100), ctx);
        fmpz_mod_poly_randtest(f, state, n_randint(state, 20), ctx);
        fmpz_mod_poly_mul(a, a, f, ctx);
        fmpz_mod_poly_mul(b, b, f, ctx);

        fmpz_mod_poly_gcd(d, a, b, ctx);
        fmpz_mod_poly_xgcd(g, s, t, a, b, ctx);

        fmpz_mod_poly_mul(v, s, a, ctx);
        fmpz_mod_poly_mul(w, t, b, ctx);
        fmpz_mod_poly_add(w, v, w, ctx);

        result = (fmpz_mod_poly_equal(d, g, ctx) && fmpz_mod_poly_equal(g, w, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(d, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(f, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(g, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(s, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(t, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(v, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(w, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(d, ctx);
        fmpz_mod_poly_clear(f, ctx);
        fmpz_mod_poly_clear(g, ctx);
        fmpz_mod_poly_clear(s, ctx);
        fmpz_mod_poly_clear(t, ctx);
        fmpz_mod_poly_clear(v, ctx);
        fmpz_mod_poly_clear(w, ctx);
        fmpz_clear(p);
    }

    fmpz_mod_ctx_clear(ctx);

    TEST_FUNCTION_END(state);
}
