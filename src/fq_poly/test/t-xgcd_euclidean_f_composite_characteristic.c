/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "fq.h"
#include "fq_poly.h"

TEST_FUNCTION_START(fq_poly_xgcd_euclidean_f_composite_characteristic, state)
{
    slong i, j;

    for (i = 0; i < 10*flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t h;
        fq_ctx_t ctx;
        fq_poly_t a, b, g, s, t, t1, t2;
        fq_t f;
        fmpz_mod_ctx_t ctxp;

        fmpz_mod_ctx_init_ui(ctxp, 2 + n_randint(state, 200));
        fmpz_mod_poly_init(h, ctxp);

        fmpz_mod_poly_set_coeff_ui(h, 2, 1, ctxp);
        fmpz_mod_poly_set_coeff_ui(h, 1, 1, ctxp);
        fmpz_mod_poly_set_coeff_ui(h, 0, 1, ctxp);

        fq_ctx_init_modulus(ctx, h, ctxp, "t");
        fq_poly_init(a, ctx);
        fq_poly_init(b, ctx);
        fq_poly_init(g, ctx);
        fq_poly_init(s, ctx);
        fq_poly_init(t, ctx);
        fq_poly_init(t1, ctx);
        fq_poly_init(t2, ctx);
        fq_init(f, ctx);

        for (j = 0; j < 999; j++)
        {
            fq_poly_randtest(a, state, n_randint(state, 10), ctx);
            fq_poly_randtest(b, state, n_randint(state, 10), ctx);
            fq_poly_xgcd_euclidean_f(f, g, s, t, a, b, ctx);
            if (fq_is_one(f, ctx))
            {
                fq_poly_mul(t1, s, a, ctx);
                fq_poly_mul(t2, t, b, ctx);
                fq_poly_add(t1, t1, t2, ctx);
                if (!fq_poly_equal(t1, g, ctx))
                {
                    flint_printf("FAIL:\n");
                    flint_printf("p: "); fmpz_print(p); flint_printf("\n");
                    flint_printf("h: "); fmpz_mod_poly_print_pretty(h, "t", ctxp); flint_printf("\n");
                    flint_printf("f: "); fq_print_pretty(f, ctx); flint_printf("\n");
                    flint_printf("a: "); fq_poly_print_pretty(a, "x", ctx); flint_printf("\n");
                    flint_printf("b: "); fq_poly_print_pretty(b, "x", ctx); flint_printf("\n");
                    flint_printf("s: "); fq_poly_print_pretty(s, "x", ctx); flint_printf("\n");
                    flint_printf("t: "); fq_poly_print_pretty(t, "x", ctx); flint_printf("\n");
                    flint_printf("g: "); fq_poly_print_pretty(g, "x", ctx); flint_printf("\n");
                    flint_printf("s*a+t*b: "); fq_poly_print_pretty(t1, "x", ctx); flint_printf("\n");
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        fmpz_mod_poly_clear(h, ctxp);
        fmpz_mod_ctx_clear(ctxp);

        fq_poly_clear(a, ctx);
        fq_poly_clear(b, ctx);
        fq_poly_clear(g, ctx);
        fq_poly_clear(s, ctx);
        fq_poly_clear(t, ctx);
        fq_poly_clear(t1, ctx);
        fq_poly_clear(t2, ctx);
        fq_clear(f, ctx);
        fq_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
