/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mod_poly.h"
#include "fq.h"
#include "fq_poly.h"
#include "fq_embed.h"

#ifdef T
#undef T
#endif
#ifdef B
#undef B
#endif

#define T fq
#define CAP_T FQ
#define B fmpz_mod

#ifdef T
#ifdef B

#include "templates.h"

TEST_FUNCTION_START(fq_embed, state)
{
    int i;

    /* Check isomorphism to self */
    for (i = 0; i < 4 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) a, b, ev_a, ev_b;
        TEMPLATE(B, poly_t) minpoly;
        TEMPLATE(T, poly_t) minpoly_fq;

        TEMPLATE(T, ctx_randtest) (ctx, state);
        TEMPLATE(T, init) (a, ctx);
        TEMPLATE(T, init) (b, ctx);
        TEMPLATE(T, init) (ev_a, ctx);
        TEMPLATE(T, init) (ev_b, ctx);
        TEMPLATE(B, poly_init) (minpoly, ctx->ctxp);
        TEMPLATE(T, poly_init) (minpoly_fq, ctx);

        TEMPLATE(T, embed_gens) (a, b, minpoly, ctx, ctx);

        TEMPLATE4(T, poly_set, B, poly)(minpoly_fq, minpoly, ctx);
        TEMPLATE3(T, poly_evaluate, T)(ev_a, minpoly_fq, a, ctx);
        TEMPLATE3(T, poly_evaluate, T)(ev_b, minpoly_fq, b, ctx);

        if (!TEMPLATE(T, is_zero)(ev_a, ctx) || !TEMPLATE(T, is_zero)(ev_b, ctx))
        {
            flint_printf("FAIL:\n\n");
            flint_printf("CTX\n"), TEMPLATE(T, ctx_print)(ctx), flint_printf("\n");
            flint_printf("a: "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
            flint_printf("b: "), TEMPLATE(T, print_pretty)(b, ctx), flint_printf("\n");
            flint_printf("minpoly: "), TEMPLATE(B, poly_print_pretty)(minpoly, "x", ctx->ctxp), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, clear) (a, ctx);
        TEMPLATE(T, clear) (b, ctx);
        TEMPLATE(T, clear) (ev_a, ctx);
        TEMPLATE(T, clear) (ev_b, ctx);
        TEMPLATE(B, poly_clear) (minpoly, ctx->ctxp);
        TEMPLATE(T, poly_clear) (minpoly_fq, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
#endif

#undef B
#undef CAP_T
#undef T
