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
#include "fmpz_mod_mat.h"
#include "fq.h"
#include "fq_embed.h"

TEST_FUNCTION_START(fq_embed_matrices, state)
{
    slong ix, jx;

    /* Check that isomorphism to self gives identity matrices */
    for (ix = 0; ix < 20 * flint_test_multiplier(); ix++)
    {
        fq_ctx_t ctx;
        fq_t gen;
        const fmpz_mod_poly_struct *modulus;
        fmpz_mod_mat_t embed, project, one;
        slong d;

        fq_ctx_init_randtest(ctx, state, 3);
        d = fq_ctx_degree(ctx);
        modulus = fq_ctx_modulus(ctx);

        fq_init(gen, ctx);
        fq_gen(gen, ctx);
        fq_pow(gen, gen, fq_ctx_prime(ctx), ctx);

        fmpz_mod_mat_init(embed, d, d, ctx->ctxp);
        fmpz_mod_mat_init(project, d, d, ctx->ctxp);
        fmpz_mod_mat_init(one, d, d, ctx->ctxp);

        fq_embed_matrices(embed, project, gen, ctx, gen, ctx, modulus);
        fmpz_mod_mat_one(one, ctx->ctxp);

        if (!fmpz_mod_mat_equal(embed, one, ctx->ctxp) || !fmpz_mod_mat_equal(project, one, ctx->ctxp)) {
            flint_printf("FAIL:\n\n");
            flint_printf("CTX\n"), fq_ctx_print(ctx), flint_printf("\n");
            flint_printf("x^p: "), fq_print_pretty(gen, ctx), flint_printf("\n");
            flint_printf("Embed\n"),
                fmpz_mod_mat_print_pretty(embed, ctx->ctxp), flint_printf("\nProject\n"),
                fmpz_mod_mat_print_pretty(project, ctx->ctxp), flint_printf("\n");
            flint_abort();
        }

        fmpz_mod_mat_clear(embed, ctx->ctxp);
        fmpz_mod_mat_clear(project, ctx->ctxp);
        fmpz_mod_mat_clear(one, ctx->ctxp);
        fq_clear(gen, ctx);
        fq_ctx_clear(ctx);
    }

    /* Check random emebedding (degrees 1..5) */
    for (ix = 1; ix < 6; ix++) {
        for (jx = 0; jx < FLINT_MAX(flint_test_multiplier(), 1); jx++)
        {
            fq_ctx_t ctx1, ctx2;
            fq_t gen1, gen2;
            fmpz_mod_poly_t minpoly;
            fmpz_mod_poly_t modulus2;
            fmpz_mod_mat_t embed, project, comp, one;
            slong m, n;

            while (fq_ctx_init_randtest(ctx1, state, 3),
                    m = fq_ctx_degree(ctx1),
                    m == 1)
            {
                fq_ctx_clear(ctx1);
            }

            n = m * ix;

            fmpz_mod_poly_init(modulus2, ctx1->ctxp);
            fmpz_mod_poly_randtest_monic_irreducible(modulus2, state, n+1, ctx1->ctxp);
            fq_ctx_init_modulus(ctx2, modulus2, ctx1->ctxp, "X");

            fq_init(gen1, ctx1);
            fq_init(gen2, ctx2);
            fmpz_mod_poly_init(minpoly, ctx1->ctxp);
            fq_embed_gens(gen1, gen2, minpoly, ctx1, ctx2);

            fmpz_mod_mat_init(embed, n, m, ctx1->ctxp);
            fmpz_mod_mat_init(project, m, n, ctx1->ctxp);
            fmpz_mod_mat_init(comp, m, m, ctx1->ctxp);
            fmpz_mod_mat_init(one, m, m, ctx1->ctxp);

            fq_embed_matrices(embed, project, gen1, ctx1, gen2, ctx2, minpoly);

            fmpz_mod_mat_mul(comp, project, embed, ctx1->ctxp);
            fmpz_mod_mat_one(one, ctx1->ctxp);
            if (!fmpz_mod_mat_equal(comp, one, ctx1->ctxp)) {
                flint_printf("FAIL:\n\n");
                flint_printf("CTX 1\n"), fq_ctx_print(ctx1), flint_printf("\n");
                flint_printf("CTX 2\n"), fq_ctx_print(ctx2), flint_printf("\n");
                flint_printf("Embed\n"),
                    fmpz_mod_mat_print_pretty(embed, ctx1->ctxp), flint_printf("\nProject\n"),
                    fmpz_mod_mat_print_pretty(project, ctx1->ctxp), flint_printf("\nComposition\n"),
                    fmpz_mod_mat_print_pretty(comp, ctx1->ctxp), flint_printf("\n");
                flint_abort();
            }

            fmpz_mod_mat_clear(embed, ctx1->ctxp);
            fmpz_mod_mat_clear(project, ctx1->ctxp);
            fmpz_mod_mat_clear(comp, ctx1->ctxp);
            fmpz_mod_mat_clear(one, ctx1->ctxp);
            fmpz_mod_poly_clear(minpoly, ctx1->ctxp);
            fmpz_mod_poly_clear(modulus2, ctx1->ctxp);
            fq_clear(gen1, ctx1);
            fq_ctx_clear(ctx1);
            fq_clear(gen2, ctx2);
            fq_ctx_clear(ctx2);
        }
    }

    TEST_FUNCTION_END(state);
}
