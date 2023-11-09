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
    int i, j;

    /* Check that isomorphism to self gives identity matrices */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fq_ctx_t ctx;
        fq_t gen;
        const fmpz_mod_poly_struct *modulus;
        fmpz_mod_mat_t embed, project, one;
        slong d;

        fq_ctx_randtest(ctx, state);
        d = fq_ctx_degree(ctx);
        modulus = fq_ctx_modulus(ctx);

        fq_init(gen, ctx);
        fq_gen(gen, ctx);
        fq_pow(gen, gen, fq_ctx_prime(ctx), ctx);

        fmpz_mod_mat_init(embed, d, d, fq_ctx_prime(ctx));
        fmpz_mod_mat_init(project, d, d, fq_ctx_prime(ctx));
        fmpz_mod_mat_init(one, d, d, fq_ctx_prime(ctx));

        fq_embed_matrices(embed, project, gen, ctx, gen, ctx, modulus);
        fmpz_mod_mat_one(one);

        if (!fmpz_mod_mat_equal(embed, one) || !fmpz_mod_mat_equal(project, one)) {
            flint_printf("FAIL:\n\n");
            flint_printf("CTX\n"), fq_ctx_print(ctx), flint_printf("\n");
            flint_printf("x^p: "), fq_print_pretty(gen, ctx), flint_printf("\n");
            flint_printf("Embed\n"),
                fmpz_mod_mat_print_pretty(embed), flint_printf("\nProject\n"),
                fmpz_mod_mat_print_pretty(project), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_mat_clear(embed);
        fmpz_mod_mat_clear(project);
        fmpz_mod_mat_clear(one);
        fq_clear(gen, ctx);
        fq_ctx_clear(ctx);
    }

    /* Check random emebedding (degrees 1..5) */
    for (j = 1; j < 6; j++) {
        for (i = 0; i < (6 - j) * flint_test_multiplier(); i++)
        {
            fq_ctx_t ctx1, ctx2;
            fq_t gen1, gen2;
            fmpz_mod_poly_t minpoly;
            fmpz_mod_poly_t modulus2;
            fmpz_mod_mat_t embed, project, comp, one;
            slong m, n;

            while (fq_ctx_randtest(ctx1, state),
                    m = fq_ctx_degree(ctx1),
                    m == 1)
            {
                fq_ctx_clear(ctx1);
            }

            n = m*j;

            fmpz_mod_poly_init(modulus2, ctx1->ctxp);
            fmpz_mod_poly_randtest_monic_irreducible(modulus2, state, n+1, ctx1->ctxp);
            fq_ctx_init_modulus(ctx2, modulus2, ctx1->ctxp, "X");

            fq_init(gen1, ctx1);
            fq_init(gen2, ctx2);
            fmpz_mod_poly_init(minpoly, ctx1->ctxp);
            fq_embed_gens(gen1, gen2, minpoly, ctx1, ctx2);

            fmpz_mod_mat_init(embed, n, m, fq_ctx_prime(ctx1));
            fmpz_mod_mat_init(project, m, n, fq_ctx_prime(ctx1));
            fmpz_mod_mat_init(comp, m, m, fq_ctx_prime(ctx1));
            fmpz_mod_mat_init(one, m, m, fq_ctx_prime(ctx1));

            fq_embed_matrices(embed, project, gen1, ctx1, gen2, ctx2, minpoly);

            fmpz_mod_mat_mul(comp, project, embed);
            fmpz_mod_mat_one(one);
            if (!fmpz_mod_mat_equal(comp, one)) {
                flint_printf("FAIL:\n\n");
                flint_printf("CTX 1\n"), fq_ctx_print(ctx1), flint_printf("\n");
                flint_printf("CTX 2\n"), fq_ctx_print(ctx2), flint_printf("\n");
                flint_printf("Embed\n"),
                    fmpz_mod_mat_print_pretty(embed), flint_printf("\nProject\n"),
                    fmpz_mod_mat_print_pretty(project), flint_printf("\nComposition\n"),
                    fmpz_mod_mat_print_pretty(comp), flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            fmpz_mod_mat_clear(embed);
            fmpz_mod_mat_clear(project);
            fmpz_mod_mat_clear(comp);
            fmpz_mod_mat_clear(one);
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
