/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "nmod_poly.h"
#include "nmod_mat.h"
#include "fq_zech.h"
#include "fq_zech_embed.h"

TEST_FUNCTION_START(fq_zech_embed_matrices, state)
{
    int i, j;
    int primes[4] = {2, 3, 5};
    int degrees[2] = {2, 3};

    /* Check that isomorphism to self gives identity matrices */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fq_zech_ctx_t ctx;
        fq_zech_t gen;
        const nmod_poly_struct *modulus;
        nmod_mat_t embed, project, one;
        slong d;

        fq_zech_ctx_randtest(ctx, state);
        d = fq_zech_ctx_degree(ctx);
        modulus = fq_zech_ctx_modulus(ctx);

        fq_zech_init(gen, ctx);
        fq_zech_gen(gen, ctx);
        fq_zech_pow(gen, gen, fq_zech_ctx_prime(ctx), ctx);

        nmod_mat_init(embed, d, d, nmod_poly_modulus(modulus));
        nmod_mat_init(project, d, d, nmod_poly_modulus(modulus));
        nmod_mat_init(one, d, d, nmod_poly_modulus(modulus));

        fq_zech_embed_matrices(embed, project, gen, ctx, gen, ctx, modulus);
        nmod_mat_one(one);

        if (!nmod_mat_equal(embed, one) || !nmod_mat_equal(project, one)) {
            flint_printf("FAIL:\n\n");
            flint_printf("CTX\n"), fq_zech_ctx_print(ctx), flint_printf("\n");
            flint_printf("x^p: "), fq_zech_print_pretty(gen, ctx), flint_printf("\n");
            flint_printf("Embed\n"),
                nmod_mat_print_pretty(embed), flint_printf("\nProject\n"),
                nmod_mat_print_pretty(project), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_mat_clear(embed);
        nmod_mat_clear(project);
        nmod_mat_clear(one);
        fq_zech_clear(gen, ctx);
        fq_zech_ctx_clear(ctx);
    }

    /* Check random embedding (degrees 1..3) */
    for (j = 1; j < 4; j++) {
        for (i = 0; i < (4 - j) * flint_test_multiplier(); i++)
        {
            fq_zech_ctx_t ctx1, ctx2;
            fq_zech_t gen1, gen2;
            nmod_poly_t minpoly;
            const nmod_poly_struct *modulus;
            nmod_poly_t modulus2;
            nmod_mat_t embed, project, comp, one;
            slong m, n;
            fmpz_t p;

            fmpz_init(p);
            fmpz_set_si(p, primes[i % 3]);
            m = degrees[i % 2];

            fq_zech_ctx_init_random(ctx1, p, m, "a");

            n = m*j;
            if (m == 1) {
                i--;
                continue;
            }
            modulus = fq_zech_ctx_modulus(ctx1);

            nmod_poly_init(modulus2, nmod_poly_modulus(modulus));
            nmod_poly_randtest_monic_primitive(modulus2, state, n+1);
            fq_zech_ctx_init_modulus(ctx2, modulus2, "X");

            fq_zech_init(gen1, ctx1);
            fq_zech_init(gen2, ctx2);
            nmod_poly_init(minpoly, nmod_poly_modulus(modulus));
            fq_zech_embed_gens(gen1, gen2, minpoly, ctx1, ctx2);

            nmod_mat_init(embed, n, m, nmod_poly_modulus(modulus));
            nmod_mat_init(project, m, n, nmod_poly_modulus(modulus));
            nmod_mat_init(comp, m, m, nmod_poly_modulus(modulus));
            nmod_mat_init(one, m, m, nmod_poly_modulus(modulus));

            fq_zech_embed_matrices(embed, project, gen1, ctx1, gen2, ctx2, minpoly);

            nmod_mat_mul(comp, project, embed);
            nmod_mat_one(one);
            if (!nmod_mat_equal(comp, one)) {
                flint_printf("FAIL:\n\n");
                flint_printf("CTX 1\n"), fq_zech_ctx_print(ctx1), flint_printf("\n");
                flint_printf("CTX 2\n"), fq_zech_ctx_print(ctx2), flint_printf("\n");
                flint_printf("Embed\n"),
                    nmod_mat_print_pretty(embed), flint_printf("\nProject\n"),
                    nmod_mat_print_pretty(project), flint_printf("\nComposition\n"),
                    nmod_mat_print_pretty(comp), flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            nmod_mat_clear(embed);
            nmod_mat_clear(project);
            nmod_mat_clear(comp);
            nmod_mat_clear(one);
            nmod_poly_clear(minpoly);
            nmod_poly_clear(modulus2);
            fq_zech_clear(gen1, ctx1);
            fq_zech_ctx_clear(ctx1);
            fq_zech_clear(gen2, ctx2);
            fq_zech_ctx_clear(ctx2);
            fmpz_clear(p);
        }
    }

    TEST_FUNCTION_END(state);
}
