/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mod_mat.h"
#include "fq.h"
#include "fq_embed.h"

TEST_FUNCTION_START(fq_embed_composition_matrix, state)
{
    int i;

    /* Check that Mat(a^p) = Mat(x^p) * Mat(a) for random a */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fq_ctx_t ctx;
        fq_t frob, a;
        fmpz_mod_mat_t mat_frob, mat_a, mat_aq, res;
        slong d;

        while (fq_ctx_randtest(ctx, state),
               d = fq_ctx_degree(ctx),
               d == 1)
        {
            fq_ctx_clear(ctx);
        }

        fq_init(frob, ctx);
        fq_init(a, ctx);
        fmpz_mod_mat_init(mat_frob, d, d, fq_ctx_prime(ctx));
        fmpz_mod_mat_init(mat_a, d, d, fq_ctx_prime(ctx));
        fmpz_mod_mat_init(mat_aq, d, d, fq_ctx_prime(ctx));
        fmpz_mod_mat_init(res, d, d, fq_ctx_prime(ctx));

        fq_gen(frob, ctx);
        fq_pow(frob, frob, fq_ctx_prime(ctx), ctx);
        fq_embed_composition_matrix(mat_frob, frob, ctx);

        fq_randtest(a, state, ctx);
        fq_embed_composition_matrix(mat_a, a, ctx);

        fmpz_mod_mat_mul(res, mat_frob, mat_a);

        fq_pow(a, a, fq_ctx_prime(ctx), ctx);
        fq_embed_composition_matrix(mat_aq, a, ctx);

        if (!fmpz_mod_mat_equal(res, mat_aq))
        {
            flint_printf("FAIL:\n\n");
            flint_printf("CTX\n"), fq_ctx_print(ctx), flint_printf("\n");
            flint_printf("x^q: "), fq_print_pretty(frob, ctx), flint_printf("\n");
            flint_printf("M(x^q)*M(a) = M(a^q)\n"),
                fmpz_mod_mat_print_pretty(mat_frob), flint_printf("\n"),
                fmpz_mod_mat_print_pretty(mat_a), flint_printf("\n"),
                fmpz_mod_mat_print_pretty(mat_aq), flint_printf("\n"),
                fmpz_mod_mat_print_pretty(res), flint_printf("\n");

            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_mat_clear(mat_frob);
        fmpz_mod_mat_clear(mat_a);
        fmpz_mod_mat_clear(mat_aq);
        fmpz_mod_mat_clear(res);
        fq_clear(frob, ctx);
        fq_clear(a, ctx);
        fq_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
