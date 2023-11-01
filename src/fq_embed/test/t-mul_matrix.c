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

TEST_FUNCTION_START(fq_embed_mul_matrix, state)
{
    int i;

    /* Check that Mat(a^2) = Mat(a)^2 for random a */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fq_ctx_t ctx;
        fq_t a;
        fmpz_mod_mat_t mat_a, mat_a_sq, mat_a_a;
        slong d;

        fq_ctx_randtest(ctx, state);
        d = fq_ctx_degree(ctx);

        fq_init(a, ctx);
        fmpz_mod_mat_init(mat_a, d, d, fq_ctx_prime(ctx));
        fmpz_mod_mat_init(mat_a_sq, d, d, fq_ctx_prime(ctx));
        fmpz_mod_mat_init(mat_a_a, d, d, fq_ctx_prime(ctx));

        fq_randtest(a, state, ctx);
        fq_embed_mul_matrix(mat_a, a, ctx);

        fq_mul(a, a, a, ctx);
        fq_embed_mul_matrix(mat_a_sq, a, ctx);

        fmpz_mod_mat_mul(mat_a_a, mat_a, mat_a);

        if (!fmpz_mod_mat_equal(mat_a_a, mat_a_sq))
        {
            flint_printf("FAIL:\n\n");
            flint_printf("CTX\n"), fq_ctx_print(ctx), flint_printf("\n");
            flint_printf("a^2: "), fq_print_pretty(a, ctx), flint_printf("\n");
            flint_printf("M(a)^2 = M(a^2)\n"),
                fmpz_mod_mat_print_pretty(mat_a), flint_printf("^2\n=\n"),
                fmpz_mod_mat_print_pretty(mat_a_a), flint_printf("\n=\n"),
                fmpz_mod_mat_print_pretty(mat_a_sq), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_mat_clear(mat_a);
        fmpz_mod_mat_clear(mat_a_sq);
        fmpz_mod_mat_clear(mat_a_a);
        fq_clear(a, ctx);
        fq_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
