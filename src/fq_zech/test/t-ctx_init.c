/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fq_nmod.h"
#include "fq_zech.h"

TEST_FUNCTION_START(fq_zech_ctx_init, state)
{
    slong primes[10] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29 };
    slong exponents[10] = { 16, 10, 6, 5, 4, 4, 3, 3, 3, 3 };
    int i, j, random;
    slong d;
    fmpz_t p, e;
    fq_nmod_ctx_struct *fq_nmod_ctx;
    fq_nmod_t lhs, rhs, one;
    fq_zech_ctx_t ctx;

    fmpz_init(p);
    fmpz_init(e);

    for (i = 0; i < 10; i++)
    {
        fmpz_set_ui(p, primes[i]);

        for (d = 2; d < exponents[i]; d++)
        {
            for (random = 0; random <= 1; random++)
            {
                if (random)
                    fq_zech_ctx_init_random(ctx, p, d, "a");
                else
                    fq_zech_ctx_init_conway(ctx, p, d, "a");
                fq_nmod_ctx = ctx->fq_nmod_ctx;
                fq_nmod_init(lhs, fq_nmod_ctx);
                fq_nmod_init(rhs, fq_nmod_ctx);
                fq_nmod_init(one, fq_nmod_ctx);

                fq_nmod_one(one, fq_nmod_ctx);

                for (j = 0; j < ctx->qm1; j++)
                {
                    /* Skip the cases where a^j + 1 == 0 */
                    if (primes[i] == 2 && i == 0)
                    {
                        continue;
                    }
                    if (j == ctx->qm1 / 2)
                    {
                        continue;
                    }

                    /* lhs = a^Z(j) */
                    fmpz_set_ui(e, ctx->zech_log_table[j]);
                    fq_nmod_gen(lhs, fq_nmod_ctx);
                    fq_nmod_pow(lhs, lhs, e, fq_nmod_ctx);

                    /* rhs = a^j + 1 */
                    fmpz_set_ui(e, j);
                    fq_nmod_gen(rhs, fq_nmod_ctx);
                    fq_nmod_pow(rhs, rhs, e, fq_nmod_ctx);
                    fq_nmod_add(rhs, rhs, one, fq_nmod_ctx);

                    if (!fq_nmod_equal(lhs, rhs, fq_nmod_ctx))
                    {
                        flint_printf("FAIL:\n\n");
                        flint_printf("K = GF(%wd^%wd)\n", primes[i], d);
                        flint_printf("Z(%d) = %wd\n", j, ctx->zech_log_table[j]);
                        flint_printf("LHS: ");
                        fq_nmod_print_pretty(lhs, fq_nmod_ctx);
                        flint_printf("\n");
                        flint_printf("RHS: ");
                        fq_nmod_print_pretty(rhs, fq_nmod_ctx);
                        flint_printf("\n");
                        fflush(stdout);
                        flint_abort();
                    }
                }

                fq_nmod_clear(lhs, fq_nmod_ctx);
                fq_nmod_clear(rhs, fq_nmod_ctx);
                fq_nmod_clear(one, fq_nmod_ctx);

                fq_zech_ctx_clear(ctx);
            }
        }
    }

    fmpz_clear(p);
    fmpz_clear(e);

    TEST_FUNCTION_END(state);
}
