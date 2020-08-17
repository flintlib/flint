/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "fmpq_mpoly.h"
#include "ulong_extras.h"

int
main(void)
{
    slong i, j1, j2;
    FLINT_TEST_INIT(state);

    flint_printf("get_coeff_vars_ui....");
    fflush(stdout);

    /* check simple example */
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g;
        /* get the coefficient of y^1*x^2*/
        slong varl[2] = {1, 0};
        ulong expl[2] = {1, 2};
        const char * vars[] = {"x", "y", "z", "w"};

        fmpq_mpoly_ctx_init(ctx, 4, ORD_DEGREVLEX);
        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_set_str_pretty(f, "x^2*y*(z+w)+x+y+x*y^2+z^2+w^2", vars, ctx);
        fmpq_mpoly_set_str_pretty(g, "z+w", vars, ctx);
        fmpq_mpoly_get_coeff_vars_ui(f, f, varl, expl, 2, ctx);
        if (!fmpq_mpoly_equal(f, g, ctx))
        {
            flint_printf("FAIL\ncheck simple example\n");
            flint_abort();
        }
        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    /* check 1 variable sum of coefficients */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, h, m;
        slong nvars, len;
        ulong exp_bound;
        flint_bitcnt_t coeff_bits;
        slong vars[1];
        ulong exps[1];
        slong var1;

        nvars = 1 + n_randint(state, 20);
        fmpq_mpoly_ctx_init(ctx, nvars, mpoly_ordering_randtest(state));

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(h, ctx);
        fmpq_mpoly_init(m, ctx);

        len = n_randint(state, 200);
        exp_bound = n_randint(state, 20) + 1;
        coeff_bits = n_randint(state, 100);

        fmpq_mpoly_randtest_bound(f, state, len, coeff_bits, exp_bound, ctx);
        fmpq_mpoly_repack_bits(f, f,
                           f->zpoly->bits + n_randint(state, FLINT_BITS), ctx);

        var1 = n_randint(state, nvars);

        fmpq_mpoly_zero(h, ctx);
        for (j1 = 0; j1 < exp_bound; j1++)
        {
            vars[0] = var1;
            exps[0] = j1;
            fmpq_mpoly_get_coeff_vars_ui(g, f, vars, exps, 1, ctx);
            fmpq_mpoly_assert_canonical(g, ctx);
            fmpq_mpoly_gen(m, var1, ctx);
            fmpq_mpoly_pow_ui(m, m, j1, ctx);
            fmpq_mpoly_mul(g, g, m, ctx);
            fmpq_mpoly_add(h, h, g, ctx);
        }

        if (!fmpq_mpoly_equal(f, h, ctx))
        {
            flint_printf("FAIL\n"
                         "check 1 variable sum of coefficients\ni = %wd\n", i);
            flint_abort();
        }

        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(h, ctx);
        fmpq_mpoly_clear(m, ctx);
        fmpq_mpoly_ctx_clear(ctx);  
    }

    /* check 2 variable sum of coefficients */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, h, m;
        slong nvars, len;
        ulong exp_bound;
        flint_bitcnt_t coeff_bits;
        slong vars[2];
        ulong exps[2];
        slong var1, var2;

        nvars = 2 + n_randint(state, 20);
        fmpq_mpoly_ctx_init(ctx, nvars, mpoly_ordering_randtest(state));

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(h, ctx);
        fmpq_mpoly_init(m, ctx);

        len = n_randint(state, 200);
        exp_bound = n_randint(state, 12) + 1;
        coeff_bits = n_randint(state, 100);

        fmpq_mpoly_randtest_bound(f, state, len, coeff_bits, exp_bound, ctx);
        fmpq_mpoly_repack_bits(f, f,
                           f->zpoly->bits + n_randint(state, FLINT_BITS), ctx);

        var1 = n_randint(state, nvars - 1);
        var2 = 1 + var1 + n_randint(state, nvars - 1 - var1);

        fmpq_mpoly_zero(h, ctx);
        for (j1 = 0; j1 < exp_bound; j1++)
        {
        for (j2 = 0; j2 < exp_bound; j2++)
        {
            vars[0] = var1;
            exps[0] = j1;
            vars[1] = var2;
            exps[1] = j2;
            fmpq_mpoly_get_coeff_vars_ui(g, f, vars, exps, 2, ctx);
            fmpq_mpoly_assert_canonical(g, ctx);
            fmpq_mpoly_gen(m, var1, ctx);
            fmpq_mpoly_pow_ui(m, m, j1, ctx);
            fmpq_mpoly_mul(g, g, m, ctx);
            fmpq_mpoly_gen(m, var2, ctx);
            fmpq_mpoly_pow_ui(m, m, j2, ctx);
            fmpq_mpoly_mul(g, g, m, ctx);
            fmpq_mpoly_add(h, h, g, ctx);
        }
        }

        if (!fmpq_mpoly_equal(f, h, ctx))
        {
            flint_printf("FAIL\n"
                         "check 2 variable sum of coefficients\ni = %wd\n", i);
            flint_abort();
        }

        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(h, ctx);
        fmpq_mpoly_clear(m, ctx);
        fmpq_mpoly_ctx_clear(ctx);  
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
