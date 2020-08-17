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
#include "fmpz_mpoly.h"
#include "ulong_extras.h"


int
main(void)
{
    int result;
    slong i, j1, j2;
    FLINT_TEST_INIT(state);

    flint_printf("cmp....");
    fflush(stdout);

    /* check polynomial terms are in order */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, mf, mg;
        slong len;
        flint_bitcnt_t coeff_bits, exp_bits;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(mf, ctx);
        fmpz_mpoly_init(mg, ctx);

        len = n_randint(state, 100) + 1;
        exp_bits = n_randint(state, 200) + 2;
        exp_bits = n_randint(state, exp_bits) + 2;
        coeff_bits = n_randint(state, 100) + 1;

        fmpz_mpoly_randtest_bits(f, state, len, coeff_bits, exp_bits, ctx);
        fmpz_mpoly_repack_bits(g, f, f->bits + n_randint(state, FLINT_BITS), ctx);
        fmpz_mpoly_assert_canonical(f, ctx);
        fmpz_mpoly_assert_canonical(g, ctx);

        for (j1 = 0; j1 < f->length; j1++)
        for (j2 = 0; j2 < g->length; j2++)
        {
            fmpz_mpoly_get_term_monomial(mf, f, j1, ctx);
            fmpz_mpoly_get_term_monomial(mg, g, j2, ctx);
            result = fmpz_mpoly_cmp(mf, mg, ctx);
            result =   (result == 0 && j1 == j2)
                    || (result == +1 && j1 < j2)
                    || (result == -1 && j1 > j2);
            if (!result)
            {
                flint_printf("FAIL\n"
                             "check polynomial terms are in order\n"
                             "i = %wd, j1 = %wd, j2 = %wd\n", i, j1, j2);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(mf, ctx);
        fmpz_mpoly_clear(mg, ctx);
        fmpz_mpoly_ctx_clear(ctx);  
    }

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, c, aa, bb, cc;
        int a_a, a_b, a_c, a_aa, a_bb, a_cc;
        int b_a, b_b, b_c, b_aa, b_bb, b_cc;
        int c_a, c_b, c_c, c_aa, c_bb, c_cc;
        flint_bitcnt_t newbits;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);

        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(c, ctx);
        fmpz_mpoly_init(aa, ctx);
        fmpz_mpoly_init(bb, ctx);
        fmpz_mpoly_init(cc, ctx);

        fmpz_mpoly_randtest_bits(a, state, n_randint(state, 100) + 1,
                                           n_randint(state, 100) + 1,
                                           n_randint(state, 200) + 2, ctx);

        fmpz_mpoly_randtest_bits(b, state, n_randint(state, 100) + 1,
                                           n_randint(state, 100) + 1,
                                           n_randint(state, 200) + 2, ctx);

        fmpz_mpoly_randtest_bits(c, state, n_randint(state, 100) + 1,
                                           n_randint(state, 100) + 1,
                                           n_randint(state, 200) + 2, ctx);

        newbits = a->bits + n_randint(state, 2*FLINT_BITS);
        newbits = mpoly_fix_bits(newbits, ctx->minfo);
        fmpz_mpoly_repack_bits(aa, a, newbits, ctx);

        newbits = b->bits + n_randint(state, 2*FLINT_BITS);
        newbits = mpoly_fix_bits(newbits, ctx->minfo);
        fmpz_mpoly_repack_bits(bb, b, newbits, ctx);

        newbits = c->bits + n_randint(state, 2*FLINT_BITS);
        newbits = mpoly_fix_bits(newbits, ctx->minfo);
        fmpz_mpoly_repack_bits(cc, c, newbits, ctx);

        a_a = fmpz_mpoly_cmp(a, a, ctx);
        a_b = fmpz_mpoly_cmp(a, b, ctx);
        a_c = fmpz_mpoly_cmp(a, c, ctx);
        a_aa = fmpz_mpoly_cmp(a, aa, ctx);
        a_bb = fmpz_mpoly_cmp(a, bb, ctx);
        a_cc = fmpz_mpoly_cmp(a, cc, ctx);

        b_a = fmpz_mpoly_cmp(b, a, ctx);
        b_b = fmpz_mpoly_cmp(b, b, ctx);
        b_c = fmpz_mpoly_cmp(b, c, ctx);
        b_aa = fmpz_mpoly_cmp(b, aa, ctx);
        b_bb = fmpz_mpoly_cmp(b, bb, ctx);
        b_cc = fmpz_mpoly_cmp(b, cc, ctx);

        c_a = fmpz_mpoly_cmp(c, a, ctx);
        c_b = fmpz_mpoly_cmp(c, b, ctx);
        c_c = fmpz_mpoly_cmp(c, c, ctx);
        c_aa = fmpz_mpoly_cmp(c, aa, ctx);
        c_bb = fmpz_mpoly_cmp(c, bb, ctx);
        c_cc = fmpz_mpoly_cmp(c, cc, ctx);

        if (a_a != 0 || a_aa != 0 ||
            b_b != 0 || b_bb != 0 ||
            c_c != 0 || c_cc != 0)
        {
            flint_printf("FAIL\n"
                         "check polynomial compares equal to itself\n"
                         "i = %wd\n", i);
            flint_abort();
        }

        if (a_b != a_bb || a_c != a_cc ||
            b_a != b_aa || b_c != b_cc ||
            c_a != c_aa || c_b != c_bb)
        {
            flint_printf("FAIL\n"
                         "check polynomial comparison with differing bits\n"
                         "i = %wd\n", i);
            flint_abort();
        }

        if ((a_b*b_c == 0 && a_c != a_b + b_c) || (a_b*b_c > 0 && a_c != a_b))
        {
            flint_printf("FAIL\n"
                         "check transitivity\n"
                         "i = %wd\n", i);
            flint_abort();
        }

        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(c, ctx);
        fmpz_mpoly_clear(aa, ctx);
        fmpz_mpoly_clear(bb, ctx);
        fmpz_mpoly_clear(cc, ctx);

        fmpz_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
