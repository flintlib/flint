/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_mpoly.h"

/* Defined in t-pow_ui.c and t-pow_rmul.c */
#ifndef nmod_mpoly_pow_naive
#define nmod_mpoly_pow_naive nmod_mpoly_pow_naive
void nmod_mpoly_pow_naive(nmod_mpoly_t res, nmod_mpoly_t f,
                                                 slong n, nmod_mpoly_ctx_t ctx)
{
   if (n == 0)
      nmod_mpoly_set_ui(res, 1, ctx);
   else if (f->length == 0)
      nmod_mpoly_zero(res, ctx);
   else if (n == 1)
      nmod_mpoly_set(res, f, ctx);
   else
   {
      slong i;
      nmod_mpoly_t pow;

      nmod_mpoly_init(pow, ctx);
      nmod_mpoly_set(pow, f, ctx);

      for (i = 1; i < n - 1; i++)
         nmod_mpoly_mul_johnson(pow, pow, f, ctx);

      nmod_mpoly_mul_johnson(res, pow, f, ctx);

      nmod_mpoly_clear(pow, ctx);
   }
}
#endif

TEST_FUNCTION_START(nmod_mpoly_pow_rmul, state)
{
    slong i, j;

    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b;
        nmod_mpoly_ctx_init(ctx, 20, ORD_LEX, 4);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_set_str_pretty(a, "1", NULL, ctx);
        nmod_mpoly_set_str_pretty(b, "2 + 2*x1 + 2*x2 + 2*x3", NULL, ctx);
        nmod_mpoly_pow_rmul(a, b, 3, ctx);
        if (nmod_mpoly_is_zero(b, ctx) || !nmod_mpoly_is_zero(a, ctx))
        {
            flint_printf("FAIL: Check simple example\n");
            fflush(stdout);
            flint_abort();
        }
        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h;
        ulong pow_bound;
        slong len, len1, len2;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;
        mp_limb_t modulus;

        modulus = n_randbits(state, n_randint(state, FLINT_BITS));
        modulus = FLINT_MAX(UWORD(2), modulus);

        nmod_mpoly_ctx_init_rand(ctx, state, 10, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 10);
        len2 = n_randint(state, 10);

        exp_bits = n_randint(state, 20) + 2;
        exp_bits1 = n_randint(state, 20) + 2;
        exp_bits2 = n_randint(state, 20) + 2;

        pow_bound = 250/(len1+1);
        pow_bound = pow_bound/FLINT_MAX(WORD(1), ctx->minfo->nvars);
        pow_bound = FLINT_MAX(pow_bound, UWORD(5));

        for (j = 0; j < 10; j++)
        {
            slong pow;

            pow = n_randint(state, pow_bound);

            nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            nmod_mpoly_randtest_bits(h, state, len, exp_bits, ctx);

            nmod_mpoly_pow_rmul(g, f, pow, ctx);
            nmod_mpoly_assert_canonical(g, ctx);

            nmod_mpoly_pow_naive(h, f, pow, ctx);
            nmod_mpoly_assert_canonical(h, ctx);

            if (!nmod_mpoly_equal(g, h, ctx))
            {
                flint_printf("FAIL: Check pow_ui against pow_naive\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            nmod_mpoly_pow_rmul(f, f, pow, ctx);
            nmod_mpoly_assert_canonical(f, ctx);

            if (!nmod_mpoly_equal(g, f, ctx))
            {
                flint_printf("FAIL: Check aliasing\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
