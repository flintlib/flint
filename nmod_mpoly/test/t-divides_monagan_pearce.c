/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "nmod_mpoly.h"

int
main(void)
{
    int i, j, result, ok1, ok2;
    FLINT_TEST_INIT(state);

    flint_printf("divides_monagan_pearce....");
    fflush(stdout);

    /* Check f*g/g = f */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h, k;
        ordering_t ord;
        mp_limb_t modulus;
        slong maxbits;
        slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;
        slong exp_bits, exp_bits1, exp_bits2;

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 10) + 1;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);
        nmod_mpoly_ctx_init(ctx, nvars, ord, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);
        nmod_mpoly_init(k, ctx);

        len = n_randint(state, 4);
        len1 = n_randint(state, 4);
        len2 = n_randint(state, 4) + 1;

        maxbits = FLINT_BITS - mpoly_ordering_isdeg(ord)*FLINT_BIT_COUNT(nvars);
        exp_bits = n_randint(state, maxbits - 1) + 1;
        exp_bits1 = n_randint(state, maxbits - 2) + 1;
        exp_bits2 = n_randint(state, maxbits - 2) + 1;

        exp_bound = n_randbits(state, exp_bits);
        exp_bound1 = n_randbits(state, exp_bits1);
        exp_bound2 = n_randbits(state, exp_bits2);

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest(f, state, len1, exp_bound1, ctx);
            do {
                nmod_mpoly_randtest(g, state, len2, exp_bound2, ctx);
            } while (g->length == 0);
            nmod_mpoly_randtest(h, state, len, exp_bound, ctx);
            nmod_mpoly_randtest(k, state, len, exp_bound, ctx);

            nmod_mpoly_mul_johnson(h, f, g, ctx);
            nmod_mpoly_test(h, ctx);

            ok1 = nmod_mpoly_divides_monagan_pearce(k, h, g, ctx);
            nmod_mpoly_test(k, ctx);

            result = (ok1 && nmod_mpoly_equal(f, k, ctx));

            if (!result)
            {
                flint_printf("FAIL\n");
                flint_printf("Check f*g/g = f\ni = %wd, j = %wd\n", i, j);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_clear(k, ctx);

        nmod_mpoly_ctx_clear(ctx);

    }

    /* Check random polys don't divide */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h, k;
        ordering_t ord;
        mp_limb_t modulus;
        slong maxbits;
        slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;
        slong exp_bits, exp_bits1, exp_bits2;

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 10) + 1;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);
        nmod_mpoly_ctx_init(ctx, nvars, ord, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);
        nmod_mpoly_init(k, ctx);

        len = n_randint(state, 20);
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20) + 1;

        maxbits = 20/(nvars + mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1;
        exp_bits = n_randint(state, maxbits) + 1;
        exp_bits1 = n_randint(state, maxbits) + 1;
        exp_bits2 = n_randint(state, maxbits) + 1;

        exp_bound = n_randbits(state, exp_bits);
        exp_bound1 = n_randbits(state, exp_bits1);
        exp_bound2 = n_randbits(state, exp_bits2);

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest(f, state, len1, exp_bound1, ctx);
            do {
                nmod_mpoly_randtest(g, state, len2, exp_bound2, ctx);
            } while (g->length == 0);
            nmod_mpoly_randtest(h, state, len, exp_bound, ctx);
            nmod_mpoly_randtest(k, state, len, exp_bound, ctx);

            ok1 = nmod_mpoly_divides_monagan_pearce(h, f, g, ctx);
            nmod_mpoly_test(h, ctx);

            if (ok1)
            {
                nmod_mpoly_mul_johnson(k, h, g, ctx);
                nmod_mpoly_test(k, ctx);
            }

            result = (ok1 == 0 || nmod_mpoly_equal(f, k, ctx));

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check random polys don't divide\ni = %wd, j = %wd\n", i, j);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_clear(k, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing first argument, exact division */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h, k;
        ordering_t ord;
        mp_limb_t modulus;
        slong maxbits;
        slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;
        slong exp_bits, exp_bits1, exp_bits2;

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 10) + 1;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);
        nmod_mpoly_ctx_init(ctx, nvars, ord, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);
        nmod_mpoly_init(k, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100) + 1;

        maxbits = FLINT_BITS - mpoly_ordering_isdeg(ord)*FLINT_BIT_COUNT(nvars);
        exp_bits = n_randint(state, maxbits - 1) + 1;
        exp_bits1 = n_randint(state, maxbits - 2) + 1;
        exp_bits2 = n_randint(state, maxbits - 2) + 1;

        exp_bound = n_randbits(state, exp_bits);
        exp_bound1 = n_randbits(state, exp_bits1);
        exp_bound2 = n_randbits(state, exp_bits2);

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest(f, state, len1, exp_bound1, ctx);
            do {
                nmod_mpoly_randtest(g, state, len2, exp_bound2, ctx);
            } while (g->length == 0);
            nmod_mpoly_randtest(h, state, len, exp_bound, ctx);
            nmod_mpoly_randtest(k, state, len, exp_bound, ctx);

            nmod_mpoly_mul_johnson(h, f, g, ctx);
            nmod_mpoly_test(h, ctx);

            ok1 = nmod_mpoly_divides_monagan_pearce(k, h, g, ctx);
            nmod_mpoly_test(k, ctx);
            ok2 = nmod_mpoly_divides_monagan_pearce(h, h, g, ctx);
            nmod_mpoly_test(h, ctx);

            result = (ok1 == 1 && ok2 == 1 && nmod_mpoly_equal(h, k, ctx));

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing first argument, exact division\ni = %wd, j = %wd\n", i, j);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_clear(k, ctx);
        
        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing, first argument, random polys */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h;
        ordering_t ord;
        mp_limb_t modulus;
        slong maxbits;
        slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;
        slong exp_bits, exp_bits1, exp_bits2;

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 10) + 1;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);
        nmod_mpoly_ctx_init(ctx, nvars, ord, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);

        len = n_randint(state, 20);
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20) + 1;

        maxbits = 20/(nvars + mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1;
        exp_bits = n_randint(state, maxbits) + 1;
        exp_bits1 = n_randint(state, maxbits) + 1;
        exp_bits2 = n_randint(state, maxbits) + 1;

        exp_bound = n_randbits(state, exp_bits);
        exp_bound1 = n_randbits(state, exp_bits1);
        exp_bound2 = n_randbits(state, exp_bits2);

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest(f, state, len1, exp_bound1, ctx);
            do {
                nmod_mpoly_randtest(g, state, len2, exp_bound2, ctx);
            } while (g->length == 0);
            nmod_mpoly_randtest(h, state, len, exp_bound, ctx);

            ok1 = nmod_mpoly_divides_monagan_pearce(h, f, g, ctx);
            nmod_mpoly_test(h, ctx);
            ok2 = nmod_mpoly_divides_monagan_pearce(f, f, g, ctx);
            nmod_mpoly_test(f, ctx);

            result = ((ok1 == ok2) &&  (ok1 == 0 || nmod_mpoly_equal(f, h, ctx)));

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing, first argument, random polys\ni = %wd, j = %wd\n", i, j);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);

        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing second argument, exact division */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h, k;
        ordering_t ord;
        mp_limb_t modulus;
        slong maxbits;
        slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;
        slong exp_bits, exp_bits1, exp_bits2;

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 10) + 1;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);
        nmod_mpoly_ctx_init(ctx, nvars, ord, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);
        nmod_mpoly_init(k, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100) + 1;

        maxbits = FLINT_BITS - mpoly_ordering_isdeg(ord)*FLINT_BIT_COUNT(nvars);
        exp_bits = n_randint(state, maxbits - 1) + 1;
        exp_bits1 = n_randint(state, maxbits - 2) + 1;
        exp_bits2 = n_randint(state, maxbits - 2) + 1;

        exp_bound = n_randbits(state, exp_bits);
        exp_bound1 = n_randbits(state, exp_bits1);
        exp_bound2 = n_randbits(state, exp_bits2);

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest(f, state, len1, exp_bound1, ctx);
            do {
                nmod_mpoly_randtest(g, state, len2, exp_bound2, ctx);
            } while (g->length == 0);
            nmod_mpoly_randtest(h, state, len, exp_bound, ctx);
            nmod_mpoly_randtest(k, state, len, exp_bound, ctx);

            nmod_mpoly_mul_johnson(h, f, g, ctx);
            nmod_mpoly_test(h, ctx);

            ok1 = nmod_mpoly_divides_monagan_pearce(k, h, g, ctx);
            nmod_mpoly_test(k, ctx);
            ok2 = nmod_mpoly_divides_monagan_pearce(g, h, g, ctx);
            nmod_mpoly_test(g, ctx);

            result = (ok1 == 1 && ok2 == 1 && nmod_mpoly_equal(g, k, ctx));

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing second argument, exact division\ni = %wd, j = %wd\n", i, j);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_clear(k, ctx);

        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing, second argument, random polys */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h;
        ordering_t ord;
        mp_limb_t modulus;
        slong maxbits;
        slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;
        slong exp_bits, exp_bits1, exp_bits2;

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 10) + 1;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);
        nmod_mpoly_ctx_init(ctx, nvars, ord, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100) + 1;

        maxbits = 20/(nvars + mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1;
        exp_bits = n_randint(state, maxbits) + 1;
        exp_bits1 = n_randint(state, maxbits) + 1;
        exp_bits2 = n_randint(state, maxbits) + 1;

        exp_bound = n_randbits(state, exp_bits);
        exp_bound1 = n_randbits(state, exp_bits1);
        exp_bound2 = n_randbits(state, exp_bits2);

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest(f, state, len1, exp_bound1, ctx);
            do {
                nmod_mpoly_randtest(g, state, len2, exp_bound2, ctx);
            } while (g->length == 0);
            nmod_mpoly_randtest(h, state, len, exp_bound, ctx);

            ok1 = nmod_mpoly_divides_monagan_pearce(h, f, g, ctx);
            nmod_mpoly_test(h, ctx);
            ok2 = nmod_mpoly_divides_monagan_pearce(g, f, g, ctx);
            nmod_mpoly_test(g, ctx);

            result = ((ok1 == ok2) &&  (ok1 == 0 || nmod_mpoly_equal(g, h, ctx)));

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing, second argument, random polys\ni = %wd, j = %wd\n", i, j);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);

        nmod_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

