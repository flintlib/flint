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
#include "nmod_mpoly.h"

int
main(void)
{
    int result;
    slong i, j, w;
    FLINT_TEST_INIT(state);

    flint_printf("divrem_ideal_monagan_pearce....");
    fflush(stdout);

    /* Check f*g/g = f */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h, k, r;
        ordering_t ord;
        mp_limb_t modulus;
        slong maxbits;
        slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;
        slong exp_bits, exp_bits1, exp_bits2;
        nmod_mpoly_struct * qarr[1], * darr[1];
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
        nmod_mpoly_init(r, ctx);

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
                nmod_mpoly_randtest(g, state, len2, exp_bound2 + 1, ctx);
            } while (g->length == 0);
            nmod_mpoly_randtest(h, state, len, exp_bound, ctx);
            nmod_mpoly_randtest(k, state, len, exp_bound, ctx);
            nmod_mpoly_randtest(r, state, len, exp_bound, ctx);

            nmod_mpoly_mul_johnson(h, f, g, ctx);

            qarr[0] = k;
            darr[0] = g;

            nmod_mpoly_divrem_ideal_monagan_pearce(qarr, r, h, darr, 1, ctx);
            nmod_mpoly_test(qarr[0], ctx);
            nmod_mpoly_test(r, ctx);

            result = nmod_mpoly_equal(f, k, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*g/g = f\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_clear(k, ctx);
        nmod_mpoly_clear(r, ctx);

        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check f = g1*q1 + ... + gn*qn + r for random polys */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, r, k1, k2;
        nmod_mpoly_struct * g, * q;
        ordering_t ord;
        mp_limb_t modulus;
        slong maxbits;
        slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2, num;
        slong exp_bits, exp_bits1, exp_bits2;
        nmod_mpoly_struct * qarr[5], * darr[5];

        num = n_randint(state, 5) + 1;

        g = (nmod_mpoly_struct *) flint_malloc(num*sizeof(nmod_mpoly_struct));
        q = (nmod_mpoly_struct *) flint_malloc(num*sizeof(nmod_mpoly_struct));

        ord = mpoly_ordering_randtest(state);

        nvars = n_randint(state, 10) + 1;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);
        nmod_mpoly_ctx_init(ctx, nvars, ord, modulus);

        for (w = 0; w < num; w++)
        {
            nmod_mpoly_init(g + w, ctx);
            darr[w] = g + w;

            nmod_mpoly_init(q + w, ctx);
            qarr[w] = q + w;
        }  

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(k1, ctx);
        nmod_mpoly_init(k2, ctx);
        nmod_mpoly_init(r, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 10) + 1;

        maxbits = 14/(nvars + mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1;
        exp_bits = n_randint(state, maxbits) + 1;
        exp_bits1 = n_randint(state, maxbits) + 1;
        exp_bits2 = n_randint(state, maxbits) + 1;

        exp_bound = n_randbits(state, exp_bits);
        exp_bound1 = n_randbits(state, exp_bits1);
        exp_bound2 = n_randbits(state, exp_bits2);

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest(f, state, len1, exp_bound1, ctx);
            for (w = 0; w < num; w++)
            {
                do {
                    nmod_mpoly_randtest(darr[w], state, len2, exp_bound2 + 1, ctx);
                } while (darr[w]->length == 0);
                nmod_mpoly_randtest(qarr[w], state, len, exp_bound, ctx);
            }
            nmod_mpoly_randtest(k1, state, len, exp_bound, ctx);
            nmod_mpoly_randtest(k2, state, len, exp_bound, ctx);

            nmod_mpoly_divrem_ideal_monagan_pearce(qarr, r, f, darr, num, ctx);
            nmod_mpoly_test(r, ctx);
            for (w = 0; w < num; w++)
            {
                nmod_mpoly_test(qarr[w], ctx);
                nmod_mpoly_remainder_strongtest(r, darr[w], ctx);
            }

            nmod_mpoly_zero(k2, ctx);
            for (w = 0; w < num; w++)
            {
                nmod_mpoly_mul_johnson(k1, qarr[w], darr[w], ctx);
                nmod_mpoly_add(k2, k2, k1, ctx);
            }
            nmod_mpoly_add(k2, k2, r, ctx);

            result = nmod_mpoly_equal(f, k2, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f = g1*q1 + ... + gn*qn + r for random polys\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
        }

        for (w = 0; w < num; w++)
            nmod_mpoly_clear(qarr[w], ctx);
        for (w = 0; w < num; w++)
            nmod_mpoly_clear(darr[w], ctx);
        nmod_mpoly_clear(f, ctx);  
        nmod_mpoly_clear(k1, ctx);  
        nmod_mpoly_clear(k2, ctx);  
        nmod_mpoly_clear(r, ctx);

        flint_free(g);
        flint_free(q);  
        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, r, k1, k2;
        nmod_mpoly_struct * g, * q;
        ordering_t ord;
        mp_limb_t modulus;
        slong maxbits;
        slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2, num;
        slong exp_bits, exp_bits1, exp_bits2;
        nmod_mpoly_struct * qarr[5], * darr[5];

        num = n_randint(state, 5) + 1;

        g = (nmod_mpoly_struct *) flint_malloc(num*sizeof(nmod_mpoly_struct));
        q = (nmod_mpoly_struct *) flint_malloc(num*sizeof(nmod_mpoly_struct));

        ord = mpoly_ordering_randtest(state);

        nvars = n_randint(state, 10) + 1;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);
        nmod_mpoly_ctx_init(ctx, nvars, ord, modulus);

        for (w = 0; w < num; w++)
        {
            nmod_mpoly_init(g + w, ctx);
            darr[w] = g + w;

            nmod_mpoly_init(q + w, ctx);
            qarr[w] = q + w;
        }  

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(k1, ctx);
        nmod_mpoly_init(k2, ctx);
        nmod_mpoly_init(r, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 10) + 1;

        maxbits = 14/(nvars + mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1;
        exp_bits = n_randint(state, maxbits) + 1;
        exp_bits1 = n_randint(state, maxbits) + 1;
        exp_bits2 = n_randint(state, maxbits) + 1;

        exp_bound = n_randbits(state, exp_bits);
        exp_bound1 = n_randbits(state, exp_bits1);
        exp_bound2 = n_randbits(state, exp_bits2);

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest(f, state, len1, exp_bound1, ctx);
            for (w = 0; w < num; w++)
            {
                do {
                    nmod_mpoly_randtest(darr[w], state, len2, exp_bound2 + 1, ctx);
                } while (darr[w]->length == 0);
                nmod_mpoly_randtest(qarr[w], state, len, exp_bound, ctx);
            }
            nmod_mpoly_randtest(k1, state, len, exp_bound, ctx);
            nmod_mpoly_randtest(k2, state, len, exp_bound, ctx);

            nmod_mpoly_set(r, f, ctx);

            nmod_mpoly_divrem_ideal_monagan_pearce(qarr, f, f, darr, num, ctx);
            nmod_mpoly_test(f, ctx);
            for (w = 0; w < num; w++)
            {
                nmod_mpoly_test(qarr[w], ctx);
                nmod_mpoly_remainder_strongtest(f, darr[w], ctx);
            }

            nmod_mpoly_zero(k2, ctx);
            for (w = 0; w < num; w++)
            {
                nmod_mpoly_mul_johnson(k1, qarr[w], darr[w], ctx);
                nmod_mpoly_add(k2, k2, k1, ctx);
            }
            nmod_mpoly_add(k2, k2, f, ctx);

            result = nmod_mpoly_equal(r, k2, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
        }

        for (w = 0; w < num; w++)
            nmod_mpoly_clear(qarr[w], ctx);
        for (w = 0; w < num; w++)
            nmod_mpoly_clear(darr[w], ctx);
        nmod_mpoly_clear(f, ctx);  
        nmod_mpoly_clear(k1, ctx);  
        nmod_mpoly_clear(k2, ctx);  
        nmod_mpoly_clear(r, ctx);

        flint_free(g);
        flint_free(q);  
        nmod_mpoly_ctx_clear(ctx);
    }

    printf("PASS\n");
    FLINT_TEST_CLEANUP(state);

    return 0;
}

