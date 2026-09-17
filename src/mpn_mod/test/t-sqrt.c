/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_mod.h"
#include "fmpz.h"
#include "ulong_extras.h"

/*
    A prime of about the given number of bits. With val2 > 0 the prime is
    congruent to 1 modulo 2^val2, which is what exercises the descent in
    Tonelli and Shanks; val2 = 1 covers the p = 3 mod 4 shortcut instead.
*/
static void
_random_prime(fmpz_t p, flint_rand_t state, flint_bitcnt_t bits, slong val2)
{
    fmpz_t q;

    fmpz_init(q);

    while (1)
    {
        if (val2 <= 0)
        {
            fmpz_randprime(p, state, bits, 0);
            break;
        }

        /* p = q 2^val2 + 1 with q odd */
        fmpz_randbits(q, state, FLINT_MAX(bits - val2, 2));
        fmpz_abs(q, q);
        fmpz_setbit(q, 0);
        fmpz_mul_2exp(p, q, val2);
        fmpz_add_ui(p, p, 1);

        if (fmpz_bits(p) >= MPN_MOD_MIN_LIMBS * FLINT_BITS
                && fmpz_size(p) <= MPN_MOD_MAX_LIMBS
                && fmpz_is_probabprime(p))
            break;
    }

    fmpz_clear(q);
}

TEST_FUNCTION_START(mpn_mod_sqrt, state)
{
    slong iter;

    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx, ref;
        fmpz_t p, x;
        nn_ptr a, r, t;
        slong nlimbs, j, val2;
        truth_t sq, sq_ref;
        int status;

        /* a prime of 2 to 5 limbs, with a range of 2-adic valuations */
        nlimbs = 2 + n_randint(state, 4);
        val2 = (slong) n_randint(state, 12);

        fmpz_init(p);
        fmpz_init(x);

        _random_prime(p, state, nlimbs * FLINT_BITS, val2);

        if (gr_ctx_init_mpn_mod(ctx, p) != GR_SUCCESS)
        {
            fmpz_clear(p);
            fmpz_clear(x);
            continue;
        }

        gr_ctx_init_fmpz_mod(ref, p);

        /* the modulus is prime, and neither ring knows that by itself */
        GR_MUST_SUCCEED(gr_ctx_set_is_field(ctx, T_TRUE));
        GR_MUST_SUCCEED(gr_ctx_set_is_field(ref, T_TRUE));

        nlimbs = MPN_MOD_CTX_NLIMBS(ctx);
        a = flint_malloc(3 * nlimbs * sizeof(ulong));
        r = a + nlimbs;
        t = r + nlimbs;

        for (j = 0; j < 8; j++)
        {
            /* zero and one now and then, otherwise uniform; and squares
               often enough that the successful branch is well covered */
            if (j == 0)
                fmpz_zero(x);
            else if (j == 1)
                fmpz_one(x);
            else
            {
                fmpz_randm(x, state, p);

                if (j % 2 == 0)     /* force a square */
                {
                    fmpz_mul(x, x, x);
                    fmpz_mod(x, x, p);
                }
            }

            GR_MUST_SUCCEED(mpn_mod_set_fmpz(a, x, ctx));

            /* is_square must agree with fmpz_mod */
            sq = mpn_mod_is_square(a, ctx);
            {
                fmpz_t y;
                fmpz_init_set(y, x);
                sq_ref = gr_is_square(y, ref);
                fmpz_clear(y);
            }

            if (sq != sq_ref)
            {
                flint_printf("FAIL: is_square disagrees with fmpz_mod\n");
                flint_printf("p = "); fmpz_print(p); flint_printf("\n");
                flint_printf("x = "); fmpz_print(x); flint_printf("\n");
                flint_printf("mpn_mod %d, fmpz_mod %d\n", (int) sq, (int) sq_ref);
                fflush(stdout);
                flint_abort();
            }

            status = mpn_mod_sqrt(r, a, ctx);

            /* a root exists exactly when x is a square */
            if ((status == GR_SUCCESS) != (sq == T_TRUE))
            {
                flint_printf("FAIL: sqrt status %d but is_square %d\n",
                        status, (int) sq);
                flint_printf("p = "); fmpz_print(p); flint_printf("\n");
                flint_printf("x = "); fmpz_print(x); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            /* and when it does, it must square back to x */
            if (status == GR_SUCCESS)
            {
                GR_MUST_SUCCEED(mpn_mod_sqr(t, r, ctx));

                if (mpn_mod_equal(t, a, ctx) != T_TRUE)
                {
                    flint_printf("FAIL: sqrt(x)^2 != x\n");
                    flint_printf("p = "); fmpz_print(p); flint_printf("\n");
                    flint_printf("x = "); fmpz_print(x); flint_printf("\n");
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        flint_free(a);
        gr_ctx_clear(ref);
        gr_ctx_clear(ctx);
        fmpz_clear(p);
        fmpz_clear(x);
    }

    /* without a prime modulus neither question can be answered */
    {
        gr_ctx_t ctx;
        fmpz_t p;
        nn_ptr a, r;
        slong nlimbs;

        fmpz_init(p);
        fmpz_set_ui(p, 1);
        fmpz_mul_2exp(p, p, 130);       /* 2^130, very much not prime */

        if (gr_ctx_init_mpn_mod(ctx, p) == GR_SUCCESS)
        {
            nlimbs = MPN_MOD_CTX_NLIMBS(ctx);
            a = flint_malloc(2 * nlimbs * sizeof(ulong));
            r = a + nlimbs;

            GR_MUST_SUCCEED(mpn_mod_set_ui(a, 7, ctx));

            FLINT_TEST(mpn_mod_is_square(a, ctx) == T_UNKNOWN);
            FLINT_TEST(mpn_mod_sqrt(r, a, ctx) == GR_UNABLE);

            /* zero and one are squares whatever the modulus */
            GR_MUST_SUCCEED(mpn_mod_zero(a, ctx));
            FLINT_TEST(mpn_mod_is_square(a, ctx) == T_TRUE);
            FLINT_TEST(mpn_mod_sqrt(r, a, ctx) == GR_SUCCESS);

            flint_free(a);
            gr_ctx_clear(ctx);
        }

        fmpz_clear(p);
    }

    TEST_FUNCTION_END(state);
}
