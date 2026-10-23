/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "gr.h"
#include "mpn_mod.h"
#include "aprcl.h"

/* checks that the mpn element f equals the fmpz_mod element g, which
   must be reduced modulo the cyclotomic polynomial */
static int
_unity_zp_mpn_equals(const unity_zp_mpn f, const unity_zp g)
{
    slong i;
    slong n = MPN_MOD_CTX_NLIMBS(f->ctx);
    int result = 1;
    fmpz_t c1, c2;

    fmpz_init(c1);
    fmpz_init(c2);

    for (i = 0; i < (slong) f->d; i++)
    {
        GR_MUST_SUCCEED(mpn_mod_get_fmpz(c1, f->coeffs + i * n, f->ctx));
        fmpz_mod_poly_get_coeff_fmpz(c2, g->poly, i, g->ctx);

        if (!fmpz_equal(c1, c2))
        {
            result = 0;
            break;
        }
    }

    fmpz_clear(c1);
    fmpz_clear(c2);

    return result;
}

TEST_FUNCTION_START(aprcl_unity_zp_mpn, state)
{
    int i;

    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        ulong p, exp, pk, ind, x, pow;
        slong j, h1, h2;
        fmpz_t n, c, upow;
        unity_zp f, g, h, temp;
        unity_zp_mpn F, G, H;
        gr_ctx_t ctx;

        static const ulong pk_tab[][2] = {
            {2, 1}, {2, 2}, {2, 3}, {2, 4}, {2, 5},
            {3, 1}, {3, 2}, {3, 3},
            {5, 1}, {5, 2}, {7, 1}, {11, 1}, {13, 1}
        };

        j = n_randint(state, sizeof(pk_tab) / sizeof(pk_tab[0]));
        p = pk_tab[j][0];
        exp = pk_tab[j][1];
        pk = n_pow(p, exp);

        /* random modulus of 2 to 16 limbs */
        fmpz_init(n);
        if (n_randint(state, 4) == 0)
        {
            /* moduli of exactly FLINT_BITS (nlimbs - 1) + 1 bits, exercising
               boundary cases in the accumulator size logic */
            slong nl = 2 + n_randint(state, MPN_MOD_MAX_LIMBS - 1);
            fmpz_randbits(n, state, 3 + n_randint(state, 40));
            fmpz_abs(n, n);
            fmpz_setbit(n, FLINT_BITS * (nl - 1));
            fmpz_clrbit(n, 0);
            fmpz_setbit(n, 1);   /* keep n odd or even at random */
            if (n_randint(state, 2))
                fmpz_setbit(n, 0);
        }
        else
        {
            fmpz_randbits(n, state, FLINT_BITS + 1
                + n_randint(state, FLINT_BITS * (MPN_MOD_MAX_LIMBS - 1)));
            fmpz_abs(n, n);
            fmpz_setbit(n, FLINT_BITS);   /* ensure at least 2 limbs */
        }

        if (gr_ctx_init_mpn_mod(ctx, n) != GR_SUCCESS)
        {
            flint_printf("FAIL: cannot create mpn_mod ctx\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_init(c);
        fmpz_init(upow);

        unity_zp_init(f, p, exp, n);
        unity_zp_init(g, p, exp, n);
        unity_zp_init(h, p, exp, n);
        unity_zp_init(temp, p, exp, n);

        unity_zp_mpn_init(F, p, exp, ctx);
        unity_zp_mpn_init(G, p, exp, ctx);
        unity_zp_mpn_init(H, p, exp, ctx);

        /* random g, h with coefficients at arbitrary indices < p^k */
        for (j = 0; j < (slong) pk; j++)
        {
            fmpz_randtest_mod(c, state, n);
            unity_zp_coeff_set_fmpz(g, j, c);
            fmpz_randtest_mod(c, state, n);
            unity_zp_coeff_set_fmpz(h, j, c);
        }

        /* conversion must agree with cyclotomic reduction */
        unity_zp_mpn_set_unity_zp(G, g);
        unity_zp_mpn_set_unity_zp(H, h);

        unity_zp_copy(temp, g);
        _unity_zp_reduce_cyclotomic(temp);

        if (!_unity_zp_mpn_equals(G, temp))
            TEST_FUNCTION_FAIL("conversion, p = %wu, exp = %wu\n", p, exp);

        /* reduce g, h so that fmpz and mpn versions match */
        _unity_zp_reduce_cyclotomic(g);
        _unity_zp_reduce_cyclotomic(h);

        /* multiplication */
        unity_zp_mul(f, g, h);
        unity_zp_mpn_mul(F, G, H);

        if (!_unity_zp_mpn_equals(F, f))
            TEST_FUNCTION_FAIL("mul, p = %wu, exp = %wu\n", p, exp);

        /* aliased multiplication */
        unity_zp_mpn_copy(F, G);
        unity_zp_mpn_mul(F, F, H);
        if (!_unity_zp_mpn_equals(F, f))
            TEST_FUNCTION_FAIL("aliased mul, p = %wu, exp = %wu\n", p, exp);

        /* squaring */
        unity_zp_sqr(f, g);
        unity_zp_mpn_sqr(F, G);

        if (!_unity_zp_mpn_equals(F, f))
            TEST_FUNCTION_FAIL("sqr, p = %wu, exp = %wu\n", p, exp);

        /* scalar multiplication */
        x = n_randtest(state);
        unity_zp_mul_scalar_ui(f, g, x);
        unity_zp_mpn_mul_scalar_ui(F, G, x);

        if (!_unity_zp_mpn_equals(F, f))
            TEST_FUNCTION_FAIL("mul_scalar_ui, p = %wu, exp = %wu\n", p, exp);

        /* pow_ui */
        pow = n_randint(state, 100);
        unity_zp_pow_ui(f, g, pow);
        _unity_zp_reduce_cyclotomic(f);
        unity_zp_mpn_pow_ui(F, G, pow);

        if (!_unity_zp_mpn_equals(F, f))
            TEST_FUNCTION_FAIL("pow_ui, p = %wu, exp = %wu, pow = %wu\n",
                    p, exp, pow);

        /* pow_sliding_fmpz */
        fmpz_randtest_unsigned(upow, state, 300);
        fmpz_add_ui(upow, upow, 1);
        unity_zp_pow_sliding_fmpz(f, g, upow);
        _unity_zp_reduce_cyclotomic(f);
        unity_zp_mpn_pow_sliding_fmpz(F, G, upow);

        if (!_unity_zp_mpn_equals(F, f))
            TEST_FUNCTION_FAIL("pow_sliding_fmpz, p = %wu, exp = %wu\n",
                    p, exp);

        /* aut_inv (x must be coprime to p) */
        if (pk > 2)
        {
            do
                x = 1 + n_randint(state, pk - 1);
            while (x % p == 0);

            unity_zp_aut_inv(f, g, x);
            unity_zp_mpn_aut_inv(F, G, x);

            if (!_unity_zp_mpn_equals(F, f))
                TEST_FUNCTION_FAIL("aut_inv, p = %wu, exp = %wu, x = %wu\n",
                        p, exp, x);
        }

        /* is_unity on a root of unity */
        ind = n_randint(state, pk);
        unity_zp_set_zero(f);
        unity_zp_coeff_set_ui(f, ind, 1);
        unity_zp_mpn_set_unity_zp(F, f);

        h1 = unity_zp_is_unity(f);
        h2 = unity_zp_mpn_is_unity(F);

        if (h1 != h2 || h1 != (slong) ind)
            TEST_FUNCTION_FAIL("is_unity(zeta^%wu) = %wd vs %wd, "
                    "p = %wu, exp = %wu\n", ind, h1, h2, p, exp);

        /* is_unity on a random element */
        h1 = unity_zp_is_unity(g);
        h2 = unity_zp_mpn_is_unity(G);

        if (h1 != h2)
            TEST_FUNCTION_FAIL("is_unity random, %wd vs %wd, "
                    "p = %wu, exp = %wu\n", h1, h2, p, exp);

        unity_zp_mpn_clear(F);
        unity_zp_mpn_clear(G);
        unity_zp_mpn_clear(H);

        unity_zp_clear(f);
        unity_zp_clear(g);
        unity_zp_clear(h);
        unity_zp_clear(temp);

        gr_ctx_clear(ctx);

        fmpz_clear(n);
        fmpz_clear(c);
        fmpz_clear(upow);
    }

    TEST_FUNCTION_END(state);
}
