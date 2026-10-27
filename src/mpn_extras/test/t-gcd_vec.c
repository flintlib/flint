/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "test_helpers.h"
#include "ulong_extras.h"
#include "mpn_extras.h"
#include "gmpcompat.h"

/* random value of the form c * (random of about bits bits), as limbs */
static mp_size_t
_rand_mul(mp_ptr x, flint_rand_t state, const mpz_t c, flint_bitcnt_t bits)
{
    mpz_t t;
    mp_size_t n;

    mpz_init(t);
    if (bits == 0)
        mpz_set_ui(t, 0);
    else if (n_randint(state, 2))
    {
        mp_size_t k = (bits + FLINT_BITS - 1) / FLINT_BITS;
        flint_mpn_rrandom(mpz_limbs_write(t, k), state, k);
        mpz_limbs_finish(t, k);
        mpz_fdiv_r_2exp(t, t, bits);
    }
    else
    {
        mp_size_t k = (bits + FLINT_BITS - 1) / FLINT_BITS;
        flint_mpn_urandomb(mpz_limbs_write(t, k), state, bits);
        mpz_limbs_finish(t, k);
    }

    /* occasionally a very small cofactor */
    if (n_randint(state, 8) == 0)
        mpz_set_ui(t, n_randint(state, 4));

    mpz_mul(t, t, c);
    n = mpz_size(t);
    if (n != 0)
        flint_mpn_copyi(x, mpz_limbs_read(t), n);
    mpz_clear(t);
    return n;
}

static void
_rand_content(mpz_t c, flint_rand_t state)
{
    switch (n_randint(state, 6))
    {
        case 0:
            mpz_set_ui(c, 1);
            break;
        case 1:
            flint_mpz_set_ui(c, n_randtest_not_zero(state));
            break;
        case 2:
            mpz_set_ui(c, 1);
            mpz_mul_2exp(c, c, n_randint(state, 200));
            break;
        case 3:
            mpz_set_ui(c, n_randint(state, 100) + 1);
            mpz_mul_2exp(c, c, n_randint(state, 200));
            break;
        default:
        {
            mp_size_t k = 1 + n_randint(state, 6);
            flint_mpn_rrandom(mpz_limbs_write(c, k), state, k);
            mpz_limbs_finish(c, k);
            if (mpz_sgn(c) == 0)
                mpz_set_ui(c, 1);
            mpz_mul_2exp(c, c, n_randint(state, 3) * n_randint(state, 70));
        }
    }
}

TEST_FUNCTION_START(flint_mpn_gcd_vec, state)
{
    slong iter;

    /* single limbs */
    for (iter = 0; iter < 20000 * flint_test_multiplier(); iter++)
    {
        slong len = n_randint(state, 4) == 0 ? n_randint(state, 300) : n_randint(state, 20);
        mp_ptr x = flint_malloc((len + 1) * sizeof(mp_limb_t));
        ulong c, g, h, ref;
        slong i;

        c = n_randtest_not_zero(state);
        if (n_randint(state, 4) == 0)
            c = UWORD(1) << n_randint(state, FLINT_BITS);

        for (i = 0; i < len; i++)
        {
            ulong m = UWORD_MAX / c, r;
            r = n_randtest(state);
            if (m != UWORD_MAX)
                r %= (m + 1);
            switch (n_randint(state, 4))
            {
                case 0: x[i] = 0; break;
                case 1: x[i] = c * n_randint(state, FLINT_MIN(m, 10) + 1); break;
                default: x[i] = c * r; break;
            }
            /* occasionally break the content */
            if (n_randint(state, 50) == 0)
                x[i] = n_randtest(state);
        }

        g = n_randint(state, 3) == 0 ? c * n_randint(state, 5) : 0;

        ref = g;
        for (i = 0; i < len; i++)
            ref = n_gcd(ref, x[i]);

        h = flint_mpn_gcd_vec_1(x, len, g);

        if (h != ref)
            TEST_FUNCTION_FAIL("single limbs: len = %wd, g = %wu, h = %wu, ref = %wu\n", len, g, h, ref);

        flint_free(x);
    }

    /* multi-limb values */
    for (iter = 0; iter < 3000 * flint_test_multiplier(); iter++)
    {
        slong len, i;
        mpz_t c, ref, t;
        mp_ptr * xs;
        mp_srcptr * xp;
        mp_size_t * xn;
        mp_ptr g;
        mp_size_t gn, gn0, gmax, n;
        flint_bitcnt_t maxbits;

        len = n_randint(state, 4) == 0 ? n_randint(state, 200) : n_randint(state, 12);
        maxbits = n_randint(state, 5) == 0 ? (n_randint(state, 2) ? 6000 : 3000) : 400;

        mpz_init(c);
        mpz_init(ref);
        mpz_init(t);
        _rand_content(c, state);

        xs = flint_malloc((len + 1) * sizeof(mp_ptr));
        xp = flint_malloc((len + 1) * sizeof(mp_srcptr));
        xn = flint_malloc((len + 1) * sizeof(mp_size_t));

        gmax = 0;
        for (i = 0; i < len; i++)
        {
            flint_bitcnt_t bits = n_randint(state, maxbits + 1);
            if (n_randint(state, 4) == 0)
                bits = n_randint(state, 70);
            xs[i] = flint_malloc((mpz_size(c) + bits / FLINT_BITS + 2) * sizeof(mp_limb_t));
            xn[i] = _rand_mul(xs[i], state, c, bits);
            /* occasionally break the content */
            if (n_randint(state, 30) == 0)
            {
                mpz_t one;
                mpz_init_set_ui(one, 1);
                xn[i] = _rand_mul(xs[i], state, one, bits);
                mpz_clear(one);
            }
            xp[i] = xs[i];
            gmax = FLINT_MAX(gmax, xn[i]);
        }

        /* initial g */
        g = flint_malloc((gmax + mpz_size(c) + 10) * sizeof(mp_limb_t));
        gn0 = 0;
        if (n_randint(state, 3) == 0)
            gn0 = _rand_mul(g, state, c, n_randint(state, 300));

        mpz_set_ui(ref, 0);
        if (gn0 != 0)
            mpz_import(ref, gn0, -1, sizeof(mp_limb_t), 0, 0, g);
        for (i = 0; i < len; i++)
        {
            if (xn[i] != 0)
            {
                mpz_import(t, xn[i], -1, sizeof(mp_limb_t), 0, 0, xp[i]);
                mpz_gcd(ref, ref, t);
            }
        }

        gn = flint_mpn_gcd_vec(g, gn0, xp, xn, len);

        n = mpz_size(ref);
        if (gn != n || (n != 0 && !flint_mpn_equal_p(g, mpz_limbs_read(ref), n)))
        {
            TEST_FUNCTION_FAIL("multi-limb: len = %wd, gn0 = %wd, gn = %wd, n = %wd\n", len, (slong) gn0, (slong) gn, (slong) n);
        }

        for (i = 0; i < len; i++)
            flint_free(xs[i]);
        flint_free(xs);
        flint_free(xp);
        flint_free(xn);
        flint_free(g);
        mpz_clear(c);
        mpz_clear(ref);
        mpz_clear(t);
    }

    /* targeted shapes: two-limb content over two- and three-limb values,
       and single-limb content over longer values, with near misses */
    for (iter = 0; iter < 3000 * flint_test_multiplier(); iter++)
    {
        slong len, i;
        int two = n_randint(state, 2);
        mpz_t c, ref, t;
        mp_limb_t * X;
        mp_srcptr * xp;
        mp_size_t * xn;
        mp_limb_t g[200];
        mp_size_t gn, n, maxn;

        len = 1 + n_randint(state, 40);
        maxn = two ? 3 : 2 + n_randint(state, 100);

        mpz_init(c);
        mpz_init(ref);
        mpz_init(t);

        if (two)
        {
            flint_mpn_rrandom(mpz_limbs_write(c, 2), state, 2);
            mpz_limbs_finish(c, 2);
            mpz_setbit(c, FLINT_BITS + n_randint(state, FLINT_BITS));
            if (n_randint(state, 2))
                mpz_setbit(c, 0);
        }
        else
        {
            flint_mpz_set_ui(c, n_randtest_not_zero(state));
        }

        X = flint_malloc(len * (maxn + 3) * sizeof(mp_limb_t));
        xp = flint_malloc(len * sizeof(mp_srcptr));
        xn = flint_malloc(len * sizeof(mp_size_t));

        mpz_set_ui(ref, 0);
        for (i = 0; i < len; i++)
        {
            mp_size_t qn = 1 + n_randint(state, two ? 1 : maxn);
            flint_mpn_rrandom(mpz_limbs_write(t, qn), state, qn);
            mpz_limbs_finish(t, qn);
            mpz_mul(t, t, c);
            if (n_randint(state, 20) == 0)
            {
                if (n_randint(state, 2))
                    mpz_add_ui(t, t, 1 + n_randint(state, 3));
                else
                    mpz_add(t, t, c), mpz_sub_ui(t, t, 1);
            }
            while (mpz_size(t) > (two ? 3 : (size_t) maxn))
                mpz_tdiv_q_2exp(t, t, FLINT_BITS);
            n = mpz_size(t);
            if (n != 0)
                flint_mpn_copyi(X + i * (maxn + 3), mpz_limbs_read(t), n);
            xp[i] = X + i * (maxn + 3);
            xn[i] = n;
            mpz_gcd(ref, ref, t);
        }

        gn = flint_mpn_gcd_vec(g, 0, xp, xn, len);
        n = mpz_size(ref);
        if (gn != n || (n != 0 && !flint_mpn_equal_p(g, mpz_limbs_read(ref), n)))
            TEST_FUNCTION_FAIL("targeted: two = %d, len = %wd, gn = %wd, n = %wd\n", two, len, (slong) gn, (slong) n);

        flint_free(X);
        flint_free(xp);
        flint_free(xn);
        mpz_clear(c);
        mpz_clear(ref);
        mpz_clear(t);
    }

    TEST_FUNCTION_END(state);
}
