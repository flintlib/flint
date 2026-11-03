/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"

/*
    The general high products flint_mpn_mulhigh / flint_mpn_sqrhigh and their
    normalised variants, for all shapes: the computed n + 1 limbs (result and
    guard limb) must never exceed the exact ones, and fall short by at most
    the documented bound.
*/

/* Pick a length: mostly small, sometimes up to the FFT range. */
static mp_size_t
random_length(flint_rand_t state, mp_size_t maxn)
{
    switch (n_randint(state, 8))
    {
        case 0: return 1 + n_randint(state, maxn);
        case 1: return 1 + n_randint(state, 100);
        case 2: return 1 + n_randint(state, 40);
        default: return 1 + n_randint(state, 20);
    }
}

/* Inputs maximising the dropped partial products of the kernels (as in
   t-mulhigh_n.c), stretched over any length. */
static void
adversarial_general(mp_ptr xp, mp_size_t xn, mp_ptr yp, mp_size_t yn, flint_rand_t state)
{
    mp_size_t i;
    mp_limb_t a = -(UWORD(1) << (FLINT_BITS / 2)) - 1, b = -(UWORD(1) << (FLINT_BITS / 2)) + 1;

    for (i = 0; i < xn; i++)
        xp[i] = (yp == xp && (i & 1)) ? b : a;
    if (yp != xp)
        for (i = 0; i < yn; i++)
            yp[i] = b;

    if (n_randint(state, 2))
    {
        xp[xn - 1] = n_randtest(state) | (UWORD(1) << (FLINT_BITS - 1));
        if (yp != xp)
            yp[yn - 1] = n_randtest(state) | (UWORD(1) << (FLINT_BITS - 1));
    }
}

TEST_FUNCTION_START(flint_mpn_mulhigh, state)
{
    slong ix;
    mp_size_t maxn;

#if FLINT_HAVE_FFT_SMALL
    maxn = FLINT_MAX(FLINT_MPN_MULHIGH_FFT_SMALL_CUTOFF, FLINT_MPN_SQRHIGH_FFT_SMALL_CUTOFF) + 200;
#else
    maxn = FLINT_MPN_MULHIGH_MULDERS_CUTOFF + 200;
#endif

    for (ix = 0; ix < 3000 * flint_test_multiplier(); ix++)
    {
        mp_size_t xn, yn, n, zn, d, k;
        mp_ptr xp, yp, full, rp, exact, abuf = NULL, rp_save = NULL, xp_orig, yp_orig;
        mp_limb_t g, bound, deficit;
        int alias;
        int squaring, normalised, shifted = 0;
        mp_limb_pair_t pair;

        squaring = n_randint(state, 2);
        normalised = n_randint(state, 2);

        xn = random_length(state, maxn);
        yn = squaring ? xn : random_length(state, maxn);

        /* output lengths around the operand lengths, sometimes anything */
        switch (n_randint(state, 4))
        {
            case 0: n = 1 + n_randint(state, xn + yn + 2); break;
            case 1: n = FLINT_MIN(xn, yn) - (slong) n_randint(state, 3); break;
            case 2: n = FLINT_MAX(xn, yn) - (slong) n_randint(state, 3); break;
            default: n = 1 + n_randint(state, FLINT_MAX(xn, yn)); break;
        }
        n = FLINT_MAX(n, 1);

        zn = xn + yn;

        xp = flint_malloc(sizeof(mp_limb_t) * xn);
        yp = squaring ? xp : flint_malloc(sizeof(mp_limb_t) * yn);
        xp_orig = xp;
        yp_orig = yp;
        full = flint_malloc(sizeof(mp_limb_t) * (zn + 1));
        exact = flint_malloc(sizeof(mp_limb_t) * (n + 2));
        rp = flint_malloc(sizeof(mp_limb_t) * (n + 1));
        rp[n] = UWORD(0xdeadbeef);

        if (n_randint(state, 4) == 0)
        {
            adversarial_general(xp, xn, yp, yn, state);
        }
        else
        {
            flint_mpn_rrandom(xp, state, xn);
            if (!squaring)
                flint_mpn_rrandom(yp, state, yn);
        }
        if (n_randint(state, 2))
        {
            xp[xn - 1] |= (UWORD(1) << (FLINT_BITS - 1));
            yp[yn - 1] |= (UWORD(1) << (FLINT_BITS - 1));
        }

        if (squaring)
            flint_mpn_sqr(full, xp, xn);
        else if (xn >= yn)
            flint_mpn_mul(full, xp, xn, yp, yn);
        else
            flint_mpn_mul(full, yp, yn, xp, xn);

        d = zn - n - 1;

        /* the output may alias an operand (whose buffer must then hold
           n limbs): compute into a copy of the operand's buffer */
        alias = n_randint(state, 4);
        if (alias == 1 || (alias == 2 && !squaring))
        {
            mp_size_t len = (alias == 1) ? xn : yn;
            mp_srcptr src = (alias == 1) ? xp : yp;
            abuf = flint_malloc(sizeof(mp_limb_t) * FLINT_MAX(len, n + 1));
            flint_mpn_copyi(abuf, src, len);
            rp_save = rp;
            rp = abuf;
            if (alias == 1) { xp = abuf; if (squaring) yp = abuf; }
            else yp = abuf;
        }
        else
            alias = 0;

        if (!normalised)
        {
            if (squaring)
                g = flint_mpn_sqrhigh(rp, xp, xn, n);
            else
                g = flint_mpn_mulhigh(rp, xp, xn, yp, yn, n);
            bound = squaring ? FLINT_MPN_SQRHIGH_GUARD_ERR_BOUND(n)
                             : FLINT_MPN_MULHIGH_GUARD_ERR_BOUND(n);
        }
        else
        {
            if (squaring)
                pair = flint_mpn_sqrhigh_normalised(rp, xp, xn, n);
            else
                pair = flint_mpn_mulhigh_normalised(rp, xp, xn, yp, yn, n);
            g = pair.m1;
            shifted = pair.m2;
            bound = squaring ? FLINT_MPN_SQRHIGH_NORMALISED_GUARD_ERR_BOUND(n)
                             : FLINT_MPN_MULHIGH_NORMALISED_GUARD_ERR_BOUND(n);

            if (shifted != 1 && shifted != 0)
                TEST_FUNCTION_FAIL("shift flag = %d\n", shifted);

            /* the exact product shifted by the same amount (as in
               t-mulhigh_n_normalised.c, an exact value whose top bit the
               deficit hid is not provoked here) */
            if (shifted)
            {
                mpn_lshift(full, full, zn, 1);
            }
            else if ((rp[n - 1] >> (FLINT_BITS - 1)) == 0)
            {
                TEST_FUNCTION_FAIL("not normalised\n"
                    "squaring = %d, xn = %wd, yn = %wd, n = %wd\n",
                    squaring, xn, yn, n);
            }
        }

        /* exact n + 1 limbs: limbs [d, zn) of the product, zero padded */
        for (k = 0; k <= n; k++)
        {
            if (d + k >= 0 && d + k < zn)
                exact[k] = full[d + k];
            else
                exact[k] = 0;
        }

        if (alias)
        {
            /* restore the operand pointers for the exact product below */
            if (alias == 1) { xp = xp_orig; if (squaring) yp = xp_orig; }
            else yp = yp_orig;
        }
        else if (rp[n] != UWORD(0xdeadbeef))
            TEST_FUNCTION_FAIL("wrote past the output\n"
                "squaring = %d, xn = %wd, yn = %wd, n = %wd\n",
                squaring, xn, yn, n);

        /* deficit = exact - (rp, g), must be in [0, bound] */
        rp[n] = 0;
        {
            mp_ptr e = exact;
            mp_limb_t borrow;

            /* subtract g from the guard position then the n limbs */
            borrow = mpn_sub_1(e, e, n + 1, g);
            borrow |= mpn_sub(e + 1, e + 1, n, rp, n);

            if (borrow)
                TEST_FUNCTION_FAIL("computed value exceeds the exact one\n"
                    "squaring = %d, normalised = %d, xn = %wd, yn = %wd, n = %wd, shifted = %d\n",
                    squaring, normalised, xn, yn, n, shifted);

            for (k = 1; k <= n; k++)
                if (e[k] != 0)
                    TEST_FUNCTION_FAIL("deficit of at least a limb\n"
                        "squaring = %d, normalised = %d, xn = %wd, yn = %wd, n = %wd, shifted = %d\n",
                        squaring, normalised, xn, yn, n, shifted);

            deficit = e[0];
            if (deficit > bound)
                TEST_FUNCTION_FAIL("deficit %wu exceeds the bound %wu\n"
                    "squaring = %d, normalised = %d, xn = %wd, yn = %wd, n = %wd, shifted = %d\n",
                    deficit, bound, squaring, normalised, xn, yn, n, shifted);
        }

        if (alias)
        {
            flint_free(abuf);
            rp = rp_save;
        }
        flint_free(xp);
        if (!squaring)
            flint_free(yp);
        flint_free(full);
        flint_free(exact);
        flint_free(rp);
    }

    TEST_FUNCTION_END(state);
}
