/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod.h"

/* NMOD_RED2_BOUNDED (and NMOD_RED3_BOUNDED, below) against n_ll_mod_preinv,
   on two-word values below
   2^(nu + FLINT_BITS - 2): random values, values at and next to multiples
   of n (the quotient estimate is one too small there), values at the
   bound, products of residues and sums of them */
TEST_FUNCTION_START(nmod_red2_bounded, state)
{
    slong i;

    for (i = 0; i < 100000 * flint_test_multiplier(); i++)
    {
        nmod_t mod;
        ulong n, ninv, hi, lo, r, s, b_hi;
        flint_bitcnt_t nu = 3 + n_randint(state, FLINT_BITS - 3);

        /* n of nu bits, sometimes 2^(nu-1), 2^(nu-1) + 1 or 2^nu - 1 */
        switch (n_randint(state, 8))
        {
            case 0:
                n = UWORD(1) << (nu - 1);
                break;
            case 1:
                n = (UWORD(1) << (nu - 1)) + 1;
                break;
            case 2:
                n = (UWORD(1) << (nu - 1)) | ((UWORD(1) << (nu - 1)) - 1);
                break;
            default:
                n = n_randbits(state, nu);
        }

        nmod_init(&mod, n);
        ninv = n_preinvert_limb(n);

        /* the bound is b_hi B: the high word must be below b_hi */
        b_hi = UWORD(1) << (nu - 2);

        switch (n_randint(state, 6))
        {
            case 0:
                /* random below the bound */
                hi = n_randtest(state) % b_hi;
                lo = n_randtest(state);
                break;
            case 1:
                /* the largest values */
                hi = b_hi - 1;
                lo = UWORD_MAX - n_randint(state, 2 * n);
                break;
            case 2:
            case 3:
            {
                /* q n + d for d in {0, 1, 2, n - 1}, q random or close to the
                   largest quotient */
                ulong q, d, t[2], u[2];
                t[1] = b_hi - 1;
                t[0] = UWORD_MAX;
                mpn_divrem_1(u, 0, t, 2, n);
                q = u[0];  /* largest quotient, which fits in a word */
                if (n_randint(state, 2))
                    q -= n_randint(state, FLINT_MIN(q, 1000) + 1);
                else
                    q = n_randtest(state) % (q + 1);
                switch (n_randint(state, 4))
                {
                    case 0: d = 0; break;
                    case 1: d = 1; break;
                    case 2: d = 2; break;
                    default: d = n - 1;
                }
                umul_ppmm(hi, lo, q, n);
                add_ssaaaa(hi, lo, hi, lo, 0, d);
                if (hi >= b_hi)
                {
                    umul_ppmm(hi, lo, q, n);
                }
                break;
            }
            case 4:
            {
                /* a product of two residues, when it is below the bound */
                ulong a = n_randint(state, 2) ? n - 1 - n_randint(state, 4) : n_randtest(state) % n;
                ulong b = n_randint(state, 2) ? n - 1 - n_randint(state, 4) : n_randtest(state) % n;
                a = (a >= n) ? n - 1 : a;
                b = (b >= n) ? n - 1 : b;
                umul_ppmm(hi, lo, a, b);
                if (hi >= b_hi)
                {
                    hi = n_randtest(state) % b_hi;
                    lo = n_randtest(state);
                }
                break;
            }
            default:
            {
                /* a sum of k products of residues n - 1, k as large as the
                   bound allows */
                slong k;
                ulong p_hi, p_lo;
                umul_ppmm(p_hi, p_lo, n - 1, n - 1);
                hi = 0;
                lo = 0;
                for (k = 0; k < 64; k++)
                {
                    ulong t_hi = hi, t_lo = lo;
                    add_ssaaaa(t_hi, t_lo, t_hi, t_lo, p_hi, p_lo);
                    if (t_hi < hi || t_hi >= b_hi)
                        break;
                    hi = t_hi;
                    lo = t_lo;
                }
                break;
            }
        }

        if (hi >= b_hi)
            TEST_FUNCTION_FAIL("value above the bound\nn = %wu\nhi = %wu\nlo = %wu\n", n, hi, lo);

        s = n_ll_mod_preinv(hi, lo, n, ninv);
        NMOD_RED2_BOUNDED(r, hi, lo, mod);

        if (r != s)
            TEST_FUNCTION_FAIL("n = %wu\nhi = %wu\nlo = %wu\nr = %wu\ns = %wu\n", n, hi, lo, r, s);

        /* the result may overwrite an argument */
        NMOD_RED2_BOUNDED(lo, hi, lo, mod);
        if (lo != s)
            TEST_FUNCTION_FAIL("aliasing\nn = %wu\nr = %wu\ns = %wu\n", n, lo, s);
    }

    /* NMOD_RED3_BOUNDED against n_lll_mod_preinv, on three-word values
       below 2^(nu + 2 FLINT_BITS - 2), i.e. with a top word below
       2^(nu - 2): random values, the largest values, sums of products */
    for (i = 0; i < 100000 * flint_test_multiplier(); i++)
    {
        nmod_t mod;
        ulong n, ninv, hi, me, lo, r, s, b_hi;
        flint_bitcnt_t nu = 3 + n_randint(state, FLINT_BITS - 3);

        switch (n_randint(state, 4))
        {
            case 0:
                n = UWORD(1) << (nu - 1);
                break;
            case 1:
                n = (UWORD(1) << (nu - 1)) | ((UWORD(1) << (nu - 1)) - 1);
                break;
            default:
                n = n_randbits(state, nu);
        }

        nmod_init(&mod, n);
        ninv = n_preinvert_limb(n);
        b_hi = UWORD(1) << (nu - 2);

        switch (n_randint(state, 3))
        {
            case 0:
                hi = n_randtest(state) % b_hi;
                me = n_randtest(state);
                lo = n_randtest(state);
                break;
            case 1:
                hi = b_hi - 1;
                me = UWORD_MAX - n_randint(state, 2);
                lo = UWORD_MAX - n_randint(state, 2 * n);
                break;
            default:
            {
                /* a sum of products of residues n - 1 */
                slong k, kmax = 1 + n_randint(state, 1000);
                ulong p_hi, p_lo;
                umul_ppmm(p_hi, p_lo, n - 1, n - 1);
                hi = me = lo = 0;
                for (k = 0; k < kmax; k++)
                {
                    ulong t2 = hi, t1 = me, t0 = lo;
                    add_sssaaaaaa(t2, t1, t0, t2, t1, t0, UWORD(0), p_hi, p_lo);
                    if (t2 >= b_hi)
                        break;
                    hi = t2;
                    me = t1;
                    lo = t0;
                }
                break;
            }
        }

        s = n_lll_mod_preinv(hi % n, me, lo, n, ninv);
        NMOD_RED3_BOUNDED(r, hi, me, lo, mod);

        if (r != s)
            TEST_FUNCTION_FAIL("NMOD_RED3_BOUNDED\nn = %wu\nhi = %wu\nme = %wu\nlo = %wu\nr = %wu\ns = %wu\n",
                               n, hi, me, lo, r, s);
    }

    TEST_FUNCTION_END(state);
}
