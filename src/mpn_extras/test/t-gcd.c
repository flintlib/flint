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
#include "gmpcompat.h"

/* x = random integer of at most n limbs (uniform, or with long runs of
   ones and zeros), possibly zero */
static void
_mpz_randlimbs(mpz_t x, flint_rand_t state, slong n, int runs)
{
    mp_ptr d = mpz_limbs_write(x, n);
    if (runs)
        flint_mpn_rrandomb(d, state, n * FLINT_BITS - n_randint(state, FLINT_BITS));
    else
        flint_mpn_urandomb(d, state, n * FLINT_BITS);
    MPN_NORM(d, n);
    mpz_limbs_finish(x, n);
}

/* Random operand pair with |a| >= |b|, roughly na and nb limbs, with a
   structure chosen to hit the special cases of the gcd algorithms: large
   common factors, consecutive Fibonacci numbers (all quotients 1), huge
   quotients, nearly equal operands, long runs of ones and zeros. Shared
   by the gcd, gcdext and hgcd tests. */
static void
_flint_mpn_gcd_test_pair(mpz_t a, mpz_t b, flint_rand_t state, slong na, slong nb)
{
    mpz_t g, t;
    mpz_init(g);
    mpz_init(t);

    switch (n_randint(state, 8))
    {
        case 0:
            _mpz_randlimbs(a, state, na, 0);
            _mpz_randlimbs(b, state, nb, 0);
            break;
        case 1:
            _mpz_randlimbs(a, state, na, 1);
            _mpz_randlimbs(b, state, nb, 1);
            break;
        case 2:
        {
            /* common factor */
            slong gn = 1 + n_randint(state, FLINT_MIN(na, nb));
            _mpz_randlimbs(g, state, gn, n_randint(state, 2));
            mpz_add_ui(g, g, 1);
            _mpz_randlimbs(a, state, na - gn + 1, n_randint(state, 2));
            _mpz_randlimbs(b, state, nb - gn + 1, n_randint(state, 2));
            mpz_mul(a, a, g);
            mpz_mul(b, b, g);
            break;
        }
        case 3:
        {
            /* Fibonacci numbers times a cofactor */
            ulong k = n_randint(state, nb * FLINT_BITS * 1.44 + 10) + 2;
            mpz_fib2_ui(a, b, k);
            _mpz_randlimbs(g, state, n_randint(state, 3) + 1, 0);
            mpz_add_ui(g, g, 1);
            mpz_mul(a, a, g);
            mpz_mul(b, b, g);
            break;
        }
        case 4:
        {
            /* a huge quotient somewhere */
            _mpz_randlimbs(a, state, na, n_randint(state, 2));
            _mpz_randlimbs(t, state, n_randint(state, na) + 1, n_randint(state, 2));
            mpz_add_ui(t, t, 1);
            mpz_tdiv_q(b, a, t);
            _mpz_randlimbs(g, state, n_randint(state, nb) + 1, 0);
            mpz_tdiv_q_2exp(g, g, n_randint(state, FLINT_BITS));
            mpz_add(b, b, g);
            break;
        }
        case 5:
            /* nearly equal */
            _mpz_randlimbs(a, state, na, n_randint(state, 2));
            _mpz_randlimbs(t, state, n_randint(state, na) + 1, 0);
            mpz_tdiv_q_2exp(t, t, n_randint(state, FLINT_BITS));
            mpz_add(b, a, t);
            break;
        case 6:
            /* powers of two, small partner */
            _mpz_randlimbs(a, state, na, 1);
            mpz_mul_2exp(a, a, n_randint(state, 200));
            _mpz_randlimbs(b, state, n_randint(state, 5) + 1, 1);
            break;
        default:
            /* equal, or a multiple */
            _mpz_randlimbs(a, state, na, n_randint(state, 2));
            if (n_randint(state, 2))
                mpz_set(b, a);
            else
                flint_mpz_mul_ui(b, a, n_randtest(state) | 1);
            break;
    }

    if (mpz_sgn(a) == 0)
        mpz_set_ui(a, 1 + n_randint(state, 1000));
    if (mpz_sgn(b) == 0)
        mpz_set_ui(b, 1 + n_randint(state, 1000));
    if (mpz_cmpabs(a, b) < 0)
        mpz_swap(a, b);

    mpz_clear(g);
    mpz_clear(t);
}

/* a length: mostly small, sometimes in the subquadratic range, and
   occasionally large enough for the transform-space matrix products and
   the approximate half-gcd (n >= 3 FLINT_MPN_HGCD_REDUCE_CUTOFF) */
static slong
_flint_mpn_gcd_test_len(flint_rand_t state)
{
    switch (n_randint(state, 16))
    {
        case 0: case 1: case 2: case 3: case 4:
            return 1 + n_randint(state, 8);
        case 5: case 6: case 7: case 8:
            return 1 + n_randint(state, 60);
        case 9: case 10: case 11:
            return 1 + n_randint(state, 600);
        case 12: case 13: case 14:
            return 1 + n_randint(state, 3000);
        default:
            return (n_randint(state, 8) == 0) ? 1 + n_randint(state, 12000)
                                              : 1 + n_randint(state, 3000);
    }
}

TEST_FUNCTION_START(_flint_mpn_gcd_inplace, state)
{
    slong iter;
    mpz_t a, b, ref;

    mpz_init(a);
    mpz_init(b);
    mpz_init(ref);

    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        slong na, nb, an, bn, gn;
        mp_ptr ap, bp, gp;
        int alias;

        na = _flint_mpn_gcd_test_len(state);
        nb = n_randint(state, 2) ? na : 1 + n_randint(state, na);
        _flint_mpn_gcd_test_pair(a, b, state, na, nb);

        /* at least one operand must be odd */
        if (mpz_even_p(a) && mpz_even_p(b))
            mpz_add_ui(b, b, 1);
        if (mpz_cmp(a, b) < 0)
            mpz_swap(a, b);

        an = mpz_size(a);
        bn = mpz_size(b);

        ap = flint_malloc(an * sizeof(mp_limb_t));
        bp = flint_malloc(an * sizeof(mp_limb_t));
        flint_mpn_copyi(ap, a->_mp_d, an);
        flint_mpn_copyi(bp, b->_mp_d, bn);

        mpz_gcd(ref, a, b);

        /* the output may be written over either input */
        alias = n_randint(state, 3);
        if (alias == 0)
            gp = flint_malloc(bn * sizeof(mp_limb_t));
        else if (alias == 1)
            gp = ap;
        else
            gp = bp;

        gn = _flint_mpn_gcd_inplace(gp, ap, an, bp, bn);

        if (gn != (slong) mpz_size(ref) || mpn_cmp(gp, ref->_mp_d, gn) != 0)
            TEST_FUNCTION_FAIL("an = %wd, bn = %wd, alias = %d\n"
                "a = %{mpz}\nb = %{mpz}\ng = %{ulong*}\nexpected %{mpz}\n",
                an, bn, alias, a, b, gp, gn, ref);

        if (alias == 0)
            flint_free(gp);
        flint_free(ap);
        flint_free(bp);
    }

    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(ref);

    TEST_FUNCTION_END(state);
}

TEST_FUNCTION_START(flint_mpn_gcd, state)
{
    slong iter;
    mpz_t a, b, ref;

    mpz_init(a);
    mpz_init(b);
    mpz_init(ref);

    for (iter = 0; iter < 4000 * flint_test_multiplier(); iter++)
    {
        slong na, nb, an, bn, gn, gmax;
        mp_ptr ap, bp, gp;
        int alias;

        /* mostly short operands, which have their own code paths */
        if (n_randint(state, 2))
        {
            na = 1 + n_randint(state, 4);
            nb = 1 + n_randint(state, 4);
        }
        else
        {
            na = _flint_mpn_gcd_test_len(state);
            nb = 1 + n_randint(state, na);
        }
        _flint_mpn_gcd_test_pair(a, b, state, FLINT_MAX(na, nb), FLINT_MIN(na, nb));

        /* any order, and any powers of two (common or not) */
        if (n_randint(state, 2))
            mpz_swap(a, b);
        if (n_randint(state, 2))
        {
            ulong e = n_randint(state, 3 * FLINT_BITS);
            mpz_mul_2exp(a, a, e + n_randint(state, 2) * n_randint(state, FLINT_BITS));
            mpz_mul_2exp(b, b, e + n_randint(state, 2) * n_randint(state, FLINT_BITS));
        }

        an = mpz_size(a);
        bn = mpz_size(b);
        gmax = FLINT_MIN(an, bn);

        ap = flint_malloc(an * sizeof(mp_limb_t));
        bp = flint_malloc(bn * sizeof(mp_limb_t));
        flint_mpn_copyi(ap, a->_mp_d, an);
        flint_mpn_copyi(bp, b->_mp_d, bn);

        mpz_gcd(ref, a, b);

        /* the output may be written over either input, which are
           otherwise not modified */
        alias = n_randint(state, 3);
        if (alias == 0)
            gp = flint_malloc(gmax * sizeof(mp_limb_t));
        else if (alias == 1)
            gp = ap;
        else
            gp = bp;

        switch (n_randint(state, 8))
        {
            /* the deprecated aliases */
            case 0: gn = flint_mpn_gcd_full(gp, ap, an, bp, bn); break;
            case 1: gn = flint_mpn_gcd_full2(gp, ap, an, bp, bn, NULL); break;
            default: gn = flint_mpn_gcd(gp, ap, an, bp, bn);
        }

        if (gn != (slong) mpz_size(ref) || mpn_cmp(gp, ref->_mp_d, gn) != 0)
            TEST_FUNCTION_FAIL("an = %wd, bn = %wd, alias = %d\n"
                "a = %{mpz}\nb = %{mpz}\ng = %{ulong*}\nexpected %{mpz}\n",
                an, bn, alias, a, b, gp, gn, ref);

        if ((alias != 1 && mpn_cmp(ap, a->_mp_d, an) != 0) ||
            (alias != 2 && mpn_cmp(bp, b->_mp_d, bn) != 0))
            TEST_FUNCTION_FAIL("input modified: an = %wd, bn = %wd, alias = %d\n",
                an, bn, alias);

        if (alias == 0)
            flint_free(gp);
        flint_free(ap);
        flint_free(bp);
    }

    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(ref);

    TEST_FUNCTION_END(state);
}
