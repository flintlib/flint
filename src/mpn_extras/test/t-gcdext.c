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
#include "thread_support.h"

/* uses _flint_mpn_gcd_test_pair and _flint_mpn_gcd_test_len (t-gcd.c) */

TEST_FUNCTION_START(_flint_mpn_gcdext_inplace, state)
{
    slong iter;
    mpz_t a, b, g, s, t, ref_g, ref_s;

    mpz_init(a);
    mpz_init(b);
    mpz_init(g);
    mpz_init(s);
    mpz_init(t);
    mpz_init(ref_g);
    mpz_init(ref_s);

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong na, nb, an, bn, gn;
        mp_size_t sn;
        mp_ptr ap, bp, gp, sp;

        /* the matrix operations run on several threads for large operands */
        flint_set_num_threads(1 + n_randint(state, 3));

        na = _flint_mpn_gcd_test_len(state);
        nb = n_randint(state, 2) ? na : 1 + n_randint(state, na);
        _flint_mpn_gcd_test_pair(a, b, state, na, nb);

        an = mpz_size(a);
        bn = mpz_size(b);

        /* as documented for GMP's mpn_gcdext: room for an + 1 and
           bn + 1 limbs in the inputs, bn and bn + 1 in the outputs */
        ap = flint_malloc((an + 1) * sizeof(mp_limb_t));
        bp = flint_malloc((bn + 1) * sizeof(mp_limb_t));
        gp = flint_malloc(bn * sizeof(mp_limb_t));
        sp = flint_malloc((bn + 1) * sizeof(mp_limb_t));
        flint_mpn_copyi(ap, a->_mp_d, an);
        flint_mpn_copyi(bp, b->_mp_d, bn);

        gn = _flint_mpn_gcdext_inplace(gp, sp, &sn, ap, an, bp, bn);

        mpz_gcd(ref_g, a, b);

        /* g is the gcd */
        if (gn != (slong) mpz_size(ref_g) || mpn_cmp(gp, ref_g->_mp_d, gn) != 0)
            TEST_FUNCTION_FAIL("gcd: an = %wd, bn = %wd\na = %{mpz}\nb = %{mpz}\n",
                an, bn, a, b);

        /* s a = g mod b */
        mpz_import(s, FLINT_ABS(sn), -1, sizeof(mp_limb_t), 0, 0, sp);
        if (sn < 0)
            mpz_neg(s, s);
        mpz_mul(t, s, a);
        mpz_sub(t, t, ref_g);
        if (!mpz_divisible_p(t, b))
            TEST_FUNCTION_FAIL("cofactor: an = %wd, bn = %wd\na = %{mpz}\nb = %{mpz}\ns = %{mpz}\n",
                an, bn, a, b, s);

        /* GMP documents that its cofactor is normalized so as to be
           unique; it agrees with GMP's */
        mpz_gcdext(g, ref_s, NULL, a, b);
        if (mpz_cmp(s, ref_s) != 0)
            TEST_FUNCTION_FAIL("cofactor differs from GMP: an = %wd, bn = %wd\n"
                "a = %{mpz}\nb = %{mpz}\ns = %{mpz}\nexpected %{mpz}\n",
                an, bn, a, b, s, ref_s);

        flint_free(ap);
        flint_free(bp);
        flint_free(gp);
        flint_free(sp);
    }

    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(g);
    mpz_clear(s);
    mpz_clear(t);
    mpz_clear(ref_g);
    mpz_clear(ref_s);

    flint_set_num_threads(1);

    TEST_FUNCTION_END(state);
}

TEST_FUNCTION_START(flint_mpn_gcdext, state)
{
    slong iter;
    mpz_t a, b, s, ref_g, ref_s;

    mpz_init(a);
    mpz_init(b);
    mpz_init(s);
    mpz_init(ref_g);
    mpz_init(ref_s);

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        slong na, nb, an, bn, gn;
        mp_size_t sn;
        mp_ptr ap, bp, gp, sp;
        int alias;

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

        /* any order, and even operands */
        if (n_randint(state, 2))
            mpz_swap(a, b);
        if (n_randint(state, 4) == 0)
        {
            mpz_mul_2exp(a, a, n_randint(state, 2 * FLINT_BITS));
            mpz_mul_2exp(b, b, n_randint(state, 2 * FLINT_BITS));
        }

        an = mpz_size(a);
        bn = mpz_size(b);

        ap = flint_malloc(an * sizeof(mp_limb_t));
        bp = flint_malloc(bn * sizeof(mp_limb_t));
        flint_mpn_copyi(ap, a->_mp_d, an);
        flint_mpn_copyi(bp, b->_mp_d, bn);

        /* g and s may be written over the inputs (not over each other) */
        alias = n_randint(state, 4);
        gp = (alias == 1) ? ap : (alias == 2) ? bp : flint_malloc(FLINT_MIN(an, bn) * sizeof(mp_limb_t));
        /* s needs room for bn limbs */
        sp = (alias == 3) ? bp : (alias == 2 && an >= bn) ? ap : flint_malloc(bn * sizeof(mp_limb_t));

        gn = flint_mpn_gcdext(gp, sp, &sn, ap, an, bp, bn);

        /* the same g and s as mpz_gcdext */
        mpz_gcdext(ref_g, ref_s, NULL, a, b);
        mpz_import(s, FLINT_ABS(sn), -1, sizeof(mp_limb_t), 0, 0, sp);
        if (sn < 0)
            mpz_neg(s, s);

        if (gn != (slong) mpz_size(ref_g) || mpn_cmp(gp, ref_g->_mp_d, gn) != 0 ||
            mpz_cmp(s, ref_s) != 0)
            TEST_FUNCTION_FAIL("an = %wd, bn = %wd, alias = %d\n"
                "a = %{mpz}\nb = %{mpz}\ns = %{mpz}\nexpected %{mpz}\n",
                an, bn, alias, a, b, s, ref_s);

        if ((gp != ap && sp != ap && mpn_cmp(ap, a->_mp_d, an) != 0) ||
            (gp != bp && sp != bp && mpn_cmp(bp, b->_mp_d, bn) != 0))
            TEST_FUNCTION_FAIL("input modified: an = %wd, bn = %wd\n", an, bn);

        if (gp != ap && gp != bp)
            flint_free(gp);
        if (sp != ap && sp != bp)
            flint_free(sp);
        flint_free(ap);
        flint_free(bp);
    }

    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(s);
    mpz_clear(ref_g);
    mpz_clear(ref_s);

    TEST_FUNCTION_END(state);
}
