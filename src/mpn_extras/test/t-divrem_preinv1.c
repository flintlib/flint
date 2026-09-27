/*
    Copyright (C) 2013 William Hart
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"

TEST_FUNCTION_START(flint_mpn_divrem_preinv1, state)
{
    int i, result;
    mpz_t a, a2, b, q, r, q2;
    gmp_randstate_t st;
    mp_limb_t d1, d2, inv;
    slong s1, s2;

    mpz_init(a);
    mpz_init(a2);
    mpz_init(b);
    mpz_init(q);
    mpz_init(r);

    gmp_randinit_default(st);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        /* also sizes using the divide and conquer division */
        slong maxbits = n_randint(state, 4) ? 200 : 12000;

        do {
            mpz_rrandomb(a, st, n_randint(state, maxbits));
            do {
                mpz_rrandomb(b, st, n_randint(state, maxbits));
            } while (mpz_sgn(b) == 0);

            s1 = a->_mp_size;
            s2 = b->_mp_size;
        } while (s1 < s2 || s2 < 2);

        mpz_set(a2, a);

        /* normalise b */
        b->_mp_d[b->_mp_size - 1] |= ((mp_limb_t) 1 << (GMP_LIMB_BITS - 1));

        d1 = b->_mp_d[b->_mp_size - 1];
        d2 = b->_mp_d[b->_mp_size - 2];

        mpz_fdiv_qr(q, r, a, b);

        inv = flint_mpn_preinv1(d1, d2);

        q2->_mp_d = flint_malloc((s1 - s2 + 1)*sizeof(mp_limb_t));

        q2->_mp_d[s1 - s2] = flint_mpn_divrem_preinv1(q2->_mp_d, a2->_mp_d, a2->_mp_size, b->_mp_d, b->_mp_size, inv);

        /* normalise */
        s1 -= (s2 - 1);
        while (s1 && q2->_mp_d[s1 - 1] == 0) s1--;
        q2->_mp_size = s1;
        q2->_mp_alloc = s1;

        while (s2 && a2->_mp_d[s2 - 1] == 0) s2--;
        a2->_mp_size = s2;

        result = (mpz_cmp(q, q2) == 0 && mpz_cmp(a2, r) == 0);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "%{mpz}\n"
                    "%{mpz}\n"
                    "%{mpz}\n"
                    "%{mpz}\n"
                    "%{mpz}\n"
                    "%{mpz}\n",
                    a, b, q, r, q2, a2);

        flint_free(q2->_mp_d);
    }

    /* the exported schoolbook and divide and conquer steps, against
       mpn_tdiv_qr, including shapes that the dispatch of flint_mpn_tdiv_qr
       and flint_mpn_divapprox does not pass them (dn = 2, short quotients,
       divisors below the cutoffs) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        mp_ptr np, np0, dp, qp, q2, r2, tp;
        mp_size_t nn, dn, qn;
        mp_limb_t qh, dinv;
        int which, ok;

        which = n_randint(state, 7);

        if (which <= 2)
        {
            /* schoolbook: dn >= 3 */
            dn = 3 + n_randint(state, 40);
            nn = dn + n_randint(state, 60);
        }
        else if (which <= 4)
        {
            /* dn >= 2, below and above the divide and conquer cutoffs,
               quotients shorter and longer than the divisor */
            dn = 2 + n_randint(state, n_randint(state, 4) ? 30 : 250);
            nn = dn + n_randint(state, n_randint(state, 2) ? dn + 2 : 300);
            if (n_randint(state, 4) == 0)
                nn = 2 * dn - 1;   /* qn = dn - 1 */
        }
        else
        {
            /* balanced 2n / n steps: n >= 6 */
            dn = 6 + n_randint(state, n_randint(state, 4) ? 40 : 200);
            nn = 2 * dn;
        }

        qn = nn - dn;

        np = flint_malloc(nn * sizeof(mp_limb_t));
        np0 = flint_malloc(nn * sizeof(mp_limb_t));
        dp = flint_malloc(dn * sizeof(mp_limb_t));
        qp = flint_malloc((qn + 1) * sizeof(mp_limb_t));
        q2 = flint_malloc((qn + 1) * sizeof(mp_limb_t));
        r2 = flint_malloc(dn * sizeof(mp_limb_t));
        tp = flint_malloc(dn * sizeof(mp_limb_t));

        flint_mpn_rrandom(np, state, nn);
        flint_mpn_rrandom(dp, state, dn);
        dp[dn - 1] |= UWORD(1) << (FLINT_BITS - 1);

        /* remainders 0 and d - 1 (the approximate quotients are then most
           likely to be one too large) */
        if (n_randint(state, 3) == 0)
        {
            mpn_tdiv_qr(q2, r2, 0, np, nn, dp, dn);
            mpn_sub(np, np, nn, r2, dn);
            if (n_randint(state, 2) && mpn_add(np0, np, nn, dp, dn) == 0)
            {
                mpn_sub_1(np0, np0, nn, 1);
                flint_mpn_copyi(np, np0, nn);
            }
        }

        flint_mpn_copyi(np0, np, nn);
        mpn_tdiv_qr(q2, r2, 0, np, nn, dp, dn);
        dinv = flint_mpn_preinv1(dp[dn - 1], dp[dn - 2]);

        switch (which)
        {
            case 0: qh = _flint_mpn_divrem_basecase_preinv1(qp, np, nn, dp, dn, dinv); break;
            case 1: qh = _flint_mpn_div_basecase_preinv1(qp, np, nn, dp, dn, dinv); break;
            case 2: qh = _flint_mpn_divapprox_basecase_preinv1(qp, np, nn, dp, dn, dinv); break;
            case 3: qh = _flint_mpn_divrem_preinv1(qp, np, nn, dp, dn, dinv, tp); break;
            case 4: qh = _flint_mpn_divapprox_preinv1(qp, np, nn, dp, dn, dinv, tp); break;
            case 5: qh = _flint_mpn_divrem_n_divconquer_preinv1(qp, np, dp, dn, dinv, tp); break;
            default: qh = _flint_mpn_divapprox_n_divconquer_preinv1(qp, np, dp, dn, dinv, tp); break;
        }
        qp[qn] = qh;

        ok = (mpn_cmp(qp, q2, qn + 1) == 0);

        /* the approximate quotients may be one too large, the balanced
           divide and conquer step a few units */
        if (!ok && (which == 2 || which == 4))
        {
            mpn_add_1(q2, q2, qn + 1, 1);
            ok = (mpn_cmp(qp, q2, qn + 1) == 0);
        }
        else if (!ok && which == 6 && mpn_cmp(qp, q2, qn + 1) > 0)
        {
            mpn_sub_n(q2, qp, q2, qn + 1);
            ok = flint_mpn_zero_p(q2 + 1, qn) && q2[0] <= 32;
        }

        /* the divisions with remainder leave it in {np, dn} */
        if (ok && (which == 0 || which == 3 || which == 5))
            ok = (mpn_cmp(np, r2, dn) == 0);

        if (!ok)
            TEST_FUNCTION_FAIL("helper %d, nn = %wd, dn = %wd\nn = %{ulong*}\nd = %{ulong*}\n",
                which, nn, dn, np0, nn, dp, dn);

        flint_free(np);
        flint_free(np0);
        flint_free(dp);
        flint_free(qp);
        flint_free(q2);
        flint_free(r2);
        flint_free(tp);
    }

    mpz_clear(a);
    mpz_clear(a2);
    mpz_clear(b);
    mpz_clear(q);
    mpz_clear(r);
    /* don't clear g */
    gmp_randclear(st);

    TEST_FUNCTION_END(state);
}
