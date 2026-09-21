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
#include "fixed.h"

/* random sizes covering balanced, unbalanced and large shapes */
static void
_random_sizes(flint_rand_t state, mp_size_t * an, mp_size_t * bn)
{
    switch (n_randint(state, 9))
    {
        case 7:
            /* register-based division for short divisors */
            *bn = 1 + n_randint(state, FLINT_MPN_DIV_SMALL_BN + 1);
            *an = *bn + n_randint(state, 3 * (*bn) + 3);
            break;
        case 8:
            /* short quotients */
            *bn = 2 + n_randint(state, 200);
            *an = *bn + n_randint(state, FLINT_MIN(*bn, 12));
            break;
        case 6:
            /* above the Newton cutoff of the dispatcher */
            *bn = 1024 + n_randint(state, 300);
            *an = *bn + 1024 + n_randint(state, 1500);
            break;
        case 0:
            *bn = 1 + n_randint(state, 10);
            *an = *bn + n_randint(state, 20);
            break;
        case 1:
            *bn = 1 + n_randint(state, 40);
            *an = *bn + n_randint(state, 400);
            break;
        case 2:
            *bn = 1 + n_randint(state, 500);
            *an = *bn + n_randint(state, 1500);
            break;
        case 3:
            *bn = 1 + n_randint(state, 500);
            *an = *bn + n_randint(state, 500);
            break;
        case 4:
            *bn = 300 + n_randint(state, 200);
            *an = *bn + 1 + n_randint(state, 4);
            break;
        default:
            *bn = 3 + n_randint(state, 60);
            *an = *bn + n_randint(state, 1000);
            break;
    }
}

TEST_FUNCTION_START(flint_mpn_tdiv_qr, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        mp_ptr a, b, q, r, q2, r2, binv;
        mp_size_t an, bn, n;
        int alg;

        _random_sizes(state, &an, &bn);
        /* one-limb quotients (GMP's mpn_div_q does not support them) */
        if (n_randint(state, 10) == 0)
            an = bn;
        n = an - bn + 1;

        a = flint_malloc(an * sizeof(mp_limb_t));
        b = flint_malloc(bn * sizeof(mp_limb_t));
        q = flint_malloc(n * sizeof(mp_limb_t));
        r = flint_malloc(bn * sizeof(mp_limb_t));
        q2 = flint_malloc(n * sizeof(mp_limb_t));
        r2 = flint_malloc(bn * sizeof(mp_limb_t));

        flint_mpn_rrandom(a, state, an);
        flint_mpn_rrandom(b, state, bn);
        if (n_randint(state, 4) == 0)
            b[bn - 1] |= UWORD(1) << (FLINT_BITS - 1);
        if (b[bn - 1] == 0)
            b[bn - 1] = 1;

        /* top limbs of a equal to those of b (or one less), which triggers
           the special cases of the 3/2 division steps */
        if (n_randint(state, 5) == 0 && an > bn)
        {
            mp_size_t i, off = an - bn - n_randint(state, FLINT_MIN(an - bn, 3) + 1);
            for (i = 0; i < bn; i++)
                a[off + i] = b[i];
            if (n_randint(state, 2))
                a[an - 1] = b[bn - 1] - (b[bn - 1] > 1);
        }

        /* sometimes make the quotient nearly exact or a divides b */
        if (n_randint(state, 4) == 0)
        {
            mp_ptr t = flint_malloc((n + bn) * sizeof(mp_limb_t));
            if (n >= bn)
                flint_mpn_mul(t, a, n, b, bn);
            else
                flint_mpn_mul(t, b, bn, a, n);
            flint_mpn_copyi(a, t, an);
            if (n_randint(state, 2))
                mpn_sub_1(a, a, an, n_randint(state, 3));
            else
                mpn_add_1(a, a, an, n_randint(state, 3));
            if (a[an - 1] == 0 && an > bn)
                a[an - 1] = 1;
            flint_free(t);
        }

        mpn_tdiv_qr(q2, r2, 0, a, an, b, bn);

        alg = n_randint(state, 13);
        if (alg == 0)
            flint_mpn_tdiv_qr(q, r, a, an, b, bn);
        else if (alg == 1)
            _flint_mpn_tdiv_qr_newton(q, r, a, an, b, bn);
        else if (alg == 2 && an > 2 * bn && bn >= 3)
            _flint_mpn_tdiv_qr_unbalanced(q, r, a, an, b, bn);
        else if (alg == 3 && bn >= 3 && an >= n + 2)
        {
            binv = flint_malloc((n + 4) * sizeof(mp_limb_t));
            fixed_inv_newton(binv, b, bn, n + 2);
            _flint_mpn_tdiv_qr_preinv(q, r, a, an, b, bn, binv, n + 2);
            flint_free(binv);
        }
        else if (alg == 4)
        {
            flint_mpn_tdiv_q(q, a, an, b, bn);
            flint_mpn_tdiv_r(r, a, an, b, bn);
        }
        else if (alg == 6)
        {
            _flint_mpn_tdiv_qr_preinvn(q, r, a, an, b, bn);
        }
        else if (alg == 7 && bn >= 2)
        {
            _flint_mpn_tdiv_qr_divconquer(q, r, a, an, b, bn);
        }
        else if (alg == 8 && bn >= 2)
        {
            _flint_mpn_tdiv_qr_divconquer(q, NULL, a, an, b, bn);
            flint_mpn_copyi(r, r2, bn);
        }
        else if (alg == 9 && bn >= 2)
        {
            _flint_mpn_tdiv_q_divconquer(q, a, an, b, bn);
            flint_mpn_copyi(r, r2, bn);
        }
        else if (alg == 12)
        {
            /* the remainder may overwrite the dividend */
            mp_ptr ra = flint_malloc(an * sizeof(mp_limb_t));
            flint_mpn_copyi(ra, a, an);
            flint_mpn_tdiv_qr(q, ra, ra, an, b, bn);
            flint_mpn_copyi(r, ra, bn);
            flint_free(ra);
        }
        else if (alg == 10 && bn <= FLINT_MPN_DIV_SMALL_BN)
        {
            _flint_mpn_tdiv_qr_small(q, r, a, an, b, bn);
        }
        else if (alg == 11 && bn <= FLINT_MPN_DIV_SMALL_BN)
        {
            _flint_mpn_tdiv_qr_small(q, NULL, a, an, b, bn);
            flint_mpn_copyi(r, r2, bn);
        }
        else
        {
            alg = 5;
            _flint_mpn_tdiv_qr_newton(q, NULL, a, an, b, bn);
            flint_mpn_copyi(r, r2, bn);
        }

        if (mpn_cmp(q, q2, n) != 0 || mpn_cmp(r, r2, bn) != 0)
            TEST_FUNCTION_FAIL("alg = %d, an = %wd, bn = %wd\na = %{ulong*}\nb = %{ulong*}\nq = %{ulong*}\nq2 = %{ulong*}\nr = %{ulong*}\nr2 = %{ulong*}\n",
                alg, an, bn, a, an, b, bn, q, n, q2, n, r, bn, r2, bn);

        flint_free(a);
        flint_free(b);
        flint_free(q);
        flint_free(r);
        flint_free(q2);
        flint_free(r2);
    }

    /* schoolbook routines with a normalized divisor */
    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        mp_limb_t a[200], w[200], d[64], q[200], q2[200], r2[64], t[300], m[200];
        mp_limb_t dinv, qh, qh2;
        mp_size_t nn, dn, qn, i;
        int alg;

        dn = 3 + n_randint(state, 40);
        nn = dn + n_randint(state, (n_randint(state, 2) ? 8 : 100));
        qn = nn - dn;

        flint_mpn_rrandom(a, state, nn);
        flint_mpn_rrandom(d, state, dn);
        d[dn - 1] |= UWORD(1) << (FLINT_BITS - 1);

        if (n_randint(state, 3) == 0)
        {
            /* near multiples of d */
            mp_size_t mn = FLINT_MAX(qn, 1);
            flint_mpn_rrandom(m, state, mn);
            if (mn >= dn)
                flint_mpn_mul(t, m, mn, d, dn);
            else
                flint_mpn_mul(t, d, dn, m, mn);
            flint_mpn_zero(a, nn);
            flint_mpn_copyi(a, t, FLINT_MIN(nn, mn + dn));
            if (n_randint(state, 2))
                mpn_sub_1(a, a, nn, n_randint(state, 3));
        }
        else if (n_randint(state, 3) == 0 && qn > 0)
        {
            /* top limbs equal to d */
            for (i = 0; i < dn; i++)
                a[nn - dn - 1 + i] = d[i];
        }

        mpn_tdiv_qr(q2, r2, 0, a, nn, d, dn);
        dinv = flint_mpn_preinv1(d[dn - 1], d[dn - 2]);
        flint_mpn_copyi(w, a, nn);

        alg = n_randint(state, 3);
        if (alg == 0)
            qh = _flint_mpn_divrem_basecase_preinv1(q, w, nn, d, dn, dinv);
        else if (alg == 1)
            qh = _flint_mpn_div_basecase_preinv1(q, w, nn, d, dn, dinv);
        else
            qh = _flint_mpn_divapprox_basecase_preinv1(q, w, nn, d, dn, dinv);

        q[qn] = qh;
        qh2 = q2[qn];

        if (alg == 2)
        {
            /* correct or one too large */
            if (mpn_cmp(q, q2, qn + 1) != 0)
            {
                mpn_sub_1(q, q, qn + 1, 1);
                if (mpn_cmp(q, q2, qn + 1) != 0)
                    TEST_FUNCTION_FAIL("divapprox_basecase_preinv1: nn = %wd, dn = %wd\n", nn, dn);
            }
        }
        else if (mpn_cmp(q, q2, qn + 1) != 0 || (alg == 0 && mpn_cmp(w, r2, dn) != 0))
            TEST_FUNCTION_FAIL("sb alg = %d: nn = %wd, dn = %wd, qh = %wu, qh2 = %wu\n", alg, nn, dn, qh, qh2);
    }

    TEST_FUNCTION_END(state);
}
