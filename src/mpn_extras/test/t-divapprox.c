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

TEST_FUNCTION_START(flint_mpn_divapprox, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        mp_ptr a, b, q, q2, r2, t;
        mp_size_t an, bn, n, i;
        int ok;

        /* the shapes with operands of a thousand limbs and more cost a
           hundred times as much as the others, so they are only used in a
           few iterations */
        if (n_randint(state, 20) == 0)
        {
            if (n_randint(state, 2))
            {
                /* above the Newton cutoff */
                bn = FLINT_MPN_DIVAPPROX_NEWTON_CUTOFF + n_randint(state, 400);
                an = bn + FLINT_MPN_DIVAPPROX_NEWTON_CUTOFF + n_randint(state, 800);
            }
            else
            {
                bn = 1 + n_randint(state, 1200);
                an = bn + n_randint(state, 1200);
            }
        }
        else
        {
            switch (n_randint(state, 4))
            {
                case 0:
                    bn = 1 + n_randint(state, 10);
                    an = bn + n_randint(state, 20);
                    break;
                case 1:
                    bn = 1 + n_randint(state, 150);
                    an = bn + n_randint(state, 400);
                    break;
                case 2:
                    /* short quotient, long divisor */
                    bn = 1 + n_randint(state, 600);
                    an = bn + n_randint(state, 5);
                    break;
                default:
                    /* long quotient, divisor around the divide and conquer cutoffs */
                    bn = 3 + n_randint(state, 3 * FLINT_MPN_DIVAPPR_DC_CUTOFF);
                    an = bn + n_randint(state, 1000);
                    break;
            }
        }
        n = an - bn + 1;

        a = flint_malloc(an * sizeof(mp_limb_t));
        b = flint_malloc(bn * sizeof(mp_limb_t));
        q = flint_malloc((n + 1) * sizeof(mp_limb_t));
        q2 = flint_malloc((n + 1) * sizeof(mp_limb_t));
        r2 = flint_malloc(bn * sizeof(mp_limb_t));
        t = flint_malloc((an + 1) * sizeof(mp_limb_t));

        flint_mpn_rrandom(a, state, an);
        flint_mpn_rrandom(b, state, bn);
        if (n_randint(state, 4) == 0)
            b[bn - 1] |= UWORD(1) << (FLINT_BITS - 1);
        if (b[bn - 1] == 0)
            b[bn - 1] = 1;

        /* quotients close to exact, remainders close to b */
        if (n_randint(state, 3) == 0)
        {
            flint_mpn_rrandom(q2, state, n);
            if (n >= bn)
                flint_mpn_mul(t, q2, n, b, bn);
            else
                flint_mpn_mul(t, b, bn, q2, n);
            flint_mpn_copyi(a, t, an);
            switch (n_randint(state, 3))
            {
                case 0: mpn_sub_1(a, a, an, n_randint(state, 2)); break;
                case 1: mpn_add_1(a, a, an, n_randint(state, 2)); break;
                default:
                    /* + b - 1 */
                    mpn_add(a, a, an, b, bn);
                    mpn_sub_1(a, a, an, 1);
            }
            if (a[an - 1] == 0)
                a[an - 1] = 1;
        }

        /* largest quotient */
        if (n_randint(state, 30) == 0)
        {
            for (i = 0; i < an; i++)
                a[i] = UWORD_MAX;
            flint_mpn_zero(b, bn);
            b[bn - 1] = 1 + n_randint(state, 2);
        }

        mpn_tdiv_qr(q2, r2, 0, a, an, b, bn);
        flint_mpn_divapprox(q, a, an, b, bn);

        /* q = q2 or q = q2 + 1 */
        ok = (mpn_cmp(q, q2, n) == 0);
        if (!ok)
        {
            flint_mpn_copyi(t, q2, n);
            ok = (mpn_add_1(t, t, n, 1) == 0) && (mpn_cmp(q, t, n) == 0);
        }

        if (!ok)
            TEST_FUNCTION_FAIL("an = %wd, bn = %wd\na = %{ulong*}\nb = %{ulong*}\nq = %{ulong*}\nq2 = %{ulong*}\n",
                an, bn, a, an, b, bn, q, n, q2, n);

        /* with f fraction limbs (or the top an - bn + 1 + f < n limbs):
           floor(a B^f / b) or one more */
        {
            mp_size_t f, qn, pn;
            mp_ptr pa, pq, pq2, pr;

            /* the reference division of the padded numerator is what costs
               here, so f is only made large for short operands */
            if (n_randint(state, 2))
                f = n_randint(state, (an <= 100) ? (3 * an + 5) : 10);
            else
                f = -(mp_size_t) n_randint(state, n);
            qn = an - bn + 1 + f;
            pn = an + FLINT_MAX(f, 0);

            pa = flint_malloc(pn * sizeof(mp_limb_t));
            pq = flint_malloc((qn + 1) * sizeof(mp_limb_t));
            pq2 = flint_malloc((pn + 1) * sizeof(mp_limb_t));
            pr = flint_malloc(bn * sizeof(mp_limb_t));

            flint_mpn_zero(pa, pn - an);
            flint_mpn_copyi(pa + pn - an, a, an);
            mpn_tdiv_qr(pq2, pr, 0, pa, pn, b, bn);
            if (f < 0)
            {
                /* floor(a / (b B^-f)) = floor(floor(a / B^-f) / b) */
                mpn_tdiv_qr(pq2, pr, 0, a - f, an + f, b, bn);
            }

            flint_mpn_divapprox_fraction(pq, a, an, b, bn, f);

            ok = (mpn_cmp(pq, pq2, qn) == 0);
            if (!ok)
            {
                ok = (mpn_add_1(pq2, pq2, qn, 1) == 0) && (mpn_cmp(pq, pq2, qn) == 0);
            }

            if (!ok)
                TEST_FUNCTION_FAIL("fraction: an = %wd, bn = %wd, f = %wd\n", an, bn, f);

            flint_free(pa);
            flint_free(pq);
            flint_free(pq2);
            flint_free(pr);
        }

        flint_free(a);
        flint_free(b);
        flint_free(q);
        flint_free(q2);
        flint_free(r2);
        flint_free(t);
    }

    /* for the shapes of each quotient algorithm (short division, divide
       and conquer for balanced and for long quotients, short quotients):
       quotients next to a power of B, where the approximate quotient with
       its guard limb can round up to B^(qn + 1), and remainders next to 0
       and b; flint_mpn_tdiv_q and flint_mpn_divapprox share this code */
    for (iter = 0; iter < 400 * flint_test_multiplier(); iter++)
    {
        mp_ptr a, b, q, q2, r2, t;
        mp_size_t an, bn, n, i;
        int ok;

        switch (n_randint(state, 4))
        {
            case 0:
                n = FLINT_MPN_DIVAPPROX_SHORT_CUTOFF + n_randint(state, 100);
                bn = FLINT_MAX(2, n - 2) + n_randint(state, 100);
                break;
            case 1:
                bn = FLINT_MPN_TDIV_Q_DC_CUTOFF + n_randint(state, 100);
                n = bn + 3 + n_randint(state, 2 * bn);
                break;
            case 2:
                bn = FLINT_MAX(2, FLINT_MPN_DIV_DC_CUTOFF) + n_randint(state, 60);
                n = 3 * bn + n_randint(state, bn);
                break;
            default:
                n = 1 + n_randint(state, FLINT_MPN_DIVAPPROX_SHORT_CUTOFF);
                bn = n + 6 + n_randint(state, 200);
                break;
        }
        an = bn + n - 1;

        a = flint_malloc((an + 1) * sizeof(mp_limb_t));
        b = flint_malloc(bn * sizeof(mp_limb_t));
        q = flint_malloc((n + 1) * sizeof(mp_limb_t));
        q2 = flint_malloc((n + 1) * sizeof(mp_limb_t));
        r2 = flint_malloc(bn * sizeof(mp_limb_t));
        t = flint_malloc((an + 1) * sizeof(mp_limb_t));

        if (n_randint(state, 3) == 0 && bn >= n + 3)
        {
            /* a = B^an - 1 and b = B^(bn - 1) + beta with beta < B^(bn - n - 2),
               so that floor(a B / b) = B^(n + 1) - 1: the quotient with its
               guard limb can come out as B^(n + 1) */
            for (i = 0; i < an; i++)
                a[i] = UWORD_MAX;
            flint_mpn_zero(b, bn);
            flint_mpn_rrandom(b, state, bn - n - 2);
            b[bn - 1] = 1;
        }
        else if (n_randint(state, 2))
        {
            /* a = B^an - 1 - c for small c, b with a short top limb above
               zero, random or all-ones limbs */
            for (i = 0; i < an; i++)
                a[i] = UWORD_MAX;
            a[0] -= n_randint(state, 3);
            switch (n_randint(state, 3))
            {
                case 0: flint_mpn_zero(b, bn); break;
                case 1: flint_mpn_rrandom(b, state, bn); break;
                default: flint_mpn_store(b, bn, UWORD_MAX); break;
            }
            b[bn - 1] = n_randint(state, 4) ? 1 + n_randint(state, 3) : UWORD(1) << (FLINT_BITS - 1);
        }
        else
        {
            /* a = q b + r with r in {0, 1, b - 1, b - 2}, sometimes with the
               low half of q all ones (with r = b - 1, the partial remainders
               of the short division are then next to the divisor) */
            flint_mpn_rrandom(b, state, bn);
            if (n_randint(state, 2))
                b[bn - 1] |= UWORD(1) << (FLINT_BITS - 1);
            if (b[bn - 1] == 0)
                b[bn - 1] = 1;
            flint_mpn_rrandom(q2, state, n);
            if (n_randint(state, 2))
                flint_mpn_store(q2, n / 2 + n_randint(state, 2), UWORD_MAX);
            if (n >= bn)
                flint_mpn_mul(t, q2, n, b, bn);
            else
                flint_mpn_mul(t, b, bn, q2, n);
            an = bn + n;
            switch (n_randint(state, 4))
            {
                case 0: break;
                case 1: mpn_add_1(t, t, an, 1); break;
                default:
                    mpn_add(t, t, an, b, bn);
                    mpn_sub_1(t, t, an, 1 + (n_randint(state, 2)));
            }
            while (an > bn && t[an - 1] == 0)
                an--;
            flint_mpn_copyi(a, t, an);
            if (a[an - 1] == 0)
                a[an - 1] = 1;
            n = an - bn + 1;
        }

        mpn_tdiv_qr(q2, r2, 0, a, an, b, bn);

        flint_mpn_tdiv_q(q, a, an, b, bn);
        if (mpn_cmp(q, q2, n) != 0)
            TEST_FUNCTION_FAIL("tdiv_q: an = %wd, bn = %wd\na = %{ulong*}\nb = %{ulong*}\n",
                an, bn, a, an, b, bn);

        flint_mpn_divapprox(q, a, an, b, bn);
        ok = (mpn_cmp(q, q2, n) == 0);
        if (!ok)
        {
            flint_mpn_copyi(t, q2, n);
            ok = (mpn_add_1(t, t, n, 1) == 0) && (mpn_cmp(q, t, n) == 0);
        }
        if (!ok)
            TEST_FUNCTION_FAIL("divapprox: an = %wd, bn = %wd\na = %{ulong*}\nb = %{ulong*}\n",
                an, bn, a, an, b, bn);

        flint_free(a);
        flint_free(b);
        flint_free(q);
        flint_free(q2);
        flint_free(r2);
        flint_free(t);
    }

    TEST_FUNCTION_END(state);
}
