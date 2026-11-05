/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"
#include "acb_mzeta.h"
#include "acb_mzeta/impl.h"

/* the dual composition: reverse the word and swap the letters */
static slong
dual(ulong * t, const ulong * s, slong r)
{
    unsigned char w[512], v[512];
    slong k, i, rr = 0, z = 0;
    k = _acb_mzeta_mzv_word(w, s, r);
    for (i = 0; i < k; i++)
        v[i] = 1 - w[k - 1 - i];
    for (i = 0; i < k; i++)
    {
        if (v[i] == 0)
            z++;
        else
        {
            t[rr++] = z + 1;
            z = 0;
        }
    }
    return rr;
}

TEST_FUNCTION_START(acb_mzeta_mzv, state)
{
    slong iter;

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        ulong s[16], t[64];
        slong r, rt, i, prec1, prec2, wt;
        arb_t a, b, c;

        arb_init(a);
        arb_init(b);
        arb_init(c);

        r = 1 + n_randint(state, 5);
        wt = 0;
        for (i = 0; i < r; i++)
        {
            s[i] = 1 + n_randint(state, n_randint(state, 2) ? 3 : 6);
            wt += s[i];
        }
        if (s[0] < 2)
            s[0] = 2;
        prec1 = 2 + n_randint(state, n_randint(state, 8) ? 300 : 3000);
        prec2 = 2 + n_randint(state, 400);

        switch (n_randint(state, 3))
        {
            case 0:
                /* duality */
                rt = dual(t, s, r);
                acb_mzeta_mzv(a, s, r, prec1);
                acb_mzeta_mzv(b, t, rt, prec2);
                break;

            case 1:
                /* stuffle zeta(x) zeta(y) = zeta(x,y) + zeta(y,x) + zeta(x+y)
                   for depth-1 x, y >= 2 inserted at the front */
                {
                    ulong x = 2 + n_randint(state, 4), y = 2 + n_randint(state, 4);
                    ulong u[2];
                    arb_zeta_ui(a, x, prec1);
                    arb_zeta_ui(c, y, prec1);
                    arb_mul(a, a, c, prec1);
                    u[0] = x; u[1] = y;
                    acb_mzeta_mzv(b, u, 2, prec2);
                    u[0] = y; u[1] = x;
                    acb_mzeta_mzv(c, u, 2, prec2);
                    arb_add(b, b, c, prec2);
                    arb_zeta_ui(c, x + y, prec2);
                    arb_add(b, b, c, prec2);
                }
                break;

            default:
                /* sum formula: sum over depth-2 indices of weight w of
                   zeta(s1, s2) = zeta(w) */
                {
                    ulong u[2];
                    slong w2 = 3 + n_randint(state, 8);
                    arb_zeta_ui(a, w2, prec1);
                    arb_zero(b);
                    for (i = 2; i < w2; i++)
                    {
                        u[0] = i;
                        u[1] = w2 - i;
                        acb_mzeta_mzv(c, u, 2, prec2);
                        arb_add(b, b, c, prec2);
                    }
                }
                break;
        }

        if (!arb_overlaps(a, b))
        {
            flint_printf("FAIL (iter %wd)\n", iter);
            flint_printf("s = "); for (i = 0; i < r; i++) flint_printf("%wu ", s[i]); flint_printf("\n");
            flint_printf("prec1 = %wd, prec2 = %wd\n", prec1, prec2);
            flint_printf("a = "); arb_printd(a, 30); flint_printf("\n");
            flint_printf("b = "); arb_printd(b, 30); flint_printf("\n");
            flint_abort();
        }

        acb_mzeta_mzv(a, s, r, prec1);
        if (arb_rel_accuracy_bits(a) < prec1 - 4)
        {
            flint_printf("FAIL (accuracy)\n");
            flint_printf("s = "); for (i = 0; i < r; i++) flint_printf("%wu ", s[i]); flint_printf("\n");
            flint_printf("prec1 = %wd\n", prec1);
            flint_printf("a = "); arb_printd(a, 30); flint_printf("\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
    }

    /* the two algorithms */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        ulong s[16];
        unsigned char w[256];
        slong r, i, k, prec = 2 + n_randint(state, 1000);
        arb_t a, b;

        arb_init(a);
        arb_init(b);
        r = 1 + n_randint(state, 5);
        for (i = 0; i < r; i++)
            s[i] = 1 + n_randint(state, 4);
        if (s[0] < 2)
            s[0] = 2;
        k = _acb_mzeta_mzv_word(w, s, r);
        _acb_mzeta_mzv_akhilesh(a, w, k, prec);
        _acb_mzeta_mzv_bsplit(b, w, k, prec);
        if (!arb_overlaps(a, b) || arb_rel_accuracy_bits(b) < prec - 20)
        {
            flint_printf("FAIL (bsplit)\n");
            flint_printf("s = "); for (i = 0; i < r; i++) flint_printf("%wu ", s[i]); flint_printf("\n");
            flint_printf("prec = %wd\n", prec);
            arb_printd(a, 30); flint_printf("\n");
            arb_printd(b, 30); flint_printf("\n");
            flint_abort();
        }
        arb_clear(a);
        arb_clear(b);
    }

    /* star values, interpolation, all values, indices */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        ulong s[64], k, m, nw, i;
        slong r, prec = 2 + n_randint(state, 200);
        arb_ptr v, v2;
        arb_t a, b;
        acb_t t, z;
        int flags;

        arb_init(a);
        arb_init(b);
        acb_init(t);
        acb_init(z);

        /* zeta*({2}^n) = 2 (1 - 2^(1-2n)) zeta(2n) */
        r = 1 + n_randint(state, 5);
        for (i = 0; i < (ulong) r; i++)
            s[i] = 2;
        acb_mzeta_mzv_star(a, s, r, prec);
        arb_zeta_ui(b, 2 * r, prec);
        arb_mul_2exp_si(b, b, 1);
        {
            arb_t c;
            arb_init(c);
            arb_one(c);
            arb_mul_2exp_si(c, c, 1 - 2 * r);
            arb_sub_ui(c, c, 1, prec);
            arb_neg(c, c);
            arb_mul(b, b, c, prec);
            arb_clear(c);
        }
        if (!arb_overlaps(a, b))
        {
            flint_printf("FAIL (star of 2s), r = %wd\n", r);
            flint_abort();
        }

        /* all values of weight <= k against single values, t = 0, 1 */
        k = 2 + n_randint(state, 6);
        flags = n_randint(state, 2) ? ACB_MZETA_MZV_STAR : 0;
        if (n_randint(state, 2))
            flags |= ACB_MZETA_MZV_WEIGHT;
        nw = (flags & ACB_MZETA_MZV_WEIGHT) ? (UWORD(1) << (k - 2)) : (UWORD(1) << (k - 1)) - 1;
        v = _arb_vec_init(nw);
        v2 = _arb_vec_init(nw);
        acb_mzeta_mzv_all(v, k, flags, prec);
        for (i = 0; i < nw; i++)
        {
            m = (flags & ACB_MZETA_MZV_WEIGHT) ? (UWORD(1) << (k - 2)) + i : i + 1;
            r = acb_mzeta_mzv_index_composition(s, m);
            if (acb_mzeta_mzv_index(s, r) != m)
            {
                flint_printf("FAIL (index), m = %wu\n", m);
                flint_abort();
            }
            acb_set_si(t, (flags & ACB_MZETA_MZV_STAR) ? 1 : 0);
            acb_mzeta_mzv_t(z, s, r, t, prec);
            if (flags & ACB_MZETA_MZV_STAR)
                acb_mzeta_mzv_star(a, s, r, prec);
            else
                acb_mzeta_mzv(a, s, r, prec);
            if (!arb_overlaps(v + i, a) || !arb_overlaps(acb_realref(z), a)
                || arb_rel_accuracy_bits(v + i) < prec - 4)
            {
                flint_printf("FAIL (all), k = %wu, flags = %d, m = %wu\n", k, flags, m);
                arb_printd(v + i, 30); flint_printf("\n");
                arb_printd(a, 30); flint_printf("\n");
                acb_printd(z, 30); flint_printf("\n");
                flint_abort();
            }
        }
        _arb_vec_clear(v, nw);
        _arb_vec_clear(v2, nw);

        arb_clear(a);
        arb_clear(b);
        acb_clear(t);
        acb_clear(z);
    }

    /* large first entries (the empty-word values then grow by about
       2^(s_1 - 1) through the double tail recursion): Euler's
       zeta(n, 1) = (n/2) zeta(n+1) - (1/2) sum_{k=1}^{n-2} zeta(n-k) zeta(k+1) */
    for (iter = 0; iter < flint_test_multiplier(); iter++)
    {
        ulong s[2], n, k;
        slong prec = 32 + n_randint(state, 300);
        arb_t a, b, c, d;

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arb_init(d);

        n = 40 + n_randint(state, 60);
        s[0] = n;
        s[1] = 1;
        acb_mzeta_mzv(a, s, 2, prec);

        arb_zeta_ui(b, n + 1, prec + 20);
        arb_mul_ui(b, b, n, prec + 20);
        for (k = 1; k + 2 <= n; k++)
        {
            arb_zeta_ui(c, n - k, prec + 20);
            arb_zeta_ui(d, k + 1, prec + 20);
            arb_submul(b, c, d, prec + 20);
        }
        arb_mul_2exp_si(b, b, -1);

        if (!arb_overlaps(a, b) || arb_rel_accuracy_bits(a) < prec - 4)
        {
            flint_printf("FAIL (zeta(n, 1)), n = %wu, prec = %wd\n", n, prec);
            arb_printd(a, 30); flint_printf("\n");
            arb_printd(b, 30); flint_printf("\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        arb_clear(d);
    }

    TEST_FUNCTION_END(state);
}
