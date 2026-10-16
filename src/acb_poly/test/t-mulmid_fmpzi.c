/*
    Copyright (C) 2026 Fredrik Johansson

    Developed using Claude Opus 5.5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_poly.h"

/* a midpoint part of about 'bits' bits and magnitude 2^e (or zero) */
static void
_rand_mid(arf_t m, flint_rand_t state, slong bits, slong e)
{
    switch (n_randint(state, 8))
    {
        case 0:
            arf_zero(m);
            break;
        case 1:
            /* extreme mantissa */
            arf_one(m);
            arf_mul_2exp_si(m, m, bits);
            arf_sub_ui(m, m, 1, ARF_PREC_EXACT, ARF_RND_DOWN);
            arf_mul_2exp_si(m, m, e - bits);
            break;
        default:
            arf_randtest(m, state, bits, 4);
            arf_mul_2exp_si(m, m, e);
    }

    if (n_randint(state, 2))
        arf_neg(m, m);
}

/* coefficients 2^(c i) (1 + noise) times a common factor; parts
   possibly zero or much smaller; radii zero, small or crude; 'spread'
   adds exponent variation (large values make it poorly scaled) */
static void
_rand_poly(acb_ptr v, slong len, flint_rand_t state, slong bits, slong c, slong spread)
{
    slong i, e0 = (slong) n_randint(state, 200) - 100;
    int kind = n_randint(state, 6);
    int rkind = n_randint(state, 4);

    for (i = 0; i < len; i++)
    {
        slong e = e0 + c * i + (spread ? (slong) n_randint(state, spread) : 0);
        slong ei = e - ((kind == 3) ? (slong) n_randint(state, 2 * bits + 1) : 0);

        _rand_mid(arb_midref(acb_realref(v + i)), state, bits, e);
        _rand_mid(arb_midref(acb_imagref(v + i)), state, bits, ei);

        if (kind == 1)
            arf_zero(arb_midref(acb_imagref(v + i)));
        if (kind == 2)
            arf_zero(arb_midref(acb_realref(v + i)));

        mag_zero(arb_radref(acb_realref(v + i)));
        mag_zero(arb_radref(acb_imagref(v + i)));

        if (rkind == 1 || (rkind == 2 && n_randint(state, 2)))
        {
            mag_set_ui_2exp_si(arb_radref(acb_realref(v + i)), 1 + n_randint(state, 1000), e - bits);
            if (n_randint(state, 4))
                mag_set_ui_2exp_si(arb_radref(acb_imagref(v + i)), 1 + n_randint(state, 1000), e - bits);
        }
        else if (rkind == 3)
        {
            mag_set_ui_2exp_si(arb_radref(acb_realref(v + i)), 1 + n_randint(state, 1000), e - n_randint(state, 20));
            mag_set_ui_2exp_si(arb_radref(acb_imagref(v + i)), 1 + n_randint(state, 1000), e - n_randint(state, 20));
        }
    }

    /* zero tails */
    if (n_randint(state, 5) == 0)
        _acb_vec_zero(v + len / 2, len - len / 2);
}

/* a point of the ball x: the midpoint or a corner */
static void
_rand_point(acb_t p, const acb_t x, flint_rand_t state)
{
    arf_t r;
    int j;

    arf_init(r);
    acb_set(p, x);

    for (j = 0; j < 2; j++)
    {
        arb_ptr t = j ? acb_imagref(p) : acb_realref(p);

        if (n_randint(state, 3))
        {
            arf_set_mag(r, arb_radref(t));
            if (n_randint(state, 2))
                arf_neg(r, r);
            arf_add(arb_midref(t), arb_midref(t), r, ARF_PREC_EXACT, ARF_RND_DOWN);
        }
        mag_zero(arb_radref(t));
    }

    arf_clear(r);
}

/* the exact window [nlo, nhi) of the product of exact (zero radius)
   polynomials, with exact arf arithmetic */
static void
_exact_mulmid(acb_ptr z, acb_srcptr x, slong xlen, acb_srcptr y, slong ylen,
              slong nlo, slong nhi)
{
    arf_t t;
    slong i, k;

    arf_init(t);

    for (k = nlo; k < nhi; k++)
    {
        arb_ptr zr = acb_realref(z + k - nlo), zi = acb_imagref(z + k - nlo);

        acb_zero(z + k - nlo);

        for (i = FLINT_MAX(0, k - ylen + 1); i <= FLINT_MIN(k, xlen - 1); i++)
        {
            const arf_struct * ar = arb_midref(acb_realref(x + i));
            const arf_struct * ai = arb_midref(acb_imagref(x + i));
            const arf_struct * br = arb_midref(acb_realref(y + k - i));
            const arf_struct * bi = arb_midref(acb_imagref(y + k - i));

            arf_mul(t, ar, br, ARF_PREC_EXACT, ARF_RND_DOWN);
            arf_add(arb_midref(zr), arb_midref(zr), t, ARF_PREC_EXACT, ARF_RND_DOWN);
            arf_mul(t, ai, bi, ARF_PREC_EXACT, ARF_RND_DOWN);
            arf_sub(arb_midref(zr), arb_midref(zr), t, ARF_PREC_EXACT, ARF_RND_DOWN);
            arf_mul(t, ar, bi, ARF_PREC_EXACT, ARF_RND_DOWN);
            arf_add(arb_midref(zi), arb_midref(zi), t, ARF_PREC_EXACT, ARF_RND_DOWN);
            arf_mul(t, ai, br, ARF_PREC_EXACT, ARF_RND_DOWN);
            arf_add(arb_midref(zi), arb_midref(zi), t, ARF_PREC_EXACT, ARF_RND_DOWN);
        }
    }

    arf_clear(t);
}

TEST_FUNCTION_START(acb_poly_mulmid_fmpzi, state)
{
    slong iter, done = 0;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        slong len1, len2, nlo, nhi, w, prec, bits, c, spread, i, len;
        acb_ptr a, b, r, s, pa, pb, ex;
        int squaring, alias, ok;

        len1 = 1 + n_randint(state, n_randint(state, 4) ? 30 : 300);
        len2 = 1 + n_randint(state, n_randint(state, 4) ? 30 : 300);
        prec = 2 + n_randint(state, n_randint(state, 4) ? 300 : 3000);
        bits = 1 + n_randint(state, prec + 20);
        c = (slong) n_randint(state, 41) - 20;
        spread = n_randint(state, 3) ? n_randint(state, 30) : n_randint(state, 4 * prec + 2000);
        squaring = n_randint(state, 5) == 0;
        if (squaring)
            len2 = len1;

        len = len1 + len2 - 1;
        nlo = n_randint(state, len + 1);
        nhi = nlo + 1 + n_randint(state, len + 2 - nlo);
        if (n_randint(state, 3) == 0)
            nlo = 0;
        w = nhi - nlo;

        a = _acb_vec_init(len1);
        b = squaring ? a : _acb_vec_init(len2);
        r = _acb_vec_init(FLINT_MAX(w, len1));
        s = _acb_vec_init(w);
        pa = _acb_vec_init(len1);
        pb = _acb_vec_init(len2);
        ex = _acb_vec_init(w);

        _rand_poly(a, len1, state, bits, c, spread);
        if (!squaring)
            _rand_poly(b, len2, state, 1 + n_randint(state, prec + 20), c, spread);

        alias = !squaring && w <= len1 && n_randint(state, 4) == 0;

        if (alias)
        {
            _acb_vec_set(r, a, len1);
            ok = _acb_poly_mulmid_fmpzi(r, r, len1, b, len2, nlo, nhi, prec);
        }
        else
            ok = _acb_poly_mulmid_fmpzi(r, a, len1, b, len2, nlo, nhi, prec);

        if (ok)
        {
            done++;

            /* containment of the exact product of points of the balls */
            for (i = 0; i < len1; i++)
                _rand_point(pa + i, a + i, state);
            if (squaring)
                _acb_vec_set(pb, pa, len2);
            else
                for (i = 0; i < len2; i++)
                    _rand_point(pb + i, b + i, state);

            if (nlo < len)
                _exact_mulmid(ex, pa, len1, pb, len2, nlo, FLINT_MIN(nhi, len));
            for (i = FLINT_MAX(len - nlo, 0); i < w; i++)
                acb_zero(ex + i);

            for (i = 0; i < w; i++)
            {
                if (!acb_contains(r + i, ex + i))
                {
                    flint_printf("FAIL: containment\n");
                    flint_printf("len1 = %wd, len2 = %wd, nlo = %wd, nhi = %wd, prec = %wd, "
                                 "squaring = %d, alias = %d, i = %wd\n",
                                 len1, len2, nlo, nhi, prec, squaring, alias, i);
                    flint_printf("r = "); acb_printd(r + i, 30); flint_printf("\n");
                    flint_printf("ex = "); acb_printd(ex + i, 30); flint_printf("\n");
                    flint_abort();
                }
            }

            /* about as accurate as the transposed algorithm */
            if (nlo < len)
                _acb_poly_mulmid_transpose(s, a, len1, b, len2, nlo, FLINT_MIN(nhi, len), prec);
            for (i = 0; i < FLINT_MIN(nhi, len) - nlo; i++)
            {
                slong acc1 = acb_rel_accuracy_bits(r + i);
                slong acc2 = acb_rel_accuracy_bits(s + i);

                if (acc2 > -100000 && acc2 < prec && acc1 < acc2 - 3)
                {
                    flint_printf("FAIL: accuracy\n");
                    flint_printf("len1 = %wd, len2 = %wd, nlo = %wd, nhi = %wd, prec = %wd, "
                                 "squaring = %d, i = %wd, acc %wd vs %wd\n",
                                 len1, len2, nlo, nhi, prec, squaring, i, acc1, acc2);
                    flint_printf("r = "); acb_printd(r + i, 30); flint_printf("\n");
                    flint_printf("s = "); acb_printd(s + i, 30); flint_printf("\n");
                    flint_abort();
                }
            }
        }

        _acb_vec_clear(a, len1);
        if (!squaring)
            _acb_vec_clear(b, len2);
        _acb_vec_clear(r, FLINT_MAX(w, len1));
        _acb_vec_clear(s, w);
        _acb_vec_clear(pa, len1);
        _acb_vec_clear(pb, len2);
        _acb_vec_clear(ex, w);
    }

    /* the well-scaled inputs must mostly be accepted */
    if (done < 100 * flint_test_multiplier())
    {
        flint_printf("FAIL: only %wd accepted\n", done);
        flint_abort();
    }

    TEST_FUNCTION_END(state);
}
