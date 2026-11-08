/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acf.h"
#include "impl.h"
#include "approx_kernel.h"

/*
    The square root and the reciprocal square root with each part
    rounded in the mode rnd: the dyadic values exactly (the roots that
    are dyadic, and their reciprocals when the norm is a power of two),
    otherwise the compositions

        sqrt(a + bi) = (q, b / (2q)) for a >= 0, q = sqrt((r + |a|) / 2),
                       r = |a + bi|, and (|b| / (2q), sgn(b) q) for a < 0
                       (no cancellation in r + |a|),
        1 / sqrt(a + bi) = conj(sqrt(a + bi)) / r

    at m limbs with every operation rounded down (r within 2 ulps, q
    within 3, the quotient within 4; the reciprocal root within 7),
    certified by _arf_set_mpn_err, with retries.
*/

#define ACF_PART_NONFINITE(x) (arf_is_special(x) && !arf_is_zero(x))

/* the parts of sqrt(x) (rsqrt = 0) or 1 / sqrt(x) into Q1 + Q2 i at p
   bits, x with nonzero finite parts */
static void
_acf_sqrt_inner(arf_t Q1, arf_t Q2, const acf_t x, int rsqrt, slong p)
{
    arf_srcptr a = acf_realref(x), b = acf_imagref(x);
    arf_t r, t, q, u;

    arf_init(r);
    arf_init(t);
    arf_init(q);
    arf_init(u);

    arf_sosq(r, a, b, p, ARF_RND_DOWN);
    arf_sqrt(r, r, p, ARF_RND_DOWN);
    arf_abs(t, a);
    arf_add(t, t, r, p, ARF_RND_DOWN);
    arf_mul_2exp_si(t, t, -1);
    arf_sqrt(q, t, p, ARF_RND_DOWN);
    arf_abs(u, b);
    arf_mul_2exp_si(u, u, -1);
    arf_div(u, u, q, p, ARF_RND_DOWN);

    /* (q, u) for a >= 0, (u, q) for a < 0, the imaginary part with the
       sign of b */
    if (arf_sgn(a) >= 0)
    {
        arf_swap(Q1, q);
        arf_swap(Q2, u);
    }
    else
    {
        arf_swap(Q1, u);
        arf_swap(Q2, q);
    }
    if (arf_sgn(b) < 0)
        arf_neg(Q2, Q2);

    if (rsqrt)
    {
        arf_div(Q1, Q1, r, p, ARF_RND_DOWN);
        arf_div(Q2, Q2, r, p, ARF_RND_DOWN);
        arf_neg(Q2, Q2);
    }

    arf_clear(r);
    arf_clear(t);
    arf_clear(q);
    arf_clear(u);
}

static int
_acf_sqrt_rsqrt(acf_t res, const acf_t x, int rsqrt, slong prec, arf_rnd_t rnd)
{
    arf_srcptr a = acf_realref(x), b = acf_imagref(x);
    acf_t w;
    arf_t Q1, Q2;
    slong extra, m;
    int status, s1, s2, exact;

    if (ACF_PART_NONFINITE(a) || ACF_PART_NONFINITE(b))
    {
        arf_nan(acf_realref(res));
        arf_nan(acf_imagref(res));
        return 0;
    }

    /* real and imaginary arguments through the real roots (of |a|, with
       the result rotated; of |b| / 2 resp. 2 |b|, with the parts +-equal) */
    if (arf_is_zero(b))
    {
        arf_t t;
        arf_init(t);
        arf_abs(t, a);
        if (arf_sgn(a) >= 0)
        {
            status = rsqrt ? arf_rsqrt(acf_realref(res), t, prec, rnd) : arf_sqrt(acf_realref(res), t, prec, rnd);
            arf_zero(acf_imagref(res));
        }
        else
        {
            /* sqrt(-|a|) = i sqrt|a|, 1 / sqrt(-|a|) = -i / sqrt|a| */
            arf_rnd_t rr = rsqrt ? _acf_rnd_flip(rnd) : rnd;
            status = (rsqrt ? arf_rsqrt(acf_imagref(res), t, prec, rr) : arf_sqrt(acf_imagref(res), t, prec, rr)) << 1;
            if (rsqrt)
                arf_neg(acf_imagref(res), acf_imagref(res));
            arf_zero(acf_realref(res));
        }
        arf_clear(t);
        return status;
    }
    if (arf_is_zero(a))
    {
        /* sqrt(bi) = sqrt(|b| / 2) (1 + sgn(b) i), 1 / sqrt(bi) =
           (1 - sgn(b) i) / sqrt(2 |b|) */
        arf_t t;
        int neg = (arf_sgn(b) < 0) ^ rsqrt;
        arf_init(t);
        arf_abs(t, b);
        arf_mul_2exp_si(t, t, rsqrt ? 1 : -1);
        if (rsqrt)
        {
            s1 = arf_rsqrt(acf_realref(res), t, prec, rnd);
            s2 = arf_rsqrt(acf_imagref(res), t, prec, neg ? _acf_rnd_flip(rnd) : rnd);
        }
        else
        {
            s1 = arf_sqrt(acf_realref(res), t, prec, rnd);
            s2 = arf_sqrt(acf_imagref(res), t, prec, neg ? _acf_rnd_flip(rnd) : rnd);
        }
        if (neg)
            arf_neg(acf_imagref(res), acf_imagref(res));
        arf_clear(t);
        return s1 | (s2 << 1);
    }

    /* the dyadic roots */
    acf_init(w);
    exact = _acf_sqrt_exact(w, x);
    if (exact && rsqrt)
    {
        /* 1 / w = conj(w) / |w|^2, |w|^2 = |x| a power of two */
        arf_t t;
        arf_init(t);
        arf_sosq(t, acf_realref(w), acf_imagref(w), ARF_PREC_EXACT, ARF_RND_DOWN);
        exact = (ARF_IS_POW2(t));
        if (exact)
        {
            arf_neg(acf_imagref(w), acf_imagref(w));
            arf_div(acf_realref(w), acf_realref(w), t, ARF_PREC_EXACT, ARF_RND_DOWN);
            arf_div(acf_imagref(w), acf_imagref(w), t, ARF_PREC_EXACT, ARF_RND_DOWN);
        }
        arf_clear(t);
    }
    if (exact)
    {
        status = acf_set_round(res, w, prec, (rnd >= ARF_RND_FAST) ? arf_rnd_relaxed_to_strict(rnd) : rnd);
        acf_clear(w);
        return status;
    }
    acf_clear(w);

    arf_init(Q1);
    arf_init(Q2);

    for (extra = 0; ; extra = (extra == 0) ? 32 : 2 * extra)
    {
        ulong K = rsqrt ? 7 : 4;
        arf_rnd_t r1 = rnd;

        m = _arf_elem_limbs(prec, rnd, extra + FLINT_BIT_COUNT(K + 2));

        /* the kernels of arf (arf/mpn_complex.c) for moderate exponents:
           the parts at m limbs within 2^(2 - FLINT_BITS m) of the
           magnitude, which is 2^(E - E(part) + 2) units of the frame of
           the part for the larger exponent E (an error far above a
           part's own size makes _arf_set_mpn_err retry) */
        if (_acf_part_kernel_ok(a) && _acf_part_kernel_ok(b) && m <= (WORD(1) << (FLINT_BITS - 14)))
        {
            nn_srcptr ad, bd;
            slong an, bn, er, ei, E, eb1, eb2;
            ulong stack[2 * 40];
            nn_ptr re, im, y;
            TMP_INIT;

            ARF_GET_MPN_READONLY(ad, an, a);
            ARF_GET_MPN_READONLY(bd, bn, b);
            re = (m <= 40) ? stack : flint_malloc(2 * m * sizeof(ulong));
            im = re + m;
            if (rsqrt)
                _arf_mpn_crsqrt_approx(re, &er, im, &ei, ad, an, ARF_EXP(a), ARF_SGNBIT(a), bd, bn, ARF_EXP(b), FLINT_BITS * m);
            else
                _arf_mpn_csqrt_approx(re, &er, im, &ei, ad, an, ARF_EXP(a), ARF_SGNBIT(a), bd, bn, ARF_EXP(b), FLINT_BITS * m);

            /* (the parts normalised: their exponents) */
            er -= flint_clz(re[m - 1]);
            ei -= flint_clz(im[m - 1]);
            E = FLINT_MAX(er, ei);
            eb1 = E - er + 2;
            eb2 = E - ei + 2;
            /* (the relaxed modes measure the error of a smaller part
               against the larger one: _acf_elem_set) */
            if (rnd >= ARF_RND_FAST)
            {
                eb1 -= (E - er >= 2) ? E - er - 1 : 0;
                eb2 -= (E - ei >= 2) ? E - ei - 1 : 0;
            }

            TMP_START;
            y = TMP_ALLOC((m + 1) * sizeof(ulong));
            flint_mpn_copyi(y, re, m);
            s1 = _arf_set_mpn_err(Q1, y, m, er - FLINT_BITS * m, eb1, eb1, 0, prec, r1);
            if (s1 != ARF_ELEM_UNCERTAIN)
            {
                flint_mpn_copyi(y, im, m);
                s2 = _arf_set_mpn_err(Q2, y, m, ei - FLINT_BITS * m, eb2, eb2, ARF_SGNBIT(b) ^ rsqrt, prec, r1);
            }
            else
                s2 = s1;
            TMP_END;
            if (re != stack)
                flint_free(re);

            if (s1 != ARF_ELEM_UNCERTAIN && s2 != ARF_ELEM_UNCERTAIN)
            {
                arf_swap(acf_realref(res), Q1);
                arf_swap(acf_imagref(res), Q2);
                status = s1 | (s2 << 1);
                break;
            }
            if (FLINT_BITS * m <= 4 * prec + ACF_ELEM_RETRY_LIMIT_BITS)
                continue;
        }

        _acf_sqrt_inner(Q1, Q2, x, rsqrt, FLINT_BITS * m);

        if (FLINT_BITS * m > 4 * prec + ACF_ELEM_RETRY_LIMIT_BITS)
        {
            exact = 1;
            if (r1 >= ARF_RND_FAST)
                r1 = arf_rnd_relaxed_to_strict(r1);
        }

        {
            slong D1 = _acf_part_gap(Q1, Q2), D2 = _acf_part_gap(Q2, Q1);
            s1 = _acf_elem_set(Q1, Q1, D1, FLINT_BITS * m, K, WORD_MIN, exact, prec, r1);
            s2 = (s1 == ARF_ELEM_UNCERTAIN) ? s1 : _acf_elem_set(Q2, Q2, D2, FLINT_BITS * m, K, WORD_MIN, exact, prec, r1);
        }
        if (s1 != ARF_ELEM_UNCERTAIN && s2 != ARF_ELEM_UNCERTAIN)
        {
            arf_swap(acf_realref(res), Q1);
            arf_swap(acf_imagref(res), Q2);
            status = s1 | (s2 << 1);
            break;
        }
    }

    arf_clear(Q1);
    arf_clear(Q2);
    return status;
}

int
acf_sqrt(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
{
    return _acf_sqrt_rsqrt(res, x, 0, prec, rnd);
}

int
acf_rsqrt(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
{
    return _acf_sqrt_rsqrt(res, x, 1, prec, rnd);
}
