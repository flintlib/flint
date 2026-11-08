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

/*
    The reciprocal and the quotient with each part rounded in the mode
    rnd: 1 / (a + bi) = (a - bi) / s and (a + bi) / (c + di) = (u + vi) / s
    with u = ac + bd, v = bc - ad and s = c^2 + d^2 (or a^2 + b^2),
    the numerators and the denominator formed exactly (so that each
    part is a correctly rounded quotient of two exact numbers) when
    the terms of each sum lie within ACF_DIV_EXACT_GAP bits of each
    other, otherwise (terms too far apart to interact, so without
    cancellation) rounded to a working precision and the quotient
    certified by _arf_set_mpn_err, with retries.
*/

#define ACF_DIV_EXACT_GAP(prec) (2 * (prec) + 128)

/* the gap between the exponents of two nonzero finite numbers (0 when
   either is zero), clamped */
static slong
_acf_gap(const arf_t x, const arf_t y)
{
    fmpz_t d;
    slong g;

    if (arf_is_zero(x) || arf_is_zero(y))
        return 0;
    fmpz_init(d);
    fmpz_sub(d, ARF_EXPREF(x), ARF_EXPREF(y));
    fmpz_abs(d, d);
    g = fmpz_fits_si(d) ? fmpz_get_si(d) : WORD_MAX / 2;
    fmpz_clear(d);
    return g;
}

/* x y + z w exactly, or at p bits (rounded down) when the terms are
   too far apart to interact (the sum is then the larger term within
   2^-(p - 2) relative); sets *inexact */
static void
_acf_sum2(arf_t res, const arf_t x, const arf_t y, const arf_t z, const arf_t w, int subtract, slong gap, slong p, int * inexact)
{
    arf_t t, u;

    arf_init(t);
    arf_init(u);
    arf_mul(t, x, y, ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_mul(u, z, w, ARF_PREC_EXACT, ARF_RND_DOWN);
    if (_acf_gap(t, u) <= gap)
    {
        if (subtract)
            arf_sub(res, t, u, ARF_PREC_EXACT, ARF_RND_DOWN);
        else
            arf_add(res, t, u, ARF_PREC_EXACT, ARF_RND_DOWN);
    }
    else
    {
        if (subtract)
            arf_sub(res, t, u, p, ARF_RND_DOWN);
        else
            arf_add(res, t, u, p, ARF_RND_DOWN);
        *inexact = 1;
    }
    arf_clear(t);
    arf_clear(u);
}

/* (u + vi) / s with u, v, s exact or within 2^-(p - 2) relative, into
   res (whose parts may alias nothing of the inputs), with the sign of
   the imaginary part flipped (negimag); returns the inexact flags */
static int
_acf_quotient(acf_t res, const arf_t u, const arf_t v, const arf_t s, int inexact, slong p, slong prec, arf_rnd_t rnd)
{
    int s1, s2;

    if (!inexact)
    {
        s1 = arf_div(acf_realref(res), u, s, prec, rnd);
        s2 = arf_div(acf_imagref(res), v, s, prec, rnd);
        return s1 | (s2 << 1);
    }
    else
    {
        /* the quotients within 1 + 2 + 2 ulps at p bits (the division,
           and the relative errors of the operands): K = 6 */
        /* (at p bits, the precision of the sums: their relative errors
           2^(1 - p) are 2 ulps at p bits, not at the limb boundary) */
        arf_t q1, q2;
        arf_init(q1);
        arf_init(q2);
        arf_div(q1, u, s, p, ARF_RND_DOWN);
        arf_div(q2, v, s, p, ARF_RND_DOWN);
        s1 = _acf_elem_set(acf_realref(res), q1, _acf_part_gap(q1, q2), p, 6, WORD_MIN, 0, prec, rnd);
        s2 = (s1 == ARF_ELEM_UNCERTAIN) ? s1 : _acf_elem_set(acf_imagref(res), q2, _acf_part_gap(q2, q1), p, 6, WORD_MIN, 0, prec, rnd);
        arf_clear(q1);
        arf_clear(q2);
        if (s1 == ARF_ELEM_UNCERTAIN || s2 == ARF_ELEM_UNCERTAIN)
            return ARF_ELEM_UNCERTAIN;
        return s1 | (s2 << 1);
    }
}

#define ACF_PART_NONFINITE(x) (arf_is_special(x) && !arf_is_zero(x))

int
acf_inv(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
{
    arf_srcptr a = acf_realref(x), b = acf_imagref(x);
    arf_t s, nb;
    acf_t r;
    slong p;
    int status, inexact;

    if (arf_is_zero(b))
    {
        status = arf_ui_div(acf_realref(res), 1, a, prec, rnd);
        arf_zero(acf_imagref(res));
        return status;
    }
    if (arf_is_zero(a))
    {
        /* -1/b: the rounding of the negated value */
        status = arf_ui_div(acf_imagref(res), 1, b, prec, _acf_rnd_flip(rnd));
        arf_neg(acf_imagref(res), acf_imagref(res));
        arf_zero(acf_realref(res));
        return status << 1;
    }
    if (ACF_PART_NONFINITE(a) || ACF_PART_NONFINITE(b))
    {
        arf_nan(acf_realref(res));
        arf_nan(acf_imagref(res));
        return 0;
    }

    arf_init(s);
    arf_init_neg_shallow(nb, b);
    acf_init(r);

    for (p = prec + 64; ; p *= 2)
    {
        inexact = 0;
        _acf_sum2(s, a, a, b, b, 0, ACF_DIV_EXACT_GAP(prec), p, &inexact);
        status = _acf_quotient(r, a, nb, s, inexact, p, prec, rnd);
        if (status != ARF_ELEM_UNCERTAIN)
            break;
        if (p > 4 * prec + ACF_ELEM_RETRY_LIMIT_BITS)
        {
            /* (the approximations rounded as they stand) */
            status = _acf_quotient(r, a, nb, s, 0, p, prec, (rnd >= ARF_RND_FAST) ? arf_rnd_relaxed_to_strict(rnd) : rnd);
            break;
        }
    }
    acf_swap(res, r);
    arf_clear(s);
    acf_clear(r);
    return status;
}

int
acf_div(acf_t res, const acf_t x, const acf_t y, slong prec, arf_rnd_t rnd)
{
    arf_srcptr a = acf_realref(x), b = acf_imagref(x), c = acf_realref(y), d = acf_imagref(y);
    arf_t s, u, v;
    acf_t r;
    slong p;
    int status, inexact;

    if (ACF_PART_NONFINITE(a) || ACF_PART_NONFINITE(b) || ACF_PART_NONFINITE(c) || ACF_PART_NONFINITE(d))
    {
        arf_nan(acf_realref(res));
        arf_nan(acf_imagref(res));
        return 0;
    }

    /* real and imaginary divisors: the parts divided (x / (di) = (b - ai) / d) */
    if (arf_is_zero(d))
    {
        arf_t cc;
        arf_init(cc);
        arf_set(cc, c);
        status = arf_div(acf_realref(res), a, cc, prec, rnd);
        status |= arf_div(acf_imagref(res), b, cc, prec, rnd) << 1;
        arf_clear(cc);
        return status;
    }
    if (arf_is_zero(c))
    {
        arf_t dd, na;
        int s1, s2;
        arf_init(dd);
        arf_init(na);
        arf_set(dd, d);
        arf_neg(na, a);
        s1 = arf_div(acf_realref(res), b, dd, prec, rnd);
        s2 = arf_div(acf_imagref(res), na, dd, prec, rnd);
        arf_clear(dd);
        arf_clear(na);
        return s1 | (s2 << 1);
    }

    arf_init(s);
    arf_init(u);
    arf_init(v);
    acf_init(r);

    for (p = prec + 64; ; p *= 2)
    {
        inexact = 0;
        _acf_sum2(s, c, c, d, d, 0, ACF_DIV_EXACT_GAP(prec), p, &inexact);
        _acf_sum2(u, a, c, b, d, 0, ACF_DIV_EXACT_GAP(prec), p, &inexact);
        _acf_sum2(v, b, c, a, d, 1, ACF_DIV_EXACT_GAP(prec), p, &inexact);
        status = _acf_quotient(r, u, v, s, inexact, p, prec, rnd);
        if (status != ARF_ELEM_UNCERTAIN)
            break;
        if (p > 4 * prec + ACF_ELEM_RETRY_LIMIT_BITS)
        {
            status = _acf_quotient(r, u, v, s, 0, p, prec, (rnd >= ARF_RND_FAST) ? arf_rnd_relaxed_to_strict(rnd) : rnd);
            break;
        }
    }

    acf_swap(res, r);
    arf_clear(s);
    arf_clear(u);
    arf_clear(v);
    acf_clear(r);
    return status;
}
