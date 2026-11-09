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
    cancellation) rounded to prec + 64 bits and the quotient certified
    by _arf_set_mpn_err.  A part that this leaves undecided lies within
    2^-58 ulps of a rounding boundary beta (the far terms of the sums
    can put it within any distance), and its side of beta is the sign of
    N - beta s for its exact numerator N (s > 0): a sum of exact terms,
    whatever their exponents (acf/side.c).
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

/* the side of the part (imag = 0: real, 1: imaginary) of x / y (of 1 / x
   for y = NULL) against beta: sgn(N - beta s) with s = |y|^2 (|x|^2) and
   N = ac + bd resp. bc - ad (a resp. -b) */
static int
_acf_div_side(const acf_t x, const acf_t y, int imag, const arf_t beta)
{
    arf_srcptr a = acf_realref(x), b = acf_imagref(x);
    arf_srcptr c = (y == NULL) ? a : acf_realref(y), d = (y == NULL) ? b : acf_imagref(y);
    _acf_terms_t T;
    int side;

    _acf_terms_init(T);
    if (y == NULL)
    {
        if (imag)
            _acf_terms_add_monomial(T, -1, b, NULL, NULL, NULL);
        else
            _acf_terms_add_monomial(T, 1, a, NULL, NULL, NULL);
    }
    else if (imag)
    {
        _acf_terms_add_monomial(T, 1, b, c, NULL, NULL);
        _acf_terms_add_monomial(T, -1, a, d, NULL, NULL);
    }
    else
    {
        _acf_terms_add_monomial(T, 1, a, c, NULL, NULL);
        _acf_terms_add_monomial(T, 1, b, d, NULL, NULL);
    }
    _acf_terms_add_monomial(T, -1, beta, c, c, NULL);
    _acf_terms_add_monomial(T, -1, beta, d, d, NULL);
    side = _acf_sign_sum(T);
    _acf_terms_clear(T);
    return side;
}

/* (u + vi) / s into res (which aliases nothing of the inputs) with u, v,
   s exact or (inexact) within 2^-(p - 2) relative: the quotients
   correctly rounded, or within 1 + 2 + 2 ulps at p bits (the division,
   and the relative errors of the operands: K = 6) and certified, an
   undecided part from its side of the boundary (x / y, or 1 / x for y
   = NULL); returns the inexact flags */
static int
_acf_quotient(acf_t res, const arf_t u, const arf_t v, const arf_t s, int inexact, slong p, const acf_t x, const acf_t y, slong prec, arf_rnd_t rnd)
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
        arf_t q[2], beta;
        arf_rnd_t r1 = (rnd >= ARF_RND_FAST) ? arf_rnd_relaxed_to_strict(rnd) : rnd;
        int st[2], i;

        arf_init(q[0]);
        arf_init(q[1]);
        arf_init(beta);
        arf_div(q[0], u, s, p, ARF_RND_DOWN);
        arf_div(q[1], v, s, p, ARF_RND_DOWN);
        for (i = 0; i < 2; i++)
        {
            arf_ptr r = i ? acf_imagref(res) : acf_realref(res);

            st[i] = _acf_elem_set(r, q[i], _acf_part_gap(q[i], q[1 - i]), p, 6, WORD_MIN, 0, prec, rnd);
            if (st[i] == ARF_ELEM_UNCERTAIN)
            {
                /* (an undecided part is far from zero: its numerator
                   is the larger term within 2^-(p - 2)) */
                _acf_round_boundary(beta, q[i], prec, r1);
                st[i] = _acf_round_side(r, beta, _acf_div_side(x, y, i, beta), prec, r1);
            }
        }
        arf_clear(q[0]);
        arf_clear(q[1]);
        arf_clear(beta);
        return st[0] | (st[1] << 1);
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

    p = prec + 64;
    inexact = 0;
    _acf_sum2(s, a, a, b, b, 0, ACF_DIV_EXACT_GAP(prec), p, &inexact);
    status = _acf_quotient(r, a, nb, s, inexact, p, x, NULL, prec, rnd);

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

    p = prec + 64;
    inexact = 0;
    _acf_sum2(s, c, c, d, d, 0, ACF_DIV_EXACT_GAP(prec), p, &inexact);
    _acf_sum2(u, a, c, b, d, 0, ACF_DIV_EXACT_GAP(prec), p, &inexact);
    _acf_sum2(v, b, c, a, d, 1, ACF_DIV_EXACT_GAP(prec), p, &inexact);
    status = _acf_quotient(r, u, v, s, inexact, p, x, y, prec, rnd);

    acf_swap(res, r);
    arf_clear(s);
    arf_clear(u);
    arf_clear(v);
    acf_clear(r);
    return status;
}
