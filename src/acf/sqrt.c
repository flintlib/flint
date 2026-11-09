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
    rounded in the mode rnd: the dyadic roots exactly, otherwise the
    complex kernels of arf (arf/mpn_complex.c) for moderate exponents,
    and where they leave a part undecided (or do not apply) the
    compositions

        sqrt(a + bi) = (q, b / (2q)) for a >= 0, q = sqrt((r + |a|)/2),
                       r = |a + bi|, and (|b| / (2q), sgn(b) q) for a < 0
                       (no cancellation in r + |a|),
        1 / sqrt(a + bi) = conj(sqrt(a + bi)) / r

    at prec + 64 bits with every operation rounded down (r within 2
    ulps, q within 3, the quotient within 4; the reciprocal root within
    7).  A part that these leave undecided lies within 2^-58 ulps of a
    rounding boundary beta, and its side of beta is decided exactly
    (_acf_sqrt_side): the parts are q, u = |b| / (2q) and their
    quotients by r, whose comparisons with beta square into signs of
    polynomials in a, b and beta (acf/side.c).
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

/* sgn(w - g) for g > 0, w = q (form 0), u = |b| / (2q) (form 1), q / r
   (form 2) or u / r (form 3), A = |a|, n = a^2 + b^2 = r^2:

     q > g     <=> r > c = 2 g^2 - A:            c <= 0, or n - c^2 > 0
     u > g     <=> d = b^2 - 2 g^2 A > 2 g^2 r:  d > 0 and d^2 - 4 g^4 n > 0
     q / r > g <=> r > c = 2 g^2 n - A:          c <= 0, or n - c^2 > 0
     u / r > g <=> e = b^2 - 2 g^2 A n > 2 g^2 n r: e > 0 and e^2 - 4 g^4 n^3 > 0

   (from q^2 = (r + A) / 2 and u^2 = b^2 / (2 (r + A)); the equalities
   are excluded where they would make the root dyadic) */
static int
_acf_sqrt_side_w(arf_srcptr a, arf_srcptr b, int form, const arf_t g)
{
    _acf_terms_t C, N, X, Y, Z, G;
    arf_t A, g2, g4, a2, b2;
    int side;

    if (arf_sgn(g) <= 0)
        return 1;

    _acf_terms_init(C);
    _acf_terms_init(N);
    _acf_terms_init(X);
    _acf_terms_init(Y);
    _acf_terms_init(Z);
    _acf_terms_init(G);
    arf_init(A);
    arf_init(g2);
    arf_init(g4);
    arf_init(a2);
    arf_init(b2);

    arf_abs(A, a);
    arf_mul(g2, g, g, ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_mul(g4, g2, g2, ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_mul(a2, a, a, ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_mul(b2, b, b, ARF_PREC_EXACT, ARF_RND_DOWN);

    /* N = n */
    _acf_terms_add_monomial(N, 1, a2, NULL, NULL, NULL);
    _acf_terms_add_monomial(N, 1, b2, NULL, NULL, NULL);

    if (form == 0 || form == 2)
    {
        /* C = c; then n - c^2 */
        if (form == 0)
            _acf_terms_add_monomial(C, 2, g2, NULL, NULL, NULL);
        else
        {
            _acf_terms_add_monomial(C, 2, g2, a2, NULL, NULL);
            _acf_terms_add_monomial(C, 2, g2, b2, NULL, NULL);
        }
        _acf_terms_add_monomial(C, -1, A, NULL, NULL, NULL);
        if (_acf_sign_sum(C) <= 0)
            side = 1;
        else
        {
            _acf_terms_mul(X, C, C, -1);
            _acf_terms_append(X, N);
            side = _acf_sign_sum(X);
        }
    }
    else
    {
        /* C = d resp. e; then d^2 - 4 g^4 n resp. e^2 - 4 g^4 n^3 */
        _acf_terms_add_monomial(C, 1, b2, NULL, NULL, NULL);
        if (form == 1)
            _acf_terms_add_monomial(C, -2, g2, A, NULL, NULL);
        else
        {
            _acf_terms_add_monomial(C, -2, g2, A, a2, NULL);
            _acf_terms_add_monomial(C, -2, g2, A, b2, NULL);
        }
        if (_acf_sign_sum(C) <= 0)
            side = -1;
        else
        {
            _acf_terms_mul(X, C, C, 1);
            _acf_terms_add_monomial(G, -4, g4, NULL, NULL, NULL);
            if (form == 1)
                _acf_terms_mul(Z, N, G, 1);
            else
            {
                _acf_terms_mul(Y, N, N, 1);
                _acf_terms_mul(C, Y, N, 1);
                _acf_terms_mul(Z, C, G, 1);
            }
            _acf_terms_append(X, Z);
            side = _acf_sign_sum(X);
        }
    }

    _acf_terms_clear(C);
    _acf_terms_clear(N);
    _acf_terms_clear(X);
    _acf_terms_clear(Y);
    _acf_terms_clear(Z);
    _acf_terms_clear(G);
    arf_clear(A);
    arf_clear(g2);
    arf_clear(g4);
    arf_clear(a2);
    arf_clear(b2);
    return side;
}

/* the side of the real (imag = 0) or imaginary part of sqrt(x) or 1 /
   sqrt(x) against beta: the part is s w with w as above (q for the real
   part when a >= 0 and for the imaginary part when a < 0, u otherwise;
   divided by r for rsqrt) and s the sign of the imaginary part (that of
   b, negated for rsqrt), so that sgn(s w - beta) = s sgn(w - s beta) */
static int
_acf_sqrt_side(const acf_t x, int rsqrt, int imag, const arf_t beta)
{
    arf_srcptr a = acf_realref(x), b = acf_imagref(x);
    int form = (((arf_sgn(a) >= 0) ^ imag) ? 0 : 1) + (rsqrt ? 2 : 0);
    int s = imag ? (((arf_sgn(b) < 0) ^ rsqrt) ? -1 : 1) : 1;
    arf_t g;
    int side;

    arf_init_set_shallow(g, beta);
    if (s < 0)
        arf_neg(g, g);
    side = s * _acf_sqrt_side_w(a, b, form, g);
    return side;
}

static int
_acf_sqrt_rsqrt(acf_t res, const acf_t x, int rsqrt, slong prec, arf_rnd_t rnd)
{
    arf_srcptr a = acf_realref(x), b = acf_imagref(x);
    acf_t w;
    arf_t Q1, Q2;
    slong m;
    ulong K;
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

    /* the dyadic roots (a dyadic reciprocal root would need |sqrt(x)|^2
       = |x| a power of two, that is sqrt(x) = c (1 +- i) and x = +-2 c^2
       i, the imaginary axis above) */
    if (!rsqrt)
    {
        acf_init(w);
        exact = _acf_sqrt_exact(w, x);
        if (exact)
            status = acf_set_round(res, w, prec, (rnd >= ARF_RND_FAST) ? arf_rnd_relaxed_to_strict(rnd) : rnd);
        acf_clear(w);
        if (exact)
            return status;
    }

    arf_init(Q1);
    arf_init(Q2);

    K = rsqrt ? 7 : 4;
    m = _arf_elem_limbs(prec, rnd, FLINT_BIT_COUNT(K + 2));
    s1 = s2 = ARF_ELEM_UNCERTAIN;

    /* the kernels of arf (arf/mpn_complex.c) for moderate exponents: the
       parts at m limbs within 2^(2 - FLINT_BITS m) of the magnitude,
       which is 2^(E - E(part) + 2) units of the frame of the part for the
       larger exponent E (an error far above a part's own size leaves it
       undecided) */
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
        s1 = _arf_set_mpn_err(Q1, y, m, er - FLINT_BITS * m, eb1, eb1, 0, prec, rnd);
        if (s1 != ARF_ELEM_UNCERTAIN)
        {
            flint_mpn_copyi(y, im, m);
            s2 = _arf_set_mpn_err(Q2, y, m, ei - FLINT_BITS * m, eb2, eb2, ARF_SGNBIT(b) ^ rsqrt, prec, rnd);
        }
        TMP_END;
        if (re != stack)
            flint_free(re);
    }

    if (s1 == ARF_ELEM_UNCERTAIN || s2 == ARF_ELEM_UNCERTAIN)
    {
        /* the compositions at p bits, each part within K ulps of itself;
           an undecided part from its side of the boundary */
        slong p = prec + 64, D1, D2;
        arf_rnd_t r1 = (rnd >= ARF_RND_FAST) ? arf_rnd_relaxed_to_strict(rnd) : rnd;
        arf_t beta, P1, P2;

        arf_init(beta);
        arf_init(P1);
        arf_init(P2);
        _acf_sqrt_inner(Q1, Q2, x, rsqrt, p);
        D1 = _acf_part_gap(Q1, Q2);
        D2 = _acf_part_gap(Q2, Q1);
        s1 = _acf_elem_set(P1, Q1, D1, p, K, WORD_MIN, 0, prec, rnd);
        if (s1 == ARF_ELEM_UNCERTAIN)
        {
            _acf_round_boundary(beta, Q1, prec, r1);
            s1 = _acf_round_side(P1, beta, _acf_sqrt_side(x, rsqrt, 0, beta), prec, r1);
        }
        s2 = _acf_elem_set(P2, Q2, D2, p, K, WORD_MIN, 0, prec, rnd);
        if (s2 == ARF_ELEM_UNCERTAIN)
        {
            _acf_round_boundary(beta, Q2, prec, r1);
            s2 = _acf_round_side(P2, beta, _acf_sqrt_side(x, rsqrt, 1, beta), prec, r1);
        }
        /* (x may alias res) */
        arf_swap(acf_realref(res), P1);
        arf_swap(acf_imagref(res), P2);
        arf_clear(beta);
        arf_clear(P1);
        arf_clear(P2);
    }
    else
    {
        arf_swap(acf_realref(res), Q1);
        arf_swap(acf_imagref(res), Q2);
    }

    arf_clear(Q1);
    arf_clear(Q2);
    return s1 | (s2 << 1);
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
