/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf/elem.h"

/*
    Result assembly for the elementary functions of arf (arf/elem.h).

    An approximation Y = {y, len} of the value V with
    (Y - elo) 2^e <= V <= (Y + ehi) 2^e determines the rounding of V to
    prec bits when the bits of Y between the error and the rounding
    point decide it.  With k the number of bits of Y below the rounding
    point, h the half bit (bit k - 1), D the bits [0, k - 1) and
    err < 2^eb the two-sided bound, V's part below the rounding point is
    h 2^(k - 1) + D + (V - Y), and the rounding is decided when that
    lies strictly between h 2^(k - 1) and (h + 1) 2^(k - 1):

    - two-sided: the bits [eb, k - 1) of Y neither all zero nor all one
      (then 2^eb <= D < 2^(k - 1) - 2^eb);
    - V above Y: the bits [eb, k - 1) all zero is fine too (D < 2^eb,
      V - Y < 2^eb, their sum below 2^(k - 1) when eb <= k - 2), as is
      Y exact at prec bits (D = h = 0) with eb <= k - 1;
    - V below Y: the bits [eb, k - 1) all one is fine (the sum stays
      above h 2^(k - 1) when eb <= k - 2), and Y exact at prec bits with
      eb <= k - 1 means V rounds like Y - 1 unit, whose discarded bits
      are all ones;
    - ARF_RND_NEAR, two-sided: the window all equal to the half bit
      (h = 0 and all zero, or h = 1 and all one) with eb <= k - 4: V
      lies within 2^(eb + 1) <= 2^(k - 3) of a number representable at
      prec bits (Y truncated, resp. one unit above it), far from the
      ties at a quarter ulp or more (also across a power of two, where
      the ulp below is half as large), and rounds to it from either
      side -- a value exponentially close to a representable number
      (tanh x = 1 - 2e^(-2x) + ..., say) needs no more precision than
      any other.

    In each decided case V rounds like Y with its bits below the half
    bit replaced by a single sticky bit (resp. Y - 1 in the last case),
    which _arf_set_round_mpn then does exactly.
*/

/* y with the bits [0, k - 1) cleared and bit 0 set */
static void
_sticky(nn_ptr y, slong k)
{
    slong q = (k - 1) / FLINT_BITS, r = (k - 1) % FLINT_BITS;

    _flint_mpn_zero_small(y, q);
    if (r != 0)
        y[q] &= ~UWORD(0) << r;
    y[0] |= 1;
}

int
_arf_set_mpn_err(arf_t res, nn_ptr y, slong len, slong e, slong eblo, slong ebhi, int sgnbit, slong prec, arf_rnd_t rnd)
{
    slong fix, k, w, eb;
    int decided;

    while (len > 0 && y[len - 1] == 0)
        len--;

    /* a zero approximation decides nothing (the value is nonzero) */
    if (len == 0)
        return ARF_ELEM_UNCERTAIN;

    if (!_arf_rnd_is_strict(rnd))
    {
        /* the error bounds below an eighth resp. 2^-7 of an ulp */
        eb = FLINT_MAX(eblo, ebhi);
        if (eb >= 0)
        {
            k = len * FLINT_BITS - flint_clz(y[len - 1]) - prec;
            if (eb + ((rnd == ARF_RND_FAST) ? 3 : 7) > k)
                return ARF_ELEM_UNCERTAIN;
        }

        if (rnd == ARF_RND_FAST)
        {
            /* the upper bound truncated: within an ulp below the value
               and a sliver above it */
            if (ebhi >= 0)
            {
                slong q = ebhi / FLINT_BITS;

                FLINT_ASSERT(q < len);
                y[len] = mpn_add_1(y + q, y + q, len - q, UWORD(1) << (ebhi % FLINT_BITS));
                len += (y[len] != 0);
            }
            _arf_set_round_mpn(res, &fix, y, len, sgnbit, prec, ARF_RND_DOWN);
        }
        else
            _arf_set_round_mpn(res, &fix, y, len, sgnbit, prec, ARF_RND_NEAR);

        fmpz_set_si(ARF_EXPREF(res), e + len * FLINT_BITS + fix);
        return 1;
    }

    eb = FLINT_MAX(eblo, ebhi);

    if (eb < 0)
    {
        /* an exact value: rounded exactly */
        int inexact = _arf_set_round_mpn(res, &fix, y, len, sgnbit, prec, rnd);
        fmpz_set_si(ARF_EXPREF(res), e + len * FLINT_BITS + fix);
        return inexact;
    }

    k = len * FLINT_BITS - flint_clz(y[len - 1]) - prec;

    /* (at least two bits below the rounding point: the half bit and one
       for the sticky bit; the callers have many guard bits) */
    if (k <= 1)
        return ARF_ELEM_UNCERTAIN;

    /* the bits D = [0, k - 1) below the half bit and the window
       [eb, k - 1) above the error (the whole window: a long run of equal
       bits below the rounding point must not stall the retries) */
    w = k - 1 - eb;
    decided = 0;

    if (w >= 1 && !_arf_mpn_bits_uniform(y, eb, k - 1, 0) && !_arf_mpn_bits_uniform(y, eb, k - 1, 1))
        decided = 1;
    else if (rnd == ARF_RND_NEAR && eb <= k - 4 &&
        _arf_mpn_bits_uniform(y, eb, k - 1, (y[(k - 1) / FLINT_BITS] >> ((k - 1) % FLINT_BITS)) & 1))
        decided = 1;
    else if (k == 1 || _arf_mpn_bits_uniform(y, 0, k - 1, 0))
    {
        /* D = 0 */
        if (eb > k - 1)
            return ARF_ELEM_UNCERTAIN;

        if (eblo < 0)
            decided = 1;
        else if (ebhi < 0 && (rnd != ARF_RND_NEAR || eb <= k - 2))
        {
            /* V just below Y: rounds like Y - 1 unit (all ones below; for
               rounding to nearest, the error below the half bit of the
               shorter Y - 1 when Y is a power of two) */
            mpn_sub_1(y, y, len, 1);
            while (y[len - 1] == 0)
                len--;
            _arf_set_round_mpn(res, &fix, y, len, sgnbit, prec, rnd);
            fmpz_set_si(ARF_EXPREF(res), e + len * FLINT_BITS + fix);
            return 1;
        }
    }
    else if (_arf_mpn_bits_uniform(y, 0, k - 1, 1))
    {
        /* D = 2^(k - 1) - 1 */
        if (eblo < 0 && eb <= 0)
            decided = 1;
        else if (ebhi < 0 && eb <= k - 2)
            decided = 1;
    }
    else if (w >= 1)
    {
        /* the window alone uniform */
        if (_arf_mpn_bits_uniform(y, eb, k - 1, 0))
            decided = (eblo < 0);
        else
            decided = (ebhi < 0);
    }

    if (!decided)
        return ARF_ELEM_UNCERTAIN;

    _sticky(y, k);
    _arf_set_round_mpn(res, &fix, y, len, sgnbit, prec, rnd);
    fmpz_set_si(ARF_EXPREF(res), e + len * FLINT_BITS + fix);
    return 1;
}

int
_arf_one_plus_tiny(arf_t res, int sgnbit, int dsgn, slong prec, arf_rnd_t rnd)
{
    arf_rnd_t r = _arf_rnd_resolve(rnd, sgnbit);

    if (r == ARF_RND_UP && !dsgn)
    {
        /* 1 + 2^(1 - prec) */
        arf_set_ui_2exp_si(res, 1, 1 - prec);
        arf_add_ui(res, res, 1, prec, ARF_RND_DOWN);
    }
    else if (r == ARF_RND_DOWN && dsgn)
    {
        /* 1 - 2^-prec */
        arf_set_ui_2exp_si(res, 1, -prec);
        arf_sub_ui(res, res, 1, prec, ARF_RND_DOWN);
        arf_neg(res, res);
    }
    else
        arf_one(res);

    if (sgnbit)
        arf_neg(res, res);
    return 1;
}

int
_arf_set_tiny_rel_err(arf_t res, const arf_t x, int side, slong P, slong prec, arf_rnd_t rnd)
{
    nn_srcptr xd;
    slong xn, len, eb;
    nn_ptr t;
    int status;
    fmpz_t e;
    TMP_INIT;

    FLINT_ASSERT(P >= prec + 3);

    ARF_GET_MPN_READONLY(xd, xn, x);

    /* |x delta| < |x| 2^-P < 2^(E - P) = 2^eb units of the bottom limb of
       Y = x padded with zero limbs to at least P bits (a looser bound
       is as good once the window has a limb) */
    P = FLINT_MIN(P, prec + FLINT_BITS * (xn + 2));
    len = FLINT_MAX(xn, (P + FLINT_BITS - 1) / FLINT_BITS);
    eb = len * FLINT_BITS - P;

    TMP_START;
    t = TMP_ALLOC(sizeof(ulong) * (len + 1));
    _flint_mpn_zero_small(t, len - xn);
    flint_mpn_copyi(t + len - xn, xd, xn);

    /* the exponent of the unit, E - FLINT_BITS len (res may alias x) */
    fmpz_init(e);
    fmpz_sub_ui(e, ARF_EXPREF(x), FLINT_BITS * len);

    if (!COEFF_IS_MPZ(*e))
        status = _arf_set_mpn_err(res, t, len, *e, (side < 0) ? eb : -1,
            (side > 0) ? eb : -1, ARF_SGNBIT(x), prec, rnd);
    else
    {
        /* a bignum exponent: round with a zero unit exponent, then shift */
        status = _arf_set_mpn_err(res, t, len, 0, (side < 0) ? eb : -1,
            (side > 0) ? eb : -1, ARF_SGNBIT(x), prec, rnd);
        if (status == 1)
            fmpz_add(ARF_EXPREF(res), ARF_EXPREF(res), e);
    }

    fmpz_clear(e);
    TMP_END;
    return status;
}

int
_arf_set_tiny_recip_rel_err(arf_t res, const arf_t x, int side, slong P, slong prec, arf_rnd_t rnd)
{
    nn_srcptr xd;
    nn_ptr a, q, r;
    slong xn, wn, eb, eblo, ebhi;
    fmpz_t e;
    int status, exact;
    TMP_INIT;

    FLINT_ASSERT(P >= prec + 3);

    ARF_GET_MPN_READONLY(xd, xn, x);

    /* a looser bound is as good once the window has a couple of limbs */
    P = FLINT_MIN(P, prec + FLINT_BITS * (xn + 2));

    /* Q = floor(B^(wn + xn) / M) at wn + 1 limbs (Q <= 2 B^wn), the
       remainder telling whether 1/M = Q B^-wn exactly; the term
       |delta / M| < 2^(1 - P) in units of B^-wn is below 2^eb.  (The
       window of bits between the rounding point and 2^eb lies at the
       same positions of the expansion of 1/M for every wn, so there is
       nothing to retry: when it does not decide, the caller evaluates
       the function.) */
    wn = (P + FLINT_BITS - 1) / FLINT_BITS + 1;
    eb = FLINT_BITS * wn + 1 - P;

    TMP_START;
    a = TMP_ALLOC((2 * wn + 2 * xn + 5) * sizeof(ulong));
    q = a + wn + xn + 1;
    r = q + wn + 2;
    _flint_mpn_zero_small(a, wn + xn);
    a[wn + xn] = 1;
    flint_mpn_tdiv_qr(q, r, a, wn + xn + 1, xd, xn);
    exact = flint_mpn_zero_p(r, xn);

    /* V = (1/M)(1 + delta) B^wn against Q: 1/M B^wn in [Q, Q + 1) */
    if (side > 0)
    {
        eblo = -1;
        ebhi = exact ? eb : eb + 1;
    }
    else
    {
        eblo = eb;
        ebhi = exact ? -1 : 0;
    }

    /* the exponent of the unit, -E - FLINT_BITS wn (res may alias x) */
    fmpz_init(e);
    fmpz_neg(e, ARF_EXPREF(x));
    fmpz_sub_ui(e, e, FLINT_BITS * wn);

    if (!COEFF_IS_MPZ(*e))
        status = _arf_set_mpn_err(res, q, wn + 2, *e, eblo, ebhi, ARF_SGNBIT(x), prec, rnd);
    else
    {
        status = _arf_set_mpn_err(res, q, wn + 2, 0, eblo, ebhi, ARF_SGNBIT(x), prec, rnd);
        if (status == 1)
            fmpz_add(ARF_EXPREF(res), ARF_EXPREF(res), e);
    }

    fmpz_clear(e);
    TMP_END;
    return status;
}

/* Quotients of fixed-point results ******************************************/

int
_arf_set_mpn_err_quot(arf_t res, nn_srcptr y, slong len, slong ue, ulong err, nn_srcptr d, slong dn, slong de, int invert, int sgnbit, slong wn, slong prec, arf_rnd_t rnd)
{
    nn_ptr t, q;
    slong te, fix, qe;
    ulong et, eq;
    int status;
    TMP_INIT;

    FLINT_ASSERT(wn >= 1);
    FLINT_ASSERT(d == NULL || (dn >= 1 && (d[dn - 1] >> (FLINT_BITS - 1)) == 1));

    while (len > 0 && y[len - 1] == 0)
        len--;
    if (len == 0)
        return ARF_ELEM_UNCERTAIN;

    TMP_START;
    t = TMP_ALLOC((2 * wn + 3) * sizeof(ulong));
    q = t + wn + 1;

    /* Y 2^ue = T 2^te within et units, T normalised to wn limbs: a
       relative error below 2 et / B^wn (T >= B^wn / 2); when Y is short
       of wn limbs by more than the scaled error can take (a word: the
       left shift of a short Y, FLINT_BITS - 4 bits at most), one limb
       less, which keeps all but a few of the bits of Y */
    for (;;)
    {
        et = _mp_real_fp_normalize(t, &te, y, len, ue, err, wn);
        if (et != UWORD_MAX && et < (UWORD(1) << (FLINT_BITS - 4)))
            break;
        if (wn == 1)
        {
            TMP_END;
            return ARF_ELEM_UNCERTAIN;
        }
        wn--;
        q = t + wn + 1;
    }

    /* the quotient within an ulp of its magnitude, normalised in
       [B^wn / 2, B^wn) with the exponent correction fix: the relative
       error of T gives at most 2 et units of it, plus one for the
       quotient (and one more for luck) */
    if (!invert)
    {
        if (d == NULL)
        {
            flint_mpn_copyi(q, t, wn);
            fix = 0;
            eq = et;
        }
        else
        {
            /* (T B^-wn) / (D B^-dn) with the exponents te + FLINT_BITS wn
               and de */
            _arf_mpn_div_approx(q, &fix, t, wn, d, dn, FLINT_BITS * wn);
            eq = 2 * et + 2;
        }
        qe = fix + te + FLINT_BITS * wn - de;
    }
    else
    {
        if (d == NULL)
            _arf_mpn_inv_approx(q, &fix, t, wn, FLINT_BITS * wn);
        else
            _arf_mpn_div_approx(q, &fix, d, dn, t, wn, FLINT_BITS * wn);
        eq = 2 * et + 2;
        qe = fix - te - FLINT_BITS * wn + ((d == NULL) ? 0 : de);
    }

    /* the result Q B^-wn 2^qe, the unit FLINT_BITS wn below */
    status = _arf_set_mpn_err(res, q, wn, qe - FLINT_BITS * wn, ARF_ELEM_EB(eq), ARF_ELEM_EB(eq), sgnbit, prec, rnd);
    TMP_END;
    return status;
}

/* The short series of the small arguments: f(x) = x + c x^3 (c = -1/6,
   1/3, -1/3, 1/6, -1/6, 1/3 for sin, tan, atan, sinh, asinh, atanh)
   resp. 1 + c x^2 (c = -1/2, 1/2 for cos, cosh) for the exact
   |x| = (d, xn) 2^(e - FLINT_BITS xn), in fixed point at
   n = nk + 2 limbs: for the odd functions Y = |x| 2^-e B^n and the term
   from the top nt limbs of Y (enough for the term, 2^(2e) below Y),
   the result (Y + c Y^3 2^(2e)) 2^(e - FLINT_BITS n) within a few
   units; for the even ones B^n + c Y^2 2^(2e), relative to 1.  The
   caller bounds the tail.  Here n is the number of limbs of the frame
   (nk + 2 for the kernels' nk); Y (n + 2 limbs) receives the n + 1
   limbs of the result (the integral limb Y[n]), with its unit exponent
   ue = (e or 0) - FLINT_BITS n and the error counts elo, ehi; returns
   whether the function is even (the result then carries no sign). */
int
_arf_small_arg_series_frame(nn_ptr Y, slong * ue, ulong * elo, ulong * ehi, nn_srcptr d, slong xn, slong e, int kind, slong n)
{
    nn_ptr P, T;
    slong nt, ntn, top, s, sl, sb;
    ulong divisor;
    int even = (kind == ARF_SMALL_SERIES_COS || kind == ARF_SMALL_SERIES_COSH);
    int negative = (kind == ARF_SMALL_SERIES_SIN || kind == ARF_SMALL_SERIES_ATAN
                    || kind == ARF_SMALL_SERIES_ASINH || kind == ARF_SMALL_SERIES_COS);
    int log1p = (kind == ARF_SMALL_SERIES_LOG1P_POS || kind == ARF_SMALL_SERIES_LOG1P_NEG);
    nn_ptr S = NULL;
    TMP_INIT;

    FLINT_ASSERT(n >= 2);
    FLINT_ASSERT(xn >= 1 && (d[xn - 1] >> (FLINT_BITS - 1)) == 1);

    s = -2 * e;     /* the shift of the term (at least (prec + 2) / 2) */
    FLINT_ASSERT(s >= 2);

    /* limbs for the term: a relative error 2^(1 - FLINT_BITS ntn) on
       the term 2^-s below the unit of Y stays below a unit */
    ntn = n - s / FLINT_BITS;
    ntn = FLINT_MAX(ntn, 1);
    nt = FLINT_MIN(ntn, xn);    /* (an exact Y when xn is smaller) */

    TMP_START;
    P = TMP_ALLOC((3 * (n + 1) + 2 * (n + 1) + (log1p ? 3 * (n + 1) : 0)) * sizeof(ulong));
    T = P + 3 * (n + 1);
    if (log1p)
        S = T + 2 * (n + 1);     /* (S, n + 1) and the square (2 n + 2) behind it */

    /* Y = |x| 2^-e B^n: the mantissa, truncated or zero-extended */
    if (xn >= n)
        flint_mpn_copyi(Y, d + xn - n, n);
    else
    {
        flint_mpn_zero(Y, n - xn);
        flint_mpn_copyi(Y + n - xn, d, xn);
    }
    Y[n] = 0;

    /* P = the top nt limbs of Y squared (2 nt limbs), then (odd) its top
       nt2 limbs times them again (nt2 + nt limbs): top limbs of P are the
       term's mantissa before the shift (the truncations each within the
       relative error above) */
    flint_mpn_sqr(P, Y + n - nt, nt);
    top = 2 * nt;
    if (log1p)
    {
        /* the term x^2/2 in the units of Y is (y^2 B^n) 2^(e - 1),
           shifted by 1 - e >= 2: only that far below the unit, so from
           its own number of limbs of Y (as ntn for the cube) */
        slong sh = 1 - e, shl = sh / FLINT_BITS, shb = sh % FLINT_BITS;
        slong nts = FLINT_MIN(FLINT_MAX(n - sh / FLINT_BITS, 1), xn), tops = 2 * nts;
        nn_ptr Q = S + n + 1;
        flint_mpn_sqr(Q, Y + n - nts, nts);
        if (tops >= n)
            flint_mpn_copyi(S, Q + tops - n, n);
        else
        {
            flint_mpn_zero(S, n - tops);
            flint_mpn_copyi(S + n - tops, Q, tops);
        }
        if (shl >= n)
            flint_mpn_zero(S, n);
        else
        {
            if (shb)
                mpn_rshift(S, S + shl, n - shl, shb);
            else
                flint_mpn_copyi(S, S + shl, n - shl);
            flint_mpn_zero(S + n - shl, shl);
        }
    }
    if (!even)
    {
        slong nt2 = FLINT_MIN(2 * nt, ntn);
        flint_mpn_copyi(T, P + 2 * nt - nt2, nt2);
        flint_mpn_mul(P, T, nt2, Y + n - nt, nt);
        top = nt2 + nt;
    }
    /* the top n limbs of P as the term B^n 2^-s (zero-extended) */
    if (top >= n)
        flint_mpn_copyi(T, P + top - n, n);
    else
    {
        flint_mpn_zero(T, n - top);
        flint_mpn_copyi(T + n - top, P, top);
    }
    T[n] = 0;

    /* T >>= s */
    sl = s / FLINT_BITS;
    sb = s % FLINT_BITS;
    if (sl >= n)
        flint_mpn_zero(T, n);
    else
    {
        if (sb)
            mpn_rshift(T, T + sl, n - sl, sb);
        else
            flint_mpn_copyi(T, T + sl, n - sl);
        flint_mpn_zero(T + n - sl, sl);
    }

    /* T /= 6, 3 or 2 */
    divisor = (kind == ARF_SMALL_SERIES_SIN || kind == ARF_SMALL_SERIES_SINH
               || kind == ARF_SMALL_SERIES_ASINH) ? 6
            : (kind == ARF_SMALL_SERIES_COS || kind == ARF_SMALL_SERIES_COSH) ? 2 : 3;
    mpn_divrem_1(T, 0, T, n, divisor);

    if (even)
    {
        /* Y = B^n */
        flint_mpn_zero(Y, n);
        Y[n] = 1;
    }

    /* the result Y -+ T, within 12 units on either side (the truncation
       of Y, three times two units through the cube, that of the square,
       the shift and the division; measured within [-2.4, 3.4] units) */
    if (log1p)
    {
        /* |log(1 + x)| = |x| -+ x^2/2 + |x|^3/3 - ... for x > 0 resp.
           x < 0 (then every term adds): the square's term within 5 units
           below (the truncations of its limbs, the square and the
           shift), so the value lies within [-17, +17] units */
        mpn_add(Y, Y, n + 1, T, n);
        if (kind == ARF_SMALL_SERIES_LOG1P_POS)
            mpn_sub(Y, Y, n + 1, S, n);
        else
            mpn_add(Y, Y, n + 1, S, n);
        *elo = 20;
        *ehi = 20;
    }
    else
    {
        if (negative)
            mpn_sub(Y, Y, n + 1, T, n);
        else
            mpn_add(Y, Y, n + 1, T, n);
        *elo = 12;
        *ehi = 12;
    }

    *ue = (even ? 0 : e) - FLINT_BITS * n;
    TMP_END;
    return even;
}

int
_arf_small_arg_series(arf_t res, nn_srcptr d, slong xn, slong e, int sgnbit, int kind, slong tail, slong prec, arf_rnd_t rnd)
{
    nn_ptr Y;
    slong n = _arf_elem_limbs(prec, rnd, 0) + 2, ue, eb, teb;
    ulong elo, ehi;
    int even, status;
    TMP_INIT;

    TMP_START;
    Y = TMP_ALLOC((n + 2) * sizeof(ulong));
    even = _arf_small_arg_series_frame(Y, &ue, &elo, &ehi, d, xn, e, kind, n);

    /* the tail 2^(N e + tail) in units 2^ue (N = 4 for the even
       functions and log1p, 5 otherwise), on both sides with the
       truncations (below 2^5 units) */
    teb = ((even || kind == ARF_SMALL_SERIES_LOG1P_POS || kind == ARF_SMALL_SERIES_LOG1P_NEG) ? 4 : 5) * e + tail - ue;
    eb = FLINT_MAX(teb, 5);
    status = _arf_set_mpn_err(res, Y, n + 1, ue, eb, eb, even ? 0 : sgnbit, prec, rnd);
    TMP_END;
    return status;
}
