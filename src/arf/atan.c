/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_extras.h"
#include "arf/elem.h"

/*
    arf_atan and arf_atan2 (arf/elem.h explains the rounding).

    - tiny x: atan x = x (1 - delta) with delta < x^2 / 3, on the known
      side (_arf_set_tiny_rel_err); for atan2 the kernel returns the
      ratio of the mantissas with its bounds.
    - the register fast paths of arf/elem.h for mantissas of one and two
      limbs at prec <= FLINT_BITS resp. 2 FLINT_BITS (64-bit).
    - the kernels _mp_real_atan_fp and _mp_real_atan2_fp at increasing
      guard bits (and buffers, for the relative accuracy of small
      results) until the rounding is decided.
    - multiples of pi/4 (the special arguments, |y| = |x|) from the
      tabulated pi/4 at increasing precision.
*/

/* res = (-1)^sgnbit (c pi/4 - v) for c in {1, 2, 3, 4} and v = 0 or
   (below = 1) 0 < v below one unit of every working precision: the floor
   of pi/4 at w limbs is within one ulp below, c times, so V lies in
   (Y - below, Y + c) units */
static int
_arf_pi4_mul(arf_t res, int c, int below, int sgnbit, slong prec, arf_rnd_t rnd)
{
    slong extra, w;
    int status;

    for (extra = 0; ; extra = (extra == 0) ? 32 : 2 * extra)
    {
        nn_ptr y;
        TMP_INIT;

        w = _arf_elem_limbs(prec, rnd, extra);

        TMP_START;
        y = TMP_ALLOC((w + 2) * sizeof(ulong));
        y[w] = mpn_mul_1(y, _mp_real_const_ptr(MP_REAL_CONST_ID_PI4, w), w, c);
        status = _arf_set_mpn_err(res, y, w + 1, -FLINT_BITS * w, below ? ARF_ELEM_EB(1) : -1,
            ARF_ELEM_EB(c), sgnbit, prec, rnd);
        TMP_END;

        if (status != ARF_ELEM_UNCERTAIN)
            return status;
    }
}

int
arf_const_pi(arf_ptr res, slong prec, arf_rnd_t rnd)
{
    if (prec == ARF_PREC_EXACT)
        flint_throw(FLINT_ERROR, "arf_const_pi: exact precision not supported\n");
    return _arf_pi4_mul(res, 4, 0, 0, prec, rnd);
}

int
arf_const_pi_mul_2exp(arf_ptr res, slong k, int neg, slong prec, arf_rnd_t rnd)
{
    int status;
    if (prec == ARF_PREC_EXACT)
        flint_throw(FLINT_ERROR, "arf_const_pi: exact precision not supported\n");
    status = _arf_pi4_mul(res, 4, 0, neg, prec, rnd);
    arf_mul_2exp_si(res, res, k);
    return status;
}

int
arf_atan(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    nn_srcptr d;
    slong xn, e, nk, extra, w;
    ulong err;
    fmpz_t xe;
    int xsgn, status, eshift = 0;

    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
        {
            arf_zero(res);
            return 0;
        }
        if (arf_is_nan(x))
        {
            arf_nan(res);
            return 0;
        }
        return _arf_pi4_mul(res, 2, 0, arf_sgn(x) < 0, prec, rnd);
    }

    if (prec == ARF_PREC_EXACT)
        flint_throw(FLINT_ERROR, "arf_atan: exact precision not supported\n");

    xsgn = ARF_SGNBIT(x);
    e = ARF_EXP(x);
    ARF_GET_MPN_READONLY(d, xn, x);

    if (FLINT_UNLIKELY(COEFF_IS_MPZ(e)))
    {
        if (fmpz_sgn(ARF_EXPREF(x)) < 0)
        {
            /* tiny: evaluated as if the exponent were e (the bound below
               holds), the difference restored on the result */
            e = -prec - FLINT_BITS * (xn + 2);
            eshift = 1;
            fmpz_init_set(xe, ARF_EXPREF(x));
        }
        else
        {
            /* pi/2 - v with 0 < v < 1/|x| < 2^(1 - e) */
            return _arf_pi4_mul(res, 2, 1, xsgn, prec, rnd);
        }
    }

    /* atan x = x (1 - delta), 0 < delta < x^2 / 3 < 2^(2e - 1) */
    if (2 * e <= -prec - 2)
    {
        status = _arf_set_tiny_rel_err(res, x, -1, 1 - 2 * e, prec, rnd);
        if (status != ARF_ELEM_UNCERTAIN)
        {
            if (eshift)
                fmpz_clear(xe);
            return status;
        }
    }

#if FLINT_BITS == 64
    if (xn == 1 && prec <= FLINT_BITS)
    {
        ulong t2, t1, t0;
        if (_arf_small_atan_1(&t2, &t1, &t0, &err, d[0], e) == 0)
        {
            status = _arf_small_set_1(res, t2, t1, t0, 0, err, xsgn, prec, rnd);
            if (status != ARF_ELEM_UNCERTAIN)
                return status;
        }
    }
    else if (xn <= 2 && prec <= 2 * FLINT_BITS)
    {
        ulong t3, t2, t1, t0;
        if (_arf_small_atan_2(&t3, &t2, &t1, &t0, &err, d[xn - 1], (xn == 2) ? d[0] : 0, e) == 0)
        {
            status = _arf_small_set_2(res, t3, t2, t1, t0, 0, err, xsgn, prec, rnd);
            if (status != ARF_ELEM_UNCERTAIN)
                return status;
        }
    }
#endif

    /* small x (4e <= -prec - 2): atan x = x - x^3/3 within |x|^5 / 5 <
       2^(5e - 2), by the series in fixed point (arf/elem.c) */
    if (!eshift && e <= -((prec + 5) / 4))
    {
        status = _arf_small_arg_series(res, d, xn, e, xsgn, ARF_SMALL_SERIES_ATAN, -2, prec, rnd);
        if (status != ARF_ELEM_UNCERTAIN)
            return status;
    }

    for (extra = 0; ; extra = (extra == 0) ? 32 : 2 * extra)
    {
        nn_ptr y;
        slong wmax;
        TMP_INIT;

        nk = _arf_elem_limbs(prec, rnd, extra);
        /* the limbs for the leading zero bits of |x| < 1 (at most
           (prec + 1) / 2 unless a tiny argument fell through above, where
           the retries grow the buffers) */
        wmax = nk + (prec + extra) / (2 * FLINT_BITS) + extra / FLINT_BITS + 2;

        TMP_START;
        y = TMP_ALLOC((wmax + 2) * sizeof(ulong));
        _mp_real_atan_fp(y, &w, &err, d, xn, e, _arf_elem_bits(prec, rnd, extra), wmax);
        status = _arf_set_mpn_err(res, y, w + 1, -FLINT_BITS * w, ARF_ELEM_EB(err), ARF_ELEM_EB(err), xsgn, prec, rnd);
        TMP_END;

        if (status != ARF_ELEM_UNCERTAIN)
        {
            if (eshift)
            {
                /* the exponent relative to the fake e */
                fmpz_add(ARF_EXPREF(res), ARF_EXPREF(res), xe);
                fmpz_sub_si(ARF_EXPREF(res), ARF_EXPREF(res), e);
                fmpz_clear(xe);
            }
            return status;
        }
    }
}

static int
_arf_atan2(arf_t res, const arf_t y, const arf_t x, slong prec, arf_rnd_t rnd)
{
    nn_srcptr nd, dd;
    slong nn, dn, d, nk, extra, len, ue;
    ulong elo, ehi;
    arf_srcptr num, den;
    int c, sigma, ysgn, cmp, status;

    ysgn = ARF_SGNBIT(y);
    cmp = arf_cmpabs(y, x);

    /* |y| = |x|: pi/4 or 3 pi/4 */
    if (cmp == 0)
        return _arf_pi4_mul(res, ARF_SGNBIT(x) ? 3 : 1, 0, ysgn, prec, rnd);

    /* theta = c pi/4 + sigma atan(t), t = num / den < 1 */
    if (cmp < 0)
    {
        num = y;
        den = x;
        c = ARF_SGNBIT(x) ? 4 : 0;
        sigma = ARF_SGNBIT(x) ? -1 : 1;
    }
    else
    {
        num = x;
        den = y;
        c = 2;
        sigma = ARF_SGNBIT(x) ? 1 : -1;
    }

    ARF_GET_MPN_READONLY(nd, nn, num);
    ARF_GET_MPN_READONLY(dd, dn, den);

    /* t = (Mnum / Mden) 2^-d, d >= 0 (saturated) */
    d = _fmpz_sub_small(ARF_EXPREF(den), ARF_EXPREF(num));

#if FLINT_BITS == 64
    if (nn == 1 && dn == 1 && prec <= FLINT_BITS)
    {
        ulong t2, t1, t0;
        if (_arf_small_atan2_1(&t2, &t1, &t0, &elo, nd[0], dd[0], d, c, sigma) == 0)
        {
            status = _arf_small_set_1(res, t2, t1, t0, 0, elo, ysgn, prec, rnd);
            if (status != ARF_ELEM_UNCERTAIN)
                return status;
        }
    }
    else if (nn <= 2 && dn <= 2 && prec <= 2 * FLINT_BITS)
    {
        ulong t3, t2, t1, t0;
        if (_arf_small_atan2_2(&t3, &t2, &t1, &t0, &elo, nd[nn - 1], (nn == 2) ? nd[0] : 0,
                dd[dn - 1], (dn == 2) ? dd[0] : 0, d, c, sigma) == 0)
        {
            status = _arf_small_set_2(res, t3, t2, t1, t0, 0, elo, ysgn, prec, rnd);
            if (status != ARF_ELEM_UNCERTAIN)
                return status;
        }
    }
#endif

    for (extra = 0; ; extra = (extra == 0) ? 32 : 2 * extra)
    {
        nn_ptr t;
        slong wmax, dk;
        fmpz_t shift;
        TMP_INIT;

        nk = _arf_elem_limbs(prec, rnd, extra);
        wmax = nk + (nk + 1) / 2 + 4;

        /* beyond the kernel's tiny thresholds the value of d only shifts
           the exponent of a tiny theta = t (1 - delta) */
        dk = FLINT_BITS * (nk + 2) + 8;
        if (d <= dk)
            dk = d;
        else if (c == 0)
        {
            fmpz_init(shift);
            fmpz_sub(shift, ARF_EXPREF(num), ARF_EXPREF(den));
            fmpz_add_ui(shift, shift, dk);
        }

        TMP_START;
        t = TMP_ALLOC((wmax + 2) * sizeof(ulong));
        _mp_real_atan2_fp(t, &len, &ue, &elo, &ehi, nd, nn, dd, dn, dk, c, sigma, _arf_elem_bits(prec, rnd, extra), wmax);
        status = _arf_set_mpn_err(res, t, len, ue, ARF_ELEM_EB(elo), ARF_ELEM_EB(ehi), ysgn, prec, rnd);
        TMP_END;

        if (d > dk && c == 0)
        {
            if (status != ARF_ELEM_UNCERTAIN)
                fmpz_add(ARF_EXPREF(res), ARF_EXPREF(res), shift);
            fmpz_clear(shift);
        }

        if (status != ARF_ELEM_UNCERTAIN)
            return status;
    }
}

int
arf_atan2(arf_ptr res, arf_srcptr y, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    if (arf_is_special(x) || arf_is_special(y))
    {
        int xs = arf_sgn(x), ys = arf_sgn(y);

        if (arf_is_nan(x) || arf_is_nan(y))
        {
            arf_nan(res);
            return 0;
        }

        if (arf_is_zero(y))
        {
            /* 0 for x >= 0, pi for x < 0 */
            if (xs >= 0)
            {
                arf_zero(res);
                return 0;
            }
            return _arf_pi4_mul(res, 4, 0, 0, prec, rnd);
        }

        /* +-pi/2 for x = 0 or y infinite with x finite */
        if (arf_is_zero(x) || (arf_is_inf(y) && arf_is_finite(x)))
            return _arf_pi4_mul(res, 2, 0, ys < 0, prec, rnd);

        /* x infinite: 0 or pi for finite y, pi/4 or 3 pi/4 for infinite y */
        if (arf_is_finite(y))
        {
            if (xs > 0)
            {
                arf_zero(res);
                return 0;
            }
            return _arf_pi4_mul(res, 4, 0, ys < 0, prec, rnd);
        }
        return _arf_pi4_mul(res, (xs > 0) ? 1 : 3, 0, ys < 0, prec, rnd);
    }

    if (prec == ARF_PREC_EXACT)
        flint_throw(FLINT_ERROR, "arf_atan2: exact precision not supported\n");

    if (res == x || res == y)
    {
        arf_t tx, ty;
        int status;
        arf_init(tx);
        arf_init(ty);
        arf_set(tx, x);
        arf_set(ty, y);
        status = _arf_atan2(res, ty, tx, prec, rnd);
        arf_clear(tx);
        arf_clear(ty);
        return status;
    }

    return _arf_atan2(res, y, x, prec, rnd);
}
