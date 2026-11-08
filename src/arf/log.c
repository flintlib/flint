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
    arf_log, arf_log2 and arf_log1p (arf/elem.h explains the rounding).

    - exact cases: log 1 = 0, log2 of a power of two, log1p 0 = 0.
    - near 1 with |x - 1| < 2^-(prec + 3): log(1 + D) = D (1 + delta) with
      D = x - 1 exact and delta on the known side (|delta| < |D|):
      _arf_set_tiny_rel_err, no evaluation; log1p likewise for
      |x| < 2^-(prec + 3).
    - the register fast paths of arf/elem.h for mantissas of one and two
      limbs at prec <= FLINT_BITS resp. 2 FLINT_BITS (64-bit), away from 1.
    - the kernels _mp_real_log_fp / _mp_real_log1p_fp at increasing guard
      bits until the rounding is decided.
    - a bignum exponent (or |e - 1| >= 2^(FLINT_BITS - 3)): the kernel
      _mp_real_log_fp_huge, whose c log 2 term carries the exponent as a
      limb array; log2 adds the exponent to the fraction log2(2M)
      directly, and log1p of such x is log x plus a term below a unit.
*/

/* res = (-1)^sgnbit (t, tn) 2^ue + c 2^ue 2^sh (cneg: minus; the
   magnitude of the sum is positive) rounded: the bignum part of log2 */
static int
_arf_log2_huge_sum(arf_t res, nn_ptr t, slong tn, slong ue, ulong err, int tneg,
    nn_srcptr c, slong cn, int cneg, slong sh, slong prec, arf_rnd_t rnd)
{
    nn_ptr u;
    slong un, ls = sh / FLINT_BITS;
    unsigned int bs = sh % FLINT_BITS;
    int status, neg;
    TMP_INIT;

    TMP_START;
    un = FLINT_MAX(tn, ls + cn + 1) + 1;
    u = TMP_ALLOC((un + 1) * sizeof(ulong));
    _flint_mpn_zero_small(u, un);

    /* u = c 2^sh */
    if (bs == 0)
        flint_mpn_copyi(u + ls, c, cn);
    else
        u[ls + cn] = mpn_lshift(u + ls, c, cn, bs);

    /* +- (t, tn): the magnitude |c| 2^sh dominates t < 2^sh */
    if ((tneg != 0) == (cneg != 0))
        mpn_add(u, u, un, t, tn);
    else
        mpn_sub(u, u, un, t, tn);
    neg = cneg;

    status = _arf_set_mpn_err(res, u, un, ue, ARF_ELEM_EB(err), ARF_ELEM_EB(err), neg, prec, rnd);
    TMP_END;
    return status;
}

/* log(x), log2(x) or log1p(x) (plus_one: x then stands for x - 1 + 1, the
   extra term below a unit) for the bignum exponent E of the mantissa
   (d, xn) 2^(E - xn FLINT_BITS) */
static int
_arf_log_huge(arf_t res, nn_srcptr d, slong xn, const fmpz_t E, int base2, int plus_one, slong prec, arf_rnd_t rnd)
{
    fmpz_t c;
    nn_ptr carr;
    slong cn, nk, len, ue, extra;
    ulong err;
    int cneg, neg, status;

    fmpz_init(c);
    fmpz_sub_ui(c, E, 1);
    cneg = fmpz_sgn(c) < 0;
    cn = fmpz_size(c);

    carr = flint_malloc(cn * sizeof(ulong));
    fmpz_abs(c, c);
    fmpz_get_ui_array(carr, cn, c);

    for (extra = 0; ; extra = (extra == 0) ? 32 : 2 * extra)
    {
        nn_ptr y;
        slong pb = _arf_elem_bits(prec, rnd, extra);
        TMP_INIT;

        nk = _mp_real_elem_limbs(pb);
        TMP_START;

        if (!base2)
        {
            y = TMP_ALLOC(_mp_real_log_fp_huge_ylen(pb, cn) * sizeof(ulong));
            neg = _mp_real_log_fp_huge(y, &len, &ue, &err, d, xn, carr, cn, cneg, pb);
            status = _arf_set_mpn_err(res, y, len, ue, ARF_ELEM_EB(err), ARF_ELEM_EB(err + plus_one), neg, prec, rnd);
        }
        else
        {
            /* c + log2(2M), the fraction from the kernel on 2M in [1, 2) */
            y = TMP_ALLOC(_mp_real_log_fp_ylen(pb, xn) * sizeof(ulong));
            neg = _mp_real_log_fp(y, &len, &ue, &err, d, xn, 1, 1, pb);
            if (len == 0)
            {
                _flint_mpn_zero_small(y, nk + 1);
                len = nk + 1;
                ue = -FLINT_BITS * nk;
                err = 0;
            }
            status = _arf_log2_huge_sum(res, y, len, ue, err + plus_one, neg, carr, cn, cneg, -ue, prec, rnd);
        }

        TMP_END;

        if (status != ARF_ELEM_UNCERTAIN)
            break;
    }

    flint_free(carr);
    fmpz_clear(c);
    return status;
}

/* log(1 + D) for D = x - 1 exact, |D| < 2^-z with z >= prec + 3: D (1 +
   delta) with -D/2 < delta < 0 (D > 0) resp. 0 < delta < |D| (D < 0) */
static int
_arf_log_tiny(arf_t res, nn_srcptr d, slong xn, slong e, slong z, slong prec, arf_rnd_t rnd)
{
    arf_t D;
    nn_ptr t;
    slong tn;
    int status;
    TMP_INIT;

    TMP_START;
    t = TMP_ALLOC((xn + 1) * sizeof(ulong));
    arf_init(D);

    if (e == 1)
    {
        /* D = 2M - 1: the fraction limbs of 2M */
        _mp_real_fp_get_fixed(t, xn + 1, xn, d, xn, 1);
    }
    else
    {
        /* |D| = 1 - M */
        mpn_neg(t, d, xn);
    }
    tn = xn;
    while (t[tn - 1] == 0)
        tn--;
    arf_set_mpn(D, t, tn, e == 0);
    arf_mul_2exp_si(D, D, -FLINT_BITS * xn);

    status = _arf_set_tiny_rel_err(res, D, (e == 1) ? -1 : 1, z, prec, rnd);
    arf_clear(D);
    TMP_END;
    return status;
}

static int
_arf_log(arf_t res, const arf_t x, int base2, slong prec, arf_rnd_t rnd)
{
    nn_srcptr d;
    slong xn, e, extra, len, ue, z;
    ulong err;
    int neg, status;

    if (arf_is_special(x) || ARF_SGNBIT(x))
    {
        if (arf_is_zero(x))
            arf_neg_inf(res);
        else if (arf_is_pos_inf(x))
            arf_pos_inf(res);
        else
            arf_nan(res);
        return 0;
    }

    if (prec == ARF_PREC_EXACT)
        flint_throw(FLINT_ERROR, "arf_log: exact precision not supported\n");

    ARF_GET_MPN_READONLY(d, xn, x);

    if (ARF_IS_POW2(x))
    {
        if (base2)
        {
            /* log2(2^(e - 1)) = e - 1 */
            fmpz_t t;
            int inexact;
            fmpz_init(t);
            fmpz_sub_ui(t, ARF_EXPREF(x), 1);
            inexact = arf_set_round_fmpz(res, t, prec, arf_rnd_relaxed_to_strict(rnd));
            fmpz_clear(t);
            return inexact;
        }
        else if (fmpz_is_one(ARF_EXPREF(x)))
        {
            arf_zero(res);
            return 0;
        }
    }

    e = ARF_EXP(x);

    if (FLINT_UNLIKELY(COEFF_IS_MPZ(e) || e - 1 >= (WORD(1) << (FLINT_BITS - 3))
        || e - 1 <= -(WORD(1) << (FLINT_BITS - 3))))
        return _arf_log_huge(res, d, xn, ARF_EXPREF(x), base2, 0, prec, rnd);

    /* near 1: log(1 + D) = D (1 + delta), |delta| < |D| < 2^-z */
    if (e == 0 || e == 1)
    {
        z = _arf_log_lzb(d, xn, e);
        if (z >= prec + 3 && !base2)
        {
            status = _arf_log_tiny(res, d, xn, e, z, prec, rnd);
            if (status != ARF_ELEM_UNCERTAIN)
                return status;
        }
    }

#if FLINT_BITS == 64
    if (!base2)
    {
        if (xn == 1 && prec <= FLINT_BITS && (ARF_SMALL_LOG_OK_1(d[0], e)))
        {
            ulong r2, r1, r0;
            neg = _arf_small_log_2(&r2, &r1, &r0, d[0], 0, e);
            status = _arf_small_set_1(res, r2, r1, r0, 0, ARF_SMALL_LOG_ERR_2, neg, prec, rnd);
            if (status != ARF_ELEM_UNCERTAIN)
                return status;
        }
        else if (xn <= 2 && prec <= 2 * FLINT_BITS && (ARF_SMALL_LOG_OK_2(d[xn - 1], e)))
        {
            ulong r3, r2, r1, r0;
            neg = _arf_small_log_3(&r3, &r2, &r1, &r0, d[xn - 1], (xn == 2) ? d[0] : 0, 0, e);
            status = _arf_small_set_2(res, r3, r2, r1, r0, 0, ARF_SMALL_LOG_ERR_3, neg, prec, rnd);
            if (status != ARF_ELEM_UNCERTAIN)
                return status;
        }
    }
#endif

    for (extra = 0; ; extra = (extra == 0) ? 32 : 2 * extra)
    {
        nn_ptr y;
        TMP_INIT;

        TMP_START;
        y = TMP_ALLOC(_mp_real_log_fp_ylen(_arf_elem_bits(prec, rnd, extra), xn) * sizeof(ulong));
        neg = _mp_real_log_fp(y, &len, &ue, &err, d, xn, e, base2, _arf_elem_bits(prec, rnd, extra));
        FLINT_ASSERT(len != 0);
        status = _arf_set_mpn_err(res, y, len, ue, ARF_ELEM_EB(err), ARF_ELEM_EB(err), neg, prec, rnd);
        TMP_END;

        if (status != ARF_ELEM_UNCERTAIN)
            return status;
    }
}

int
arf_log(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    return _arf_log(res, x, 0, prec, rnd);
}

int
arf_log2(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    return _arf_log(res, x, 1, prec, rnd);
}

int
arf_log1p(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    nn_srcptr d;
    slong xn, e, extra, len, ue;
    ulong err;
    int sgnbit, neg, status;

    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
            arf_zero(res);
        else if (arf_is_pos_inf(x))
            arf_pos_inf(res);
        else
            arf_nan(res);
        return 0;
    }

    if (prec == ARF_PREC_EXACT)
        flint_throw(FLINT_ERROR, "arf_log1p: exact precision not supported\n");

    sgnbit = ARF_SGNBIT(x);
    e = ARF_EXP(x);
    ARF_GET_MPN_READONLY(d, xn, x);

    if (sgnbit && (COEFF_IS_MPZ(e) ? (fmpz_sgn(ARF_EXPREF(x)) > 0) : (e >= 1)))
    {
        /* x <= -1 */
        if (!COEFF_IS_MPZ(e) && e == 1 && ARF_IS_POW2(x))
            arf_neg_inf(res);
        else
            arf_nan(res);
        return 0;
    }

    if (FLINT_UNLIKELY(COEFF_IS_MPZ(e) || e >= (WORD(1) << (FLINT_BITS - 3)) || e < -prec - 2))
    {
        /* log1p(x) = x (1 + delta) with |delta| < |x| */
        if (COEFF_IS_MPZ(e) ? (fmpz_sgn(ARF_EXPREF(x)) < 0) : (e < -prec - 2))
        {
            status = _arf_set_tiny_rel_err(res, x, sgnbit ? 1 : -1,
                COEFF_IS_MPZ(e) ? prec + FLINT_BITS * (xn + 1) : -e, prec, rnd);
            if (status != ARF_ELEM_UNCERTAIN)
                return status;
            FLINT_ASSERT(!COEFF_IS_MPZ(e));
        }
        else
        {
            /* log(x) + log1p(1/x), the latter below 2^(1 - e) */
            return _arf_log_huge(res, d, xn, ARF_EXPREF(x), 0, 1, prec, rnd);
        }
    }

#if FLINT_BITS == 64
    if (xn == 1 && prec <= FLINT_BITS)
    {
        ulong r2, r1, r0;
        if (_arf_small_log1p_1(&r2, &r1, &r0, &err, &neg, d[0], e, sgnbit) == 0)
        {
            status = _arf_small_set_1(res, r2, r1, r0, 0, err, neg, prec, rnd);
            if (status != ARF_ELEM_UNCERTAIN)
                return status;
        }
    }
    else if (xn <= 2 && prec <= 2 * FLINT_BITS)
    {
        ulong r3, r2, r1, r0;
        if (_arf_small_log1p_2(&r3, &r2, &r1, &r0, &err, &neg, d[xn - 1], (xn == 2) ? d[0] : 0, e, sgnbit) == 0)
        {
            status = _arf_small_set_2(res, r3, r2, r1, r0, 0, err, neg, prec, rnd);
            if (status != ARF_ELEM_UNCERTAIN)
                return status;
        }
    }
#endif

    /* small x (3e <= -prec - 2, |x| < 1/4): log(1 + x) = x - x^2/2 + x^3/3
       within |x|^4 / 2 < 2^(4e - 1), by the series in fixed point
       (arf/elem.c) */
    if (e <= -((prec + 4) / 3) && e <= -2)
    {
        status = _arf_small_arg_series(res, d, xn, e, sgnbit,
            sgnbit ? ARF_SMALL_SERIES_LOG1P_NEG : ARF_SMALL_SERIES_LOG1P_POS, -1, prec, rnd);
        if (status != ARF_ELEM_UNCERTAIN)
            return status;
    }

    for (extra = 0; ; extra = (extra == 0) ? 32 : 2 * extra)
    {
        nn_ptr y;
        TMP_INIT;

        TMP_START;
        y = TMP_ALLOC(_mp_real_log1p_fp_ylen(_arf_elem_bits(prec, rnd, extra), xn, e) * sizeof(ulong));
        neg = _mp_real_log1p_fp(y, &len, &ue, &err, d, xn, e, sgnbit, _arf_elem_bits(prec, rnd, extra));
        FLINT_ASSERT(len != 0);
        status = _arf_set_mpn_err(res, y, len, ue, ARF_ELEM_EB(err), ARF_ELEM_EB(err), neg, prec, rnd);
        TMP_END;

        if (status != ARF_ELEM_UNCERTAIN)
            return status;
    }
}
