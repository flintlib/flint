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
    arf_exp, arf_expm1 and arf_exp2 (arf/elem.h explains the rounding).

    By the exponent e of x = (-1)^sgnbit M 2^e (|x| < 2^e):

    - e < -prec - 2: exp(x) = 1 + delta with 0 < |delta| < 2^(e + 1), half
      an ulp of 1, and expm1(x) = x (1 + delta) with |delta| < |x|
      (_arf_one_plus_tiny, _arf_set_tiny_rel_err): no evaluation.
    - e <= -(prec + 2) / 2 resp. -(prec + 3) / 3: the series 1 + x resp.
      1 + x + x^2/2 in fixed point with the tail bounded by 2^(2e) resp.
      2^(3e) on the known side (_arf_exp_series): one exact operation,
      decided by the sign of the tail without a retry for x > 0.
    - the register fast paths of arf/elem.h for mantissas of one and
      two limbs at prec <= FLINT_BITS resp. 2 FLINT_BITS (64-bit).
    - the fixed-point kernels _mp_real_exp_fp / _mp_real_expm1_fp at
      increasing guard bits until the rounding is decided (never more
      than one retry in practice).
    - MP_REAL_EXP_FP_MAX_EXP < e <= ARF_EXP_MAGLIM(prec): the kernel
      _mp_real_exp_fp_huge with a bignum quotient of log 2, which
      becomes the exponent of the result.  Beyond the limit, nan.
*/

/* exp(x) - minus_one for small x by the series 1 + x (N = 2) or
   1 + x + x^2/2 (N = 3) in fixed point at nk fraction limbs, with the
   tail bounded by 2^(N e) on the side of the sign of x^N: ARF_ELEM_UNCERTAIN
   when the window does not decide (the two-sided cases) */
static int
_arf_exp_series(arf_t res, nn_srcptr d, slong xn, slong e, int sgnbit, int minus_one, int N, slong prec, arf_rnd_t rnd)
{
    nn_ptr A, S, P;
    slong nk, nt, tb, eblo, ebhi;
    int status;
    TMP_INIT;

    nk = _arf_elem_limbs(prec, rnd, minus_one ? -e : 0);
    nt = FLINT_MIN(xn, nk + 1);

    TMP_START;
    A = TMP_ALLOC((2 * (nk + 3) + 2 * nt) * sizeof(ulong));
    S = A + nk + 3;
    P = S + nk + 3;

    /* A = |x| B^nk truncated, nk + 2 limbs */
    _mp_real_fp_get_fixed(A, nk + 2, nk, d, xn, e);

    /* the tail bound 2^(N e) in units */
    tb = FLINT_MAX(N * e + FLINT_BITS * nk, 0);

    if (N == 3)
    {
        /* S = x^2 / 2 B^nk from the top nt limbs of x (the truncation of
           x contributes below a unit), truncated */
        flint_mpn_sqr(P, d + xn - nt, nt);
        _mp_real_fp_get_fixed(S, nk + 2, nk, P, 2 * nt, 2 * e - 1);

        if (!sgnbit)
        {
            /* 1 + x + x^2/2 below the value by the tail and the
               truncations (below three units) */
            mpn_add_n(A, A, S, nk + 2);
            eblo = -1;
            ebhi = tb + 2;
        }
        else
        {
            /* 1 - |x| + x^2/2: the truncation of |x| raises, that of
               x^2/2 lowers, the tail x^3/6 lowers */
            mpn_sub_n(A, S, A, nk + 2);     /* S - |x|, negative */
            eblo = tb + 1;
            ebhi = 1;
        }
    }
    else
    {
        if (!sgnbit)
        {
            eblo = -1;
            ebhi = tb + 1;
        }
        else
        {
            /* 1 - |x|: the truncation raises by less than a unit, the
               tail raises by less than 2^tb */
            mpn_neg(A, A, nk + 2);
            eblo = 0;
            ebhi = tb;
        }
    }

    if (!minus_one)
    {
        /* + 1 (the two's complements above wrap into 1 - |x| ...) */
        mpn_add_1(A + nk, A + nk, 2, 1);
        status = _arf_set_mpn_err(res, A, nk + 2, -FLINT_BITS * nk, eblo, ebhi, 0, prec, rnd);
    }
    else if (!sgnbit)
    {
        status = _arf_set_mpn_err(res, A, nk + 2, -FLINT_BITS * nk, eblo, ebhi, 0, prec, rnd);
    }
    else
    {
        /* |expm1(x)| = |x| - x^2/2 + ...: the negation of A above, with
           the bounds swapped */
        mpn_neg(A, A, nk + 2);
        status = _arf_set_mpn_err(res, A, nk + 2, -FLINT_BITS * nk, ebhi, eblo, 1, prec, rnd);
    }

    TMP_END;
    return status;
}

int
_arf_exp_huge(arf_t res, nn_srcptr d, slong xn, slong e, int sgnbit, ulong elo, ulong ehi, slong shift, int rsgn, slong prec, arf_rnd_t rnd)
{
    nn_ptr y, q;
    fmpz_t qz;
    slong nk, qn, extra;
    ulong err;
    int qneg, status;
    TMP_INIT;

    if (e > ARF_EXP_MAGLIM(prec))
    {
        arf_nan(res);
        return 0;
    }

    fmpz_init(qz);
    TMP_START;

    for (extra = 0; ; extra = (extra == 0) ? 32 : 2 * extra)
    {
        nk = _arf_elem_limbs(prec, rnd, extra);
        y = TMP_ALLOC((nk + 2 + (e + FLINT_BITS - 1) / FLINT_BITS + 2) * sizeof(ulong));
        q = y + nk + 2;

        qneg = _mp_real_exp_fp_huge(y, &err, q, &qn, d, xn, e, sgnbit, _arf_elem_bits(prec, rnd, extra));

        status = _arf_set_mpn_err(res, y, nk + 1, shift - FLINT_BITS * nk,
            ARF_ELEM_EB(err + elo), ARF_ELEM_EB(err + ehi), rsgn, prec, rnd);

        if (status != ARF_ELEM_UNCERTAIN)
            break;
    }

    /* the exponent: (-1)^qneg q */
    fmpz_set_ui_array(qz, q, qn);
    if (qneg)
        fmpz_neg(qz, qz);
    fmpz_add(ARF_EXPREF(res), ARF_EXPREF(res), qz);

    fmpz_clear(qz);
    TMP_END;
    return status;
}

static int
_arf_exp_expm1(arf_t res, const arf_t x, int minus_one, slong prec, arf_rnd_t rnd)
{
    nn_srcptr d;
    slong xn, e, nk, k, extra, pb;
    int sgnbit, status;

    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
        {
            if (minus_one)
                arf_zero(res);
            else
                arf_one(res);
        }
        else if (arf_is_pos_inf(x))
            arf_pos_inf(res);
        else if (arf_is_neg_inf(x))
        {
            if (minus_one)
                arf_set_si(res, -1);
            else
                arf_zero(res);
        }
        else
            arf_nan(res);
        return 0;
    }

    if (prec == ARF_PREC_EXACT)
        flint_throw(FLINT_ERROR, "arf_exp: exact precision not supported\n");

    sgnbit = ARF_SGNBIT(x);
    e = ARF_EXP(x);
    ARF_GET_MPN_READONLY(d, xn, x);

    /* tiny and huge arguments (including bignum exponents) leave the
       fast paths through one comparison */
    if (FLINT_UNLIKELY(COEFF_IS_MPZ(e) || e > MP_REAL_EXP_FP_MAX_EXP || e < -prec - 2))
    {
        if (COEFF_IS_MPZ(e) ? (fmpz_sgn(ARF_EXPREF(x)) < 0) : (e < -prec - 2))
        {
            if (!minus_one)
                return _arf_one_plus_tiny(res, 0, sgnbit, prec, rnd);

            /* (the bound x of delta always decides the rounding when it
               lies below the last bit of x; a mantissa longer than that
               may leave it undecided, which the kernels handle) */
            status = _arf_set_tiny_rel_err(res, x, sgnbit ? -1 : 1,
                COEFF_IS_MPZ(e) ? prec + FLINT_BITS * (xn + 1) : -e, prec, rnd);
            if (status != ARF_ELEM_UNCERTAIN)
                return status;
            FLINT_ASSERT(!COEFF_IS_MPZ(e));
        }
        else if (COEFF_IS_MPZ(e))
        {
            if (minus_one && sgnbit)
                return _arf_one_plus_tiny(res, 1, 1, prec, rnd);
            arf_nan(res);
            return 0;
        }
        else
        {
            /* -1 + exp(x) with exp(x) < 2^(-2^(FLINT_BITS - 3)); for
               x > 0 the subtraction of 1 lies below the error of Y, whose
               unit is 2^(q - FLINT_BITS nk) with q huge: one more unit */
            if (minus_one && sgnbit)
                return _arf_one_plus_tiny(res, 1, 1, prec, rnd);
            return _arf_exp_huge(res, d, xn, e, sgnbit, minus_one, 0, 0, 0, prec, rnd);
        }
    }

    /* expm1(x) = -1 + exp(x) with exp(x) < 2^-(prec + 2) for
       x <= -(prec + 2) */
    if (minus_one && sgnbit && e > (slong) FLINT_BIT_COUNT(prec + 2))
        return _arf_one_plus_tiny(res, 1, 1, prec, rnd);

#if FLINT_BITS == 64
    if (xn == 1 && prec <= FLINT_BITS && e >= (minus_one ? -40 : -FLINT_BITS))
    {
        ulong y2, y1, y0;
        int sgn = 0;

        if (minus_one)
            k = _arf_small_expm1_1(&y2, &y1, &y0, &sgn, d[0], e, sgnbit);
        else
            k = _arf_small_exp_1(&y2, &y1, &y0, d[0], e, sgnbit);

        status = _arf_small_set_1(res, y2, y1, y0, k, ARF_SMALL_EXP_ERR_2 + minus_one, sgn, prec, rnd);
        if (status != ARF_ELEM_UNCERTAIN)
            return status;
    }
    else if (xn <= 2 && prec <= 2 * FLINT_BITS && e >= (minus_one ? -40 : -2 * FLINT_BITS))
    {
        ulong y3, y2, y1, y0, d1 = d[xn - 1], d0 = (xn == 2) ? d[0] : 0;
        int sgn = 0;

        if (minus_one)
            k = _arf_small_expm1_2(&y3, &y2, &y1, &y0, &sgn, d1, d0, e, sgnbit);
        else
            k = _arf_small_exp_2(&y3, &y2, &y1, &y0, d1, d0, e, sgnbit);

        status = _arf_small_set_2(res, y3, y2, y1, y0, k, ARF_SMALL_EXP_ERR_3 + minus_one, sgn, prec, rnd);
        if (status != ARF_ELEM_UNCERTAIN)
            return status;
    }
#endif

    /* short series for small x (expm1 needs the relative accuracy
       x^2/6 of the three-term series) */
    status = ARF_ELEM_UNCERTAIN;
    if (!minus_one)
    {
        if (2 * e + prec + 2 <= 0)
            status = _arf_exp_series(res, d, xn, e, sgnbit, 0, 2, prec, rnd);
        else if (3 * e + prec + 3 <= 0)
            status = _arf_exp_series(res, d, xn, e, sgnbit, 0, 3, prec, rnd);
    }
    else if (2 * e + prec + 4 <= 0)
        status = _arf_exp_series(res, d, xn, e, sgnbit, 1, 3, prec, rnd);

    if (status != ARF_ELEM_UNCERTAIN)
        return status;

    for (extra = 0; ; extra = (extra == 0) ? 32 : 2 * extra)
    {
        nn_ptr y;
        ulong elo, ehi;
        TMP_INIT;

        pb = _arf_elem_bits(prec, rnd, extra + ((minus_one && e <= 0) ? -e : 0));
        nk = _mp_real_elem_limbs(pb);

        TMP_START;
        y = TMP_ALLOC((nk + 2) * sizeof(ulong));

        if (minus_one)
        {
            k = _mp_real_expm1_fp(y, &elo, &ehi, d, xn, e, sgnbit, pb);
            status = _arf_set_mpn_err(res, y, nk + 1, k - FLINT_BITS * nk,
                ARF_ELEM_EB(elo), ARF_ELEM_EB(ehi), sgnbit, prec, rnd);
        }
        else
        {
            k = _mp_real_exp_fp(y, &elo, d, xn, e, sgnbit, pb);
            status = _arf_set_mpn_err(res, y, nk + 1, k - FLINT_BITS * nk,
                ARF_ELEM_EB(elo), ARF_ELEM_EB(elo), 0, prec, rnd);
        }

        TMP_END;

        if (status != ARF_ELEM_UNCERTAIN)
            return status;
    }
}

int
arf_exp(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    return _arf_exp_expm1(res, x, 0, prec, rnd);
}

int
arf_expm1(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    return _arf_exp_expm1(res, x, 1, prec, rnd);
}

int
arf_exp2(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    nn_srcptr d;
    slong xn, e, nk, k, extra;
    int sgnbit, status;

    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
            arf_one(res);
        else if (arf_is_pos_inf(x))
            arf_pos_inf(res);
        else if (arf_is_neg_inf(x))
            arf_zero(res);
        else
            arf_nan(res);
        return 0;
    }

    if (prec == ARF_PREC_EXACT)
        flint_throw(FLINT_ERROR, "arf_exp2: exact precision not supported\n");

    sgnbit = ARF_SGNBIT(x);
    e = ARF_EXP(x);

    if (FLINT_UNLIKELY(COEFF_IS_MPZ(e) || e > FLINT_BITS - 3 || e < -prec - 2))
    {
        /* 2^x = 1 + delta with |delta| < |x| log 2 */
        if (COEFF_IS_MPZ(e) ? (fmpz_sgn(ARF_EXPREF(x)) < 0) : (e < -prec - 2))
            return _arf_one_plus_tiny(res, 0, sgnbit, prec, rnd);

        if (COEFF_IS_MPZ(e) || e > ARF_EXP_MAGLIM(prec))
        {
            arf_nan(res);
            return 0;
        }

        /* 2^x = 2^k 2^f with k = floor(x) a bignum and f in [0, 1) exact */
        {
            arf_t f;
            fmpz_t kz;
            int inexact;

            arf_init(f);
            fmpz_init(kz);
            arf_floor(f, x);
            arf_get_fmpz(kz, f, ARF_RND_DOWN);
            arf_sub(f, x, f, ARF_PREC_EXACT, ARF_RND_DOWN);
            inexact = arf_exp2(res, f, prec, rnd);
            fmpz_add(ARF_EXPREF(res), ARF_EXPREF(res), kz);
            arf_clear(f);
            fmpz_clear(kz);
            return inexact;
        }
    }

    /* an integer: 2^x exactly */
    if (arf_is_int(x))
    {
        arf_set_ui_2exp_si(res, 1, arf_get_si(x, ARF_RND_DOWN));
        return 0;
    }

    ARF_GET_MPN_READONLY(d, xn, x);

    for (extra = 0; ; extra = (extra == 0) ? 32 : 2 * extra)
    {
        nn_ptr y;
        ulong err;
        TMP_INIT;

        nk = _arf_elem_limbs(prec, rnd, extra);

        TMP_START;
        y = TMP_ALLOC((nk + 2) * sizeof(ulong));
        k = _mp_real_exp2_fp(y, &err, d, xn, e, sgnbit, _arf_elem_bits(prec, rnd, extra));
        status = _arf_set_mpn_err(res, y, nk + 1, k - FLINT_BITS * nk,
            ARF_ELEM_EB(err), ARF_ELEM_EB(err), 0, prec, rnd);
        TMP_END;

        if (status != ARF_ELEM_UNCERTAIN)
            return status;
    }
}
