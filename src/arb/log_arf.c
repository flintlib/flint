/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb/elem.h"

#define ARB_LOG_SMOOTH_MIN_PREC 256

/*
    log(x) and log(1 + x) of an exact x into a ball (arb/elem.h): the
    arguments near 1 (resp. tiny) directly from x - 1 (log(1 + t) = t
    within t^2 for |t| <= 1/2), the register fast paths at one and two
    limbs, else the kernels _mp_real_log_fp (any exponent, the bignum
    ones by _mp_real_log_fp_huge) and _mp_real_log1p_fp.
*/

/* log(1 + D) for the exact D = x - 1 with |D| < 2^-z, z >= prec + 2:
   D within D^2 */
static void
_arb_log_near_one(arb_t z, const arf_t x, slong zb, slong prec)
{
    int inexact = arf_sub_ui(arb_midref(z), x, 1, prec, ARB_RND);
    mag_set_ui_2exp_si(arb_radref(z), 1, -2 * zb);
    if (inexact)
        arf_mag_add_ulp(arb_radref(z), arb_radref(z), arb_midref(z), prec);
}

/* log(x) for a bignum exponent E of the mantissa (d, xn) */
static void
_arb_log_arf_huge(arb_t z, nn_srcptr d, slong xn, const fmpz_t E, slong prec)
{
    fmpz_t c;
    nn_ptr carr, y;
    slong cn, len, ue;
    ulong err;
    int cneg, neg;
    TMP_INIT;

    fmpz_init(c);
    fmpz_sub_ui(c, E, 1);
    cneg = fmpz_sgn(c) < 0;
    cn = fmpz_size(c);
    fmpz_abs(c, c);

    TMP_START;
    carr = TMP_ALLOC((cn + _mp_real_log_fp_huge_ylen(_arb_elem_bits(prec), cn)) * sizeof(ulong));
    y = carr + cn;
    fmpz_get_ui_array(carr, cn, c);

    neg = _mp_real_log_fp_huge(y, &len, &ue, &err, d, xn, carr, cn, cneg, _arb_elem_bits(prec));
    _arb_set_mpn_err(z, y, len, ue, err, err, neg, prec);

    TMP_END;
    fmpz_clear(c);
}

/* log n for an integer n > 1 whose prime factors are among the first
   MP_REAL_LOG_P_TAB_NUM primes (2 to 41), from the shared cache of
   their logarithms in mp_real (exact floors at nv fraction limbs: the
   sum sum_j c_j log p_j is short of log n by at most sum_j c_j units):
   returns 0 if n has another factor */
static int
_arb_log_ui_smooth(arb_t z, ulong n, slong prec)
{
    static const unsigned char primes[MP_REAL_LOG_P_TAB_NUM] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41};
    ulong c[MP_REAL_LOG_P_TAB_NUM], csum = 0;
    slong j, nv;
    const _mp_real_tab_struct * T;
    nn_ptr acc;
    TMP_INIT;

    for (j = 0; j < MP_REAL_LOG_P_TAB_NUM; j++)
    {
        c[j] = 0;
        while (n % primes[j] == 0)
        {
            n /= primes[j];
            c[j]++;
        }
        csum += c[j];
    }
    if (n != 1)
        return 0;

    nv = _arb_elem_limbs(prec);
    T = _mp_real_log_primes_ensure(MP_REAL_LOG_P_TAB_NUM, nv);

    TMP_START;
    acc = TMP_ALLOC((nv + 2) * sizeof(ulong));
    flint_mpn_zero(acc, nv + 2);
    for (j = 0; j < MP_REAL_LOG_P_TAB_NUM; j++)
        if (c[j] != 0)
            acc[nv + 1] += mpn_addmul_1(acc, _mp_real_machin_entry(T, j, nv), nv + 1, c[j]);
    _arb_set_mpn_err(z, acc, nv + 2, -FLINT_BITS * nv, 0, csum, 0, prec);
    TMP_END;
    return 1;
}

void
arb_log_arf(arb_t z, const arf_t x, slong prec)
{
    nn_srcptr d;
    slong xn, e, len, ue, zb;
    ulong err;
    int neg;

    if (arf_is_special(x) || ARF_SGNBIT(x))
    {
        if (arf_is_pos_inf(x))
            arb_pos_inf(z);
        else if (arf_is_zero(x))
            arb_neg_inf(z);
        else
            arb_indeterminate(z);
        return;
    }

    e = ARF_EXP(x);
    ARF_GET_MPN_READONLY(d, xn, x);

    if (FLINT_UNLIKELY(COEFF_IS_MPZ(e) || e >= (WORD(1) << (FLINT_BITS - 3)) || e <= -(WORD(1) << (FLINT_BITS - 3))))
    {
        _arb_log_arf_huge(z, d, xn, ARF_EXPREF(x), prec);
        return;
    }

    /* a smooth integer (a few factors among the first 13 primes) at a
       precision where the sum of cached logarithms beats the kernel */
    if (xn == 1 && e >= 2 && e <= FLINT_BITS && prec >= ARB_LOG_SMOOTH_MIN_PREC
        && (e == FLINT_BITS || (d[0] << e) == 0) && _arb_log_ui_smooth(z, d[0] >> (FLINT_BITS - e), prec))
        return;

    /* near 1: the leading zero bits of |x - 1| */
    if (e == 0 || e == 1)
    {
        zb = _arf_log_lzb(d, xn, e);
        if (zb == WORD_MAX)
        {
            arb_zero(z);
            return;
        }
        if (zb >= prec + 2)
        {
            _arb_log_near_one(z, x, zb, prec);
            return;
        }
    }

#if FLINT_BITS == 64
    if (xn == 1 && prec <= FLINT_BITS && (ARF_SMALL_LOG_OK_1(d[0], e)))
    {
        ulong r2 = 0, r1 = 0, r0 = 0;
        neg = _arf_small_log_2(&r2, &r1, &r0, d[0], 0, e);
        _arb_small_set_1(z, r2, r1, r0, 0, ARF_SMALL_LOG_ERR_2, neg, prec);
        return;
    }
    if (xn <= 2 && prec <= 2 * FLINT_BITS && (ARF_SMALL_LOG_OK_2(d[xn - 1], e)))
    {
        ulong r3 = 0, r2 = 0, r1 = 0, r0 = 0;
        neg = _arf_small_log_3(&r3, &r2, &r1, &r0, d[xn - 1], (xn == 2) ? d[0] : 0, 0, e);
        _arb_small_set_2(z, r3, r2, r1, r0, 0, ARF_SMALL_LOG_ERR_3, neg, prec);
        return;
    }
#endif

    {
        nn_ptr y;
        TMP_INIT;
        TMP_START;
        y = TMP_ALLOC(_mp_real_log_fp_ylen(_arb_elem_bits(prec), xn) * sizeof(ulong));
        neg = _mp_real_log_fp(y, &len, &ue, &err, d, xn, e, 0, _arb_elem_bits(prec));
        _arb_set_mpn_err(z, y, len, ue, err, err, neg, prec);
        TMP_END;
    }
}

void
_arb_log1p_arf(arb_t z, const arf_t x, slong prec)
{
    nn_srcptr d;
    slong xn, e, len, ue;
    ulong err;
    int sgnbit, neg;

    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
            arb_zero(z);
        else if (arf_is_pos_inf(x))
            arb_pos_inf(z);
        else
            arb_indeterminate(z);
        return;
    }

    sgnbit = ARF_SGNBIT(x);
    e = ARF_EXP(x);
    ARF_GET_MPN_READONLY(d, xn, x);

    if (FLINT_UNLIKELY(COEFF_IS_MPZ(e) || e >= (WORD(1) << (FLINT_BITS - 3))))
    {
        if (fmpz_sgn(ARF_EXPREF(x)) < 0)
        {
            /* tiny: log(1 + x) = x within x^2 < 2^(2E) (relative) */
            fmpz_t t;
            int inexact;
            fmpz_init(t);
            fmpz_mul_2exp(t, ARF_EXPREF(x), 1);
            inexact = arf_set_round(arb_midref(z), x, prec, ARB_RND);
            mag_one(arb_radref(z));
            mag_mul_2exp_fmpz(arb_radref(z), arb_radref(z), t);
            if (inexact)
                arf_mag_add_ulp(arb_radref(z), arb_radref(z), arb_midref(z), prec);
            fmpz_clear(t);
            return;
        }
        else
        {
            /* x <= -1: nan (-inf at -1); else log(1 + x) = log x within
               a term below the error */
            if (sgnbit)
                arb_indeterminate(z);
            else
            {
                arb_log_arf(z, x, prec);
                /* log(1 + x) - log x < 1/x: far below an ulp */
                mag_add_ui_2exp_si(arb_radref(z), arb_radref(z), 1, -prec - 10);
            }
            return;
        }
    }

    /* tiny: |log(1 + x) - x| <= x^2 for |x| <= 1/2 */
    if (e <= -prec - 2)
    {
        int inexact = arf_set_round(arb_midref(z), x, prec, ARB_RND);
        mag_set_ui_2exp_si(arb_radref(z), 1, 2 * e);
        if (inexact)
            arf_mag_add_ulp(arb_radref(z), arb_radref(z), arb_midref(z), prec);
        return;
    }


    if (sgnbit && e >= 1)
    {
        /* x <= -1 */
        if (e == 1 && d[xn - 1] == LIMB_TOP && flint_mpn_zero_p(d, xn - 1))
            arb_neg_inf(z);
        else
            arb_indeterminate(z);
        return;
    }

#if FLINT_BITS == 64
    if (xn == 1 && prec <= FLINT_BITS)
    {
        ulong r2 = 0, r1 = 0, r0 = 0;
        if (_arf_small_log1p_1(&r2, &r1, &r0, &err, &neg, d[0], e, sgnbit) == 0)
        {
            _arb_small_set_1(z, r2, r1, r0, 0, err, neg, prec);
            return;
        }
    }
    else if (xn <= 2 && prec <= 2 * FLINT_BITS)
    {
        ulong r3 = 0, r2 = 0, r1 = 0, r0 = 0;
        if (_arf_small_log1p_2(&r3, &r2, &r1, &r0, &err, &neg, d[xn - 1], (xn == 2) ? d[0] : 0, e, sgnbit) == 0)
        {
            _arb_small_set_2(z, r3, r2, r1, r0, 0, err, neg, prec);
            return;
        }
    }
#endif

    /* small (3e <= -prec - 2, |x| < 1/4): log(1 + x) = x - x^2/2 + x^3/3
       within |x|^4 / 2 < 2^(4e - 1), by the series in fixed point
       (arf/elem.c) */
    if (e <= -((prec + 4) / 3) && e <= -2)
    {
        _arb_small_arg_series(z, d, xn, e, sgnbit, sgnbit ? ARB_SMALL_SERIES_LOG1P_NEG : ARB_SMALL_SERIES_LOG1P_POS, prec);
        arb_add_error_2exp_si(z, 4 * e - 1);
        return;
    }


    {
        nn_ptr y;
        TMP_INIT;
        TMP_START;
        y = TMP_ALLOC(_mp_real_log1p_fp_ylen(_arb_elem_bits(prec), xn, e) * sizeof(ulong));
        neg = _mp_real_log1p_fp(y, &len, &ue, &err, d, xn, e, sgnbit, _arb_elem_bits(prec));
        _arb_set_mpn_err(z, y, len, ue, err, err, neg, prec);
        TMP_END;
    }
}
