/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb/elem.h"

/*
    exp(x) and expm1(x) of an exact x into a ball (arb/elem.h): the tiny
    and huge arguments directly, the register fast paths at one and two
    limbs, else the kernels _mp_real_exp_fp, _mp_real_expm1_fp and
    _mp_real_exp_fp_huge.
*/

/* |x| >= 2^expbound: [0, 2^(-2^expbound)] for x < 0, [0 +- inf] else */
static void
_arb_exp_arf_overflow(arb_t z, slong expbound, int negative, int minus_one, slong prec)
{
    if (!negative)
    {
        arf_zero(arb_midref(z));
        mag_inf(arb_radref(z));
    }
    else
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_set_si(t, -1);
        fmpz_mul_2exp(t, t, expbound);
        arf_one(arb_midref(z));
        mag_one(arb_radref(z));
        arb_mul_2exp_fmpz(z, z, t);
        if (minus_one)
            arb_sub_ui(z, z, 1, prec);
        fmpz_clear(t);
    }
}

/* exp(x) - minus_one for small x by the series 1 + x (N = 2) or
   1 + x + x^2/2 (N = 3) in fixed point at nk fraction limbs: the
   truncations within a few units (the counterpart of _arf_exp_series),
   the tail below 2^(N e) added to the radius */
static void
_arb_exp_series(arb_t z, nn_srcptr d, slong xn, slong e, int sgnbit, int minus_one, int N, slong prec)
{
    nn_ptr A, S, P;
    slong nk, nt;
    ulong elo, ehi;
    TMP_INIT;

    nk = _mp_real_elem_limbs(_arb_elem_bits(prec + (minus_one ? -e : 0)));
    nt = FLINT_MIN(xn, nk + 1);

    TMP_START;
    A = TMP_ALLOC((2 * (nk + 3) + 2 * nt) * sizeof(ulong));
    S = A + nk + 3;
    P = S + nk + 3;

    /* A = |x| B^nk truncated, nk + 2 limbs */
    _mp_real_fp_get_fixed(A, nk + 2, nk, d, xn, e);

    if (N == 3)
    {
        /* S = x^2 / 2 B^nk from the top nt limbs of x (the truncation of
           x contributes below a unit), truncated */
        flint_mpn_sqr(P, d + xn - nt, nt);
        _mp_real_fp_get_fixed(S, nk + 2, nk, P, 2 * nt, 2 * e - 1);

        if (!sgnbit)
        {
            /* 1 + x + x^2/2 below the value by the truncations (below
               three units) */
            mpn_add_n(A, A, S, nk + 2);
            elo = 0;
            ehi = 3;
        }
        else
        {
            /* 1 - |x| + x^2/2: the truncation of |x| raises, that of
               x^2/2 lowers */
            mpn_sub_n(A, S, A, nk + 2);     /* S - |x|, negative */
            elo = 2;
            ehi = 1;
        }
    }
    else
    {
        if (!sgnbit)
        {
            elo = 0;
            ehi = 1;
        }
        else
        {
            /* 1 - |x|: the truncation raises by less than a unit */
            mpn_neg(A, A, nk + 2);
            elo = 1;
            ehi = 0;
        }
    }

    if (!minus_one)
    {
        /* + 1 (the two's complements above wrap into 1 - |x| ...) */
        mpn_add_1(A + nk, A + nk, 2, 1);
        _arb_set_mpn_err(z, A, nk + 2, -FLINT_BITS * nk, elo, ehi, 0, prec);
    }
    else if (!sgnbit)
    {
        _arb_set_mpn_err(z, A, nk + 2, -FLINT_BITS * nk, elo, ehi, 0, prec);
    }
    else
    {
        /* |expm1(x)| = |x| - x^2/2 + ...: the negation of A above, with
           the bounds swapped */
        mpn_neg(A, A, nk + 2);
        _arb_set_mpn_err(z, A, nk + 2, -FLINT_BITS * nk, ehi, elo, 1, prec);
    }

    /* the tail: x^N / N! < 2^(N e) */
    arb_add_error_2exp_si(z, N * e);

    TMP_END;
}

void
arb_exp_arf(arb_t z, const arf_t x, slong prec, int minus_one, slong maglim)
{
    nn_srcptr d;
    slong xn, e, nk, k, pb;
    ulong err;
    int sgnbit;

    if (arf_is_special(x))
    {
        if (minus_one)
        {
            if (arf_is_zero(x))
                arb_zero(z);
            else if (arf_is_pos_inf(x))
                arb_pos_inf(z);
            else if (arf_is_neg_inf(x))
                arb_set_si(z, -1);
            else
                arb_indeterminate(z);
        }
        else
        {
            if (arf_is_zero(x))
                arb_one(z);
            else if (arf_is_pos_inf(x))
                arb_pos_inf(z);
            else if (arf_is_neg_inf(x))
                arb_zero(z);
            else
                arb_indeterminate(z);
        }
        return;
    }

    sgnbit = ARF_SGNBIT(x);
    e = ARF_EXP(x);
    ARF_GET_MPN_READONLY(d, xn, x);

    if (FLINT_UNLIKELY(COEFF_IS_MPZ(e)))
    {
        if (fmpz_sgn(ARF_EXPREF(x)) > 0)
            _arb_exp_arf_overflow(z, maglim, sgnbit, minus_one, prec);
        else
        {
            /* |exp(x) - (1 + x)| <= x^2 */
            fmpz_t t;
            int inexact;
            fmpz_init(t);
            fmpz_mul_2exp(t, ARF_EXPREF(x), 1);
            inexact = arf_add_ui(arb_midref(z), x, minus_one ? 0 : 1, prec, ARB_RND);
            mag_one(arb_radref(z));
            mag_mul_2exp_fmpz(arb_radref(z), arb_radref(z), t);
            if (inexact)
                arf_mag_add_ulp(arb_radref(z), arb_radref(z), arb_midref(z), prec);
            fmpz_clear(t);
        }
        return;
    }

    /* tiny: |exp(x) - 1| <= 2|x| */
    if (!minus_one && e < -prec - 4)
    {
        arf_one(arb_midref(z));
        mag_set_ui_2exp_si(arb_radref(z), 1, e + 1);
        return;
    }
    /* |exp(x) - (1 + x)| <= x^2 */
    if (e < (minus_one ? -prec - 4 : -(prec / 2) - 4))
    {
        int inexact = arf_add_ui(arb_midref(z), x, minus_one ? 0 : 1, prec, ARB_RND);
        mag_set_ui_2exp_si(arb_radref(z), 1, 2 * e);
        if (inexact)
            arf_mag_add_ulp(arb_radref(z), arb_radref(z), arb_midref(z), prec);
        return;
    }

    if (e > maglim)
    {
        _arb_exp_arf_overflow(z, maglim, sgnbit, minus_one, prec);
        return;
    }

    /* short series for small x (expm1 needs the relative accuracy
       x^2/6 of the three-term series) */
    if (!minus_one)
    {
        if (2 * e + prec + 2 <= 0)
        {
            _arb_exp_series(z, d, xn, e, sgnbit, 0, 2, prec);
            return;
        }
        else if (3 * e + prec + 3 <= 0)
        {
            _arb_exp_series(z, d, xn, e, sgnbit, 0, 3, prec);
            return;
        }
    }
    else if (2 * e + prec + 4 <= 0)
    {
        _arb_exp_series(z, d, xn, e, sgnbit, 1, 3, prec);
        return;
    }

#if FLINT_BITS == 64
    if (xn == 1 && prec <= FLINT_BITS && e >= (minus_one ? -40 : -FLINT_BITS) && e <= MP_REAL_EXP_FP_MAX_EXP)
    {
        ulong y2, y1, y0;
        int sgn = 0;

        if (minus_one)
            k = _arf_small_expm1_1(&y2, &y1, &y0, &sgn, d[0], e, sgnbit);
        else
            k = _arf_small_exp_1(&y2, &y1, &y0, d[0], e, sgnbit);
        _arb_small_set_1(z, y2, y1, y0, k, ARF_SMALL_EXP_ERR_2 + minus_one, sgn, prec);
        return;
    }
    if (xn <= 2 && prec <= 2 * FLINT_BITS && e >= (minus_one ? -40 : -FLINT_BITS) && e <= MP_REAL_EXP_FP_MAX_EXP)
    {
        ulong y3, y2, y1, y0;
        int sgn = 0;

        if (minus_one)
            k = _arf_small_expm1_2(&y3, &y2, &y1, &y0, &sgn, d[xn - 1], (xn == 2) ? d[0] : 0, e, sgnbit);
        else
            k = _arf_small_exp_2(&y3, &y2, &y1, &y0, d[xn - 1], (xn == 2) ? d[0] : 0, e, sgnbit);
        _arb_small_set_2(z, y3, y2, y1, y0, k, ARF_SMALL_EXP_ERR_3 + minus_one, sgn, prec);
        return;
    }
#endif

    /* (expm1 of a small x keeps its relative accuracy with -e more bits) */
    pb = _arb_elem_bits(prec + ((minus_one && e <= 0) ? -e : 0));
    nk = _mp_real_elem_limbs(pb);

    if (e <= MP_REAL_EXP_FP_MAX_EXP)
    {
        nn_ptr y;
        TMP_INIT;

        TMP_START;
        y = TMP_ALLOC((nk + 2) * sizeof(ulong));

        if (!minus_one)
        {
            k = _mp_real_exp_fp(y, &err, d, xn, e, sgnbit, pb);
            _arb_set_mpn_err(z, y, nk + 1, k - FLINT_BITS * nk, err, err, 0, prec);
        }
        else
        {
            ulong elo, ehi;
            k = _mp_real_expm1_fp(y, &elo, &ehi, d, xn, e, sgnbit, pb);
            _arb_set_mpn_err(z, y, nk + 1, k - FLINT_BITS * nk, elo, ehi, sgnbit, prec);
        }
        TMP_END;
    }
    else
    {
        /* the huge exponential with its bignum exponent; expm1 for
           x > 0 within one more unit, for x < 0 it is -1 + exp(x) with
           exp(x) < 2^(-2^(FLINT_BITS - 3)) */
        nn_ptr y, q;
        slong qn;
        fmpz_t qz;
        TMP_INIT;

        if (minus_one && sgnbit)
        {
            fmpz_t t;
            fmpz_init(t);
            fmpz_set_si(t, -1);
            fmpz_mul_2exp(t, t, e - 1);
            arf_set_si(arb_midref(z), -1);
            mag_one(arb_radref(z));
            mag_mul_2exp_fmpz(arb_radref(z), arb_radref(z), t);
            fmpz_clear(t);
            return;
        }

        TMP_START;
        y = TMP_ALLOC((nk + 2 + (e + FLINT_BITS - 1) / FLINT_BITS + 2) * sizeof(ulong));
        q = y + nk + 2;
        sgnbit = _mp_real_exp_fp_huge(y, &err, q, &qn, d, xn, e, sgnbit, pb);
        _arb_set_mpn_err(z, y, nk + 1, -FLINT_BITS * nk, err + minus_one, err, 0, prec);

        fmpz_init(qz);
        fmpz_set_ui_array(qz, q, qn);
        if (sgnbit)
            fmpz_neg(qz, qz);
        arb_mul_2exp_fmpz(z, z, qz);
        fmpz_clear(qz);
        TMP_END;
    }
}
