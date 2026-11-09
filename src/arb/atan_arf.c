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
    atan(x) of an exact x into a ball (arb/elem.h): tiny arguments as
    x within x^3, huge ones as pi/2 within 1/x, the register fast paths
    at one and two limbs, else the kernel _mp_real_atan_fp.
*/

/* +-pi/2 - eps with 0 < eps < 1/|x| <= 2^(1 - e) (e the exponent of x) */
static void
_arb_atan_pi2_eps(arb_t z, int neg, const fmpz_t e, slong prec)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_neg(t, e);
    fmpz_add_ui(t, t, 1);
    arb_const_pi(z, prec);
    arb_mul_2exp_si(z, z, -1);
    if (neg)
        arb_neg(z, z);
    arb_add_error_2exp_fmpz(z, t);
    fmpz_clear(t);
}

void
arb_atan_arf(arb_t z, const arf_t x, slong prec)
{
    nn_srcptr d;
    slong xn, e, nk, w;
    ulong err;
    int sgnbit;

    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
            arb_zero(z);
        else if (arf_is_inf(x))
        {
            int neg = arf_is_neg_inf(x);    /* (z may alias x) */
            arb_const_pi(z, prec);
            arb_mul_2exp_si(z, z, -1);
            if (neg)
                arb_neg(z, z);
        }
        else
            arb_indeterminate(z);
        return;
    }

    sgnbit = ARF_SGNBIT(x);
    e = ARF_EXP(x);
    ARF_GET_MPN_READONLY(d, xn, x);

    /* tiny: atan(x) = x within x^3 (any exponent) */
    if (COEFF_IS_MPZ(e) ? (fmpz_sgn(ARF_EXPREF(x)) < 0) : (2 * e <= -prec - 2))
    {
        fmpz_t t;
        int inexact;
        fmpz_init(t);
        fmpz_mul_ui(t, ARF_EXPREF(x), 3);
        inexact = arf_set_round(arb_midref(z), x, prec, ARB_RND);
        mag_one(arb_radref(z));
        mag_mul_2exp_fmpz(arb_radref(z), arb_radref(z), t);
        if (inexact)
            arf_mag_add_ulp(arb_radref(z), arb_radref(z), arb_midref(z), prec);
        fmpz_clear(t);
        return;
    }

    /* huge: pi/2 - 1/x within the error */
    if (COEFF_IS_MPZ(e) || e > prec + 2)
    {
        _arb_atan_pi2_eps(z, sgnbit, ARF_EXPREF(x), prec);
        return;
    }

#if FLINT_BITS == 64
    if (xn == 1 && prec <= FLINT_BITS)
    {
        ulong t2 = 0, t1 = 0, t0 = 0;
        if (_arf_small_atan_1(&t2, &t1, &t0, &err, d[0], e) == 0)
        {
            _arb_small_set_1(z, t2, t1, t0, 0, err, sgnbit, prec);
            return;
        }
    }
    else if (xn <= 2 && prec <= 2 * FLINT_BITS)
    {
        ulong t3 = 0, t2 = 0, t1 = 0, t0 = 0;
        if (_arf_small_atan_2(&t3, &t2, &t1, &t0, &err, d[xn - 1], (xn == 2) ? d[0] : 0, e) == 0)
        {
            _arb_small_set_2(z, t3, t2, t1, t0, 0, err, sgnbit, prec);
            return;
        }
    }
#endif

    /* small (4e <= -prec - 2): atan(x) = x - x^3/3 (arb/elem.c) within
       |x|^5 / 5 < 2^(5e - 2) */
    if (e <= -((prec + 5) / 4))
    {
        _arb_small_arg_series(z, d, xn, e, sgnbit, ARB_SMALL_SERIES_ATAN, prec);
        arb_add_error_2exp_si(z, 5 * e - 2);
        return;
    }

    nk = _arb_elem_limbs(prec);
    {
        nn_ptr y;
        slong wmax = nk + (prec + 2) / (2 * FLINT_BITS) + 3;
        TMP_INIT;
        TMP_START;
        y = TMP_ALLOC((wmax + 2) * sizeof(ulong));
        _mp_real_atan_fp(y, &w, &err, d, xn, e, _arb_elem_bits(prec), wmax);
        _arb_set_mpn_err(z, y, w + 1, -FLINT_BITS * w, err, err, sgnbit, prec);
        TMP_END;
    }
}
