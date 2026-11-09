/*
    Copyright (C) 2012-2014, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "mp_real.h"

/* atan(p/q) resp. atanh(p/q): by the binary splitting of mp_real
   (_mp_real_atan_frac_bsplit, _mp_real_atanh_frac_bsplit: 0 < p < q
   given as limb arrays, p/q at most about 0.99) at a high precision
   relative to the ratio, else by the kernels on the fraction, which
   are faster (also for |p/q| close to 1; the atanh of such a ratio is
   indeterminate at |p| >= |q|) */
void
arb_atan_frac_bsplit(arb_t s, const fmpz_t p, const fmpz_t q, int hyperbolic, slong prec)
{
    mp_real_t t;
    nn_srcptr pd, qd;
    slong pn, qn;
    ulong ptmp, qtmp;
    int neg;

    if (fmpz_is_zero(p))
    {
        arb_zero(s);
        return;
    }

    if (fmpz_is_zero(q))
    {
        arb_indeterminate(s);
        return;
    }

    neg = (fmpz_sgn(p) < 0) ^ (fmpz_sgn(q) < 0);

    /* the limbs of |p| and |q| */
    if (COEFF_IS_MPZ(*p))
    {
        pd = COEFF_TO_PTR(*p)->_mp_d;
        pn = FLINT_ABS(COEFF_TO_PTR(*p)->_mp_size);
    }
    else
    {
        ptmp = FLINT_ABS(*p);
        pd = &ptmp;
        pn = 1;
    }
    if (COEFF_IS_MPZ(*q))
    {
        qd = COEFF_TO_PTR(*q)->_mp_d;
        qn = FLINT_ABS(COEFF_TO_PTR(*q)->_mp_size);
    }
    else
    {
        qtmp = FLINT_ABS(*q);
        qd = &qtmp;
        qn = 1;
    }

    /* the kernels of the fraction unless the precision is high for the
       ratio (binary splitting wins above about C (16 / lg(q/p))^2 bits
       with C = 10^4 for atan and 4 10^4 for atanh, lg(q/p) from the bit
       lengths), or |p/q| >= 1/2 (roughly) */
    {
        slong lg = FLINT_BITS * (qn - pn) + (slong) FLINT_BIT_COUNT(qd[qn - 1]) - (slong) FLINT_BIT_COUNT(pd[pn - 1]);
        if (lg <= 1 || (double) prec < (hyperbolic ? 4e4 : 1e4) * (16.0 / lg) * (16.0 / lg))
        {
            arb_fmpz_div_fmpz(s, p, q, prec + 10);
            if (hyperbolic)
                arb_atanh(s, s, prec);
            else
                arb_atan(s, s, prec);
            return;
        }
    }

    mp_real_init(t);
    if (hyperbolic)
        _mp_real_atanh_frac_bsplit(t, pd, pn, qd, qn, (prec + FLINT_BITS - 1) / FLINT_BITS + 1);
    else
        _mp_real_atan_frac_bsplit(t, pd, pn, qd, qn, (prec + FLINT_BITS - 1) / FLINT_BITS + 1);
    mp_real_get_arb(s, t);
    arb_set_round(s, s, prec);
    if (neg)
        arb_neg(s, s);
    mp_real_clear(t);
}
