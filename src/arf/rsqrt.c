/*
    Copyright (C) 2014, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "arf.h"

/*
    arf_rsqrt through the root kernel of arf (arf/mpn_div.c, see arf_div).
*/
int
arf_rsqrt(arf_ptr z, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    nn_srcptr xd;
    nn_ptr zd;
    slong xn, rn, fix;
    int t, inexact;
    fmpz_t e;
    ARF_KERNEL_DECL

    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
            arf_pos_inf(z);
        else if (arf_is_pos_inf(x))
            arf_zero(z);
        else
            arf_nan(z);
        return 0;
    }

    if (ARF_SGNBIT(x))
    {
        arf_nan(z);
        return 0;
    }

    if (ARF_IS_POW2(x) && fmpz_is_odd(ARF_EXPREF(x)))
    {
        arf_set(z, x);
        fmpz_neg(ARF_EXPREF(z), ARF_EXPREF(z));
        fmpz_cdiv_q_2exp(ARF_EXPREF(z), ARF_EXPREF(z), 1);
        fmpz_add_ui(ARF_EXPREF(z), ARF_EXPREF(z), 1);
        return 0;
    }

    if (prec == ARF_PREC_EXACT)
        flint_throw(FLINT_ERROR, "arf_rsqrt: exact precision not supported\n");

    t = fmpz_is_odd(ARF_EXPREF(x));
    ARF_GET_MPN_READONLY(xd, xn, x);
    rn = (prec + FLINT_BITS - 1) / FLINT_BITS;
    ARF_KERNEL_ALLOC(zd, rn, z, x, x);

    /* the exponent (1 - t) - (e - t) / 2 (e - t is even), the kernel's
       correction relative to it */
    fmpz_init(e);
    if (!COEFF_IS_MPZ(ARF_EXP(x)))
        *e = (1 - t) - ((ARF_EXP(x) - t) >> 1);
    else
    {
        fmpz_sub_ui(e, ARF_EXPREF(x), t);
        fmpz_tdiv_q_2exp(e, e, 1);
        fmpz_neg(e, e);
        fmpz_add_ui(e, e, 1 - t);
    }

    inexact = _arf_mpn_rsqrt(zd, &fix, xd, xn, t, prec, _arf_rnd_resolve(rnd, 0));

    ARF_KERNEL_FINISH(z, zd, rn, 0);
    _fmpz_add_fast(ARF_EXPREF(z), e, fix);
    fmpz_clear(e);
    return inexact;
}
