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

static void
arf_div_special(arf_t z, const arf_t x, const arf_t y)
{
    if ((arf_is_zero(x) && !arf_is_zero(y) && !arf_is_nan(y)) ||
        (arf_is_inf(y) && !arf_is_special(x)))
    {
        arf_zero(z);
    }
    else if (arf_is_zero(y) || (arf_is_special(x) && arf_is_special(y)) ||
        arf_is_nan(x) || arf_is_nan(y))
    {
        arf_nan(z);
    }
    else if (arf_sgn(x) == arf_sgn(y))
        arf_pos_inf(z);
    else
        arf_neg_inf(z);
}

/*
    arf_div through the quotient kernel of arf (arf/mpn_div.c): the
    strict modes correctly rounded (the quotient with a guard limb when
    its guard window certifies the rounding, the exact quotient with the
    remainder otherwise), the relaxed modes as the strict modes of
    arf_rnd_relaxed_to_strict (the inexact flag is reliable in all
    modes).
*/
int
arf_div(arf_ptr z, arf_srcptr x, arf_srcptr y, slong prec, arf_rnd_t rnd)
{
    nn_srcptr xd, yd;
    nn_ptr zd;
    slong xn, yn, rn, fix;
    int inexact, sgnbit;
    ARF_KERNEL_DECL

    if (arf_is_special(x) || arf_is_special(y))
    {
        arf_div_special(z, x, y);
        return 0;
    }

    /* Division by a power of two (exact: allowed at ARF_PREC_EXACT) */
    if (ARF_IS_POW2(y))
    {
        fmpz_t t;
        fmpz_init_set(t, ARF_EXPREF(y));

        if (rnd >= ARF_RND_FAST)
            rnd = arf_rnd_relaxed_to_strict(rnd);

        if (ARF_SGNBIT(y))
            inexact = arf_neg_round(z, x, prec, rnd);
        else
            inexact = arf_set_round(z, x, prec, rnd);

        _fmpz_sub2_fast(ARF_EXPREF(z), ARF_EXPREF(z), t, 1);
        fmpz_clear(t);
        return inexact;
    }

    if (prec == ARF_PREC_EXACT)
        flint_throw(FLINT_ERROR, "arf_div: exact precision not supported\n");

    ARF_GET_MPN_READONLY(xd, xn, x);
    ARF_GET_MPN_READONLY(yd, yn, y);
    sgnbit = ARF_SGNBIT(x) ^ ARF_SGNBIT(y);
    rn = (prec + FLINT_BITS - 1) / FLINT_BITS;
    ARF_KERNEL_ALLOC(zd, rn, z, x, y);

    inexact = _arf_mpn_div(zd, &fix, xd, xn, yd, yn, prec, _arf_rnd_resolve(rnd, sgnbit));

    ARF_KERNEL_FINISH(z, zd, rn, sgnbit);
    _fmpz_sub2_fast(ARF_EXPREF(z), ARF_EXPREF(x), ARF_EXPREF(y), fix);
    return inexact;
}
