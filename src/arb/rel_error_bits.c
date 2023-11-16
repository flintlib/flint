/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

slong
arb_rel_error_bits(const arb_t x)
{
    fmpz_t t;
    slong result;

    /* fast path for small exponents */
    if (ARB_IS_LAGOM(x))
    {
        if (mag_is_zero(arb_radref(x)))
            return -ARF_PREC_EXACT;
        else if (arf_is_special(arb_midref(x)))
            return ARF_PREC_EXACT;
        else
            return MAG_EXP(arb_radref(x)) + 1 - ARF_EXP(arb_midref(x));
    }

    if (mag_is_zero(arb_radref(x)))
    {
        if (arf_is_nan(arb_midref(x)))
            return ARF_PREC_EXACT;
        else
            return -ARF_PREC_EXACT;
    }

    if (arf_is_special(arb_midref(x)) || mag_is_inf(arb_radref(x)))
        return ARF_PREC_EXACT;

    fmpz_init(t);
    fmpz_add_ui(t, MAG_EXPREF(arb_radref(x)), 1);
    result = _fmpz_sub_small(t, ARF_EXPREF(arb_midref(x)));
    fmpz_clear(t);

    return result;
}

slong arb_rel_one_accuracy_bits(const arb_t x)
{
    if (arf_cmpabs_2exp_si(arb_midref(x), -1) < 0)
    {
        arb_t t;
        arf_init(arb_midref(t));
        arf_one(arb_midref(t));
        *arb_radref(t) = *arb_radref(x);
        return arb_rel_accuracy_bits(t);
    }
    else
    {
        return arb_rel_accuracy_bits(x);
    }
}

