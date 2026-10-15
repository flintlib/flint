/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "acb.h"
#include "fexpr.h"
#include "fexpr_builtin.h"

void
fexpr_set_arb_decimal(fexpr_t res, const arb_t x, slong digits)
{
    if (arb_is_zero(x))
    {
        fexpr_zero(res);
    }
    else if (!arb_is_finite(x))
    {
        if (arf_is_nan(arb_midref(x)) || mag_is_inf(arb_radref(x)))
            fexpr_set_symbol_builtin(res, FEXPR_Undefined);
        else
            fexpr_set_arf(res, arb_midref(x));
    }
    else
    {
        fexpr_t Decimal, s;
        char * str;

        fexpr_init(Decimal);
        fexpr_init(s);
        fexpr_set_symbol_builtin(Decimal, FEXPR_Decimal);
        str = arb_get_str(x, FLINT_MAX(digits, 1), ARB_STR_NO_RADIUS);
        fexpr_set_string(s, str);
        flint_free(str);
        fexpr_call1(res, Decimal, s);
        fexpr_clear(Decimal);
        fexpr_clear(s);
    }
}

void
fexpr_set_acb_decimal(fexpr_t res, const acb_t z, slong digits)
{
    if (!acb_is_finite(z) && (arf_is_nan(arb_midref(acb_realref(z))) || arf_is_nan(arb_midref(acb_imagref(z)))
            || mag_is_inf(arb_radref(acb_realref(z))) || mag_is_inf(arb_radref(acb_imagref(z)))))
    {
        fexpr_set_symbol_builtin(res, FEXPR_Undefined);
    }
    else if (arb_is_zero(acb_imagref(z)))
    {
        fexpr_set_arb_decimal(res, acb_realref(z), digits);
    }
    else
    {
        fexpr_t a, b, I;

        fexpr_init(a);
        fexpr_init(b);
        fexpr_init(I);
        fexpr_set_arb_decimal(a, acb_realref(z), digits);
        fexpr_set_arb_decimal(b, acb_imagref(z), digits);
        fexpr_set_symbol_builtin(I, FEXPR_NumberI);
        fexpr_mul(b, b, I);
        if (fexpr_is_zero(a))
            fexpr_swap(res, b);
        else
            fexpr_add(res, a, b);
        fexpr_clear(a);
        fexpr_clear(b);
        fexpr_clear(I);
    }
}
