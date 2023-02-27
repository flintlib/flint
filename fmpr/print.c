/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpr.h"

void
fmpr_print(const fmpr_t x)
{
    if (fmpr_is_normal(x))
    {
        flint_printf("(");
        fmpz_print(fmpr_manref(x));
        flint_printf(" * 2^");
        fmpz_print(fmpr_expref(x));
        flint_printf(")");
    }
    else
    {
        if (fmpr_is_zero(x)) flint_printf("(0)");
        else if (fmpr_is_pos_inf(x)) flint_printf("(+inf)");
        else if (fmpr_is_neg_inf(x)) flint_printf("(-inf)");
        else flint_printf("(nan)");
    }
}
