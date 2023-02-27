/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpr.h"

void
fmpr_ulp(fmpr_t u, const fmpr_t x, slong prec)
{
    if (fmpr_is_special(x))
    {
        fmpr_abs(u, x);
    }
    else
    {
        fmpz_t e;
        fmpz_init(e);
        fmpz_add_ui(e, fmpr_expref(x), fmpz_bits(fmpr_manref(x)));
        fmpz_sub_ui(e, e, prec);
        fmpz_one(fmpr_manref(u));
        fmpz_set(fmpr_expref(u), e);
        fmpz_clear(e);
    }
}

