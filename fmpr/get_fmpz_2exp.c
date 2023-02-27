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
fmpr_get_fmpz_2exp(fmpz_t man, fmpz_t exp, const fmpr_t x)
{
    if (fmpr_is_zero(x))
    {
        fmpz_zero(man);
        fmpz_zero(exp);
    }
    else
    {
        fmpz_set(man, fmpr_manref(x));
        fmpz_set(exp, fmpr_expref(x));
    }
}
