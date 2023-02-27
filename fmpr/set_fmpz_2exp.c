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
fmpr_set_fmpz_2exp(fmpr_t x, const fmpz_t man, const fmpz_t exp)
{
    if (fmpz_is_zero(man))
    {
        fmpr_zero(x);
    }
    else
    {
        fmpr_set_fmpz(x, man);
        fmpz_add(fmpr_expref(x), fmpr_expref(x), exp);
    }
}
