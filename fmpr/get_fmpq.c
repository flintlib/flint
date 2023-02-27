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
fmpr_get_fmpq(fmpq_t y, const fmpr_t x)
{
    if (fmpr_is_zero(x))
    {
        fmpq_zero(y);
    }
    else if (fmpr_is_special(x) || COEFF_IS_MPZ(*fmpr_expref(x)))
    {
        flint_printf("exception: fmpr_get_fmpq: cannot convert to rational\n");
        flint_abort();
    }
    else
    {
        slong exp = *fmpr_expref(x);

        fmpz_set_ui(fmpq_denref(y), UWORD(1));

        if (exp >= 0)
        {
            fmpz_mul_2exp(fmpq_numref(y), fmpr_manref(x), exp);
        }
        else
        {
            fmpz_set(fmpq_numref(y), fmpr_manref(x));
            fmpz_mul_2exp(fmpq_denref(y), fmpq_denref(y), -exp);
        }
    }
}
