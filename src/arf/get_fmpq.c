/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "arf.h"

void
arf_get_fmpq(fmpq_t y, const arf_t x)
{
    if (arf_is_zero(x))
    {
        fmpq_zero(y);
    }
    else if (arf_is_special(x) || !ARF_IS_LAGOM(x))
    {
        flint_throw(FLINT_ERROR, "exception: arf_get_fmpq: cannot convert to rational\n");
    }
    else
    {
        fmpz_t man, exp;
        slong e;

        fmpz_init(man);
        fmpz_init(exp);

        arf_get_fmpz_2exp(man, exp, x);

        e = *exp;

        fmpz_set_ui(fmpq_denref(y), UWORD(1));

        if (e >= 0)
        {
            fmpz_mul_2exp(fmpq_numref(y), man, e);
        }
        else
        {
            fmpz_set(fmpq_numref(y), man);
            fmpz_mul_2exp(fmpq_denref(y), fmpq_denref(y), -e);
        }

        fmpz_clear(man);
        fmpz_clear(exp);
    }
}
