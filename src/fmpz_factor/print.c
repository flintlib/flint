/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"

void
fmpz_factor_print(const fmpz_factor_t factor)
{
    slong i;

    if (factor->sign == 0)
    {
        flint_printf("0");
        return;
    }

    if (factor->sign == -1)
    {
        if (factor->num)
            flint_printf("-1 * ");
        else
            flint_printf("-1");
    }

    for (i = 0; i < factor->num; i++)
    {
        fmpz_print(factor->p + i);

        if (factor->exp[i] != UWORD(1))
            flint_printf("^%wu", factor->exp[i]);

        if (i != factor->num - 1)
            flint_printf(" * ");
    }
}
