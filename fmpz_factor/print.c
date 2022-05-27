/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#undef ulong
#define ulong ulongxx /* ensure vendor doesn't typedef ulong */
#include <stdio.h>
#undef ulong
#include "fmpz_factor.h"

void
fmpz_factor_print(const fmpz_factor_t factor)
{
    slong i;

    if (factor->sign == 0)
    {
        printf("0");
        return;
    }

    if (factor->sign == -1)
    {
        if (factor->num)
            printf("-1 * ");
        else
            printf("-1");
    }

    for (i = 0; i < factor->num; i++)
    {
        fmpz_print(factor->p + i);

        if (factor->exp[i] != UWORD(1))
            printf("^" WORD_FMT "u", factor->exp[i]);

        if (i != factor->num - 1)
            printf(" * ");
    }
}
