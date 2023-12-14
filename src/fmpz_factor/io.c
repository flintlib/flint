/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "fmpz.h"
#include "fmpz_factor.h"

int fmpz_factor_fprint(FILE * fs, const fmpz_factor_t factor)
{
    slong i;
    int res = 0;

    if (factor->sign == 0)
        return (fputc('0', fs) != EOF);

    if (factor->sign == -1)
        res += fwrite("-1 * ", sizeof(char), 2 + 3 * (factor->num != 0), fs);
    else if (factor->num == 0)
        return (fputc('1', fs) != EOF);

    for (i = 0; i < factor->num; i++)
    {
        res += fmpz_fprint(fs, factor->p + i);

        if (factor->exp[i] != UWORD(1))
            res += fprintf(fs, "^" WORD_FMT "u", factor->exp[i]);

        if (i != factor->num - 1)
            res += fwrite(" * ", sizeof(char), 3, fs);
    }

    return res;
}

int fmpz_factor_print(const fmpz_factor_t factor)
{
    return fmpz_factor_fprint(stdout, factor);
}
