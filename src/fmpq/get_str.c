/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

char * _fmpq_get_str(char * str, int b, const fmpz_t num, const fmpz_t den)
{
    char *s;

    if (str == NULL)
    {
        str = flint_malloc(fmpz_sizeinbase(num, b) + fmpz_sizeinbase(den, b) + 3);

        if (str == NULL)
        {
            flint_throw(FLINT_ERROR, "Exception (_fmpq_get_str). Not enough memory.\n");
        }
    }

    fmpz_get_str(str, b, num);

    if (!fmpz_is_one(den))
    {
        s = str;
        while (*s != '\0')
            s++;

        *s = '/';
        s++;
        fmpz_get_str(s, b, den);
    }

    return str;
}

char * fmpq_get_str(char * str, int b, const fmpq_t f)
{
    return _fmpq_get_str(str, b, fmpq_numref(f), fmpq_denref(f));
}

