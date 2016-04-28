/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

void _fmpq_fprint(FILE * file, const fmpz_t num, const fmpz_t den)
{
    if (fmpz_is_one(den))
    {
        fmpz_fprint(file, num);
    }
    else
    {
        fmpz_fprint(file, num);
        fputc('/', file);
        fmpz_fprint(file, den);
    }
}

void fmpq_fprint(FILE * file, const fmpq_t x)
{
    _fmpq_fprint(file, &x->num, &x->den);
}
