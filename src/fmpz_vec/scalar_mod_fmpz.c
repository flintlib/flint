/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"

void _fmpz_vec_scalar_mod_fmpz(fmpz *res, const fmpz *vec, slong len, const fmpz_t p)
{
    slong i;

    flint_printf("modulus "); fmpz_print(p); flint_printf("\n\n");
    flint_printf("length %wd\n\n", len);

    for (i = 0; i < len; i++)
    {
        flint_printf("vec[%wd] = ", i); fmpz_print(vec + i); flint_printf("\n\n"); fflush(stdout);
        fmpz_mod(res + i, vec + i, p);
    }
}

