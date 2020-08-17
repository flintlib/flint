/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"

void _fmpz_vec_scalar_mod_fmpz(fmpz *res, const fmpz *vec, slong len, const fmpz_t p)
{
    slong i;

    for (i = 0; i < len; i++)
        fmpz_mod(res + i, vec + i, p);
}

