/*
    Copyright (C) 2011 Fredrik Johansson

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

void
_fmpz_vec_prod(fmpz_t res, const fmpz * vec, slong len)
{
    if (len <= 1)
    {
        if (len == 1)
            fmpz_set(res, vec);
        else
            fmpz_one(res);
    }
    else if (len <= 3)
    {
        slong i;

        fmpz_mul(res, vec, vec + 1);
        for (i = 2; i < len; i++)
            fmpz_mul(res, res, vec + i);
    }
    else
    {
        slong m = len / 2;
        fmpz_t tmp;
        fmpz_init(tmp);
        _fmpz_vec_prod(res, vec, m);
        _fmpz_vec_prod(tmp, vec + m, len - m);
        fmpz_mul(res, res, tmp);
        fmpz_clear(tmp);
    }
}
