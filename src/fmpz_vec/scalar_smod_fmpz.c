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

void _fmpz_vec_scalar_smod_fmpz(fmpz *res, const fmpz *vec, slong len, const fmpz_t p)
{
    slong i;
    fmpz_t pdiv2;

    fmpz_init(pdiv2);
    fmpz_fdiv_q_2exp(pdiv2, p, 1);

    for (i = 0; i < len; i++)
    {
        fmpz_mod(res + i, vec + i, p);

        if (fmpz_cmp(res + i, pdiv2) > 0)
        {
            fmpz_sub(res + i, res + i, p);
        }
    }

    fmpz_clear(pdiv2);
}

