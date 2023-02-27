/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"

void
fmpz_factor_expand_iterative(fmpz_t n, const fmpz_factor_t factor)
{
    slong i;
    fmpz_t tmp;

    fmpz_set_si(n, factor->sign);

    fmpz_init(tmp);
    for (i = 0; i < factor->num; i++)
    {
        fmpz_pow_ui(tmp, factor->p + i, factor->exp[i]);
        fmpz_mul(n, n, tmp);
    }
    fmpz_clear(tmp);
}
