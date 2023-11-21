/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_modular.h"

int
psl2z_is_correct(const psl2z_t g)
{
    int res;
    fmpz_t t;

    if (fmpz_sgn(&g->c) < 0)
        return 0;

    if (fmpz_is_zero(&g->c) && fmpz_sgn(&g->d) <= 0)
        return 0;

    fmpz_init(t);
    fmpz_mul(t, &g->a, &g->d);
    fmpz_submul(t, &g->b, &g->c);
    res = fmpz_is_one(t);
    fmpz_clear(t);
    return res;
}

