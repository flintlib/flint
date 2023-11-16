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
psl2z_is_one(const psl2z_t g)
{
    return fmpz_is_one(&g->a) && fmpz_is_zero(&g->b) &&
            fmpz_is_zero(&g->c) && fmpz_is_one(&g->d);
}

