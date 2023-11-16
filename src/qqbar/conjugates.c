/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "qqbar.h"

void
qqbar_conjugates(qqbar_ptr res, const qqbar_t x)
{
    if (qqbar_degree(x) == 1)
        qqbar_set(res, x);
    else
        qqbar_roots_fmpz_poly(res, QQBAR_POLY(x), QQBAR_ROOTS_IRREDUCIBLE);
}

