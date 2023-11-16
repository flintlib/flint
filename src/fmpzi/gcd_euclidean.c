/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpzi.h"

void
fmpzi_gcd_euclidean(fmpzi_t res, const fmpzi_t X, const fmpzi_t Y)
{
    fmpzi_t x, y, q, r;

    if (fmpzi_is_zero(X))
    {
        fmpzi_canonicalise_unit(res, Y);
        return;
    }

    if (fmpzi_is_zero(Y))
    {
        fmpzi_canonicalise_unit(res, X);
        return;
    }

    fmpzi_init(x);
    fmpzi_init(y);
    fmpzi_init(q);
    fmpzi_init(r);

    fmpzi_set(x, X);
    fmpzi_set(y, Y);

    while (!fmpzi_is_zero(y))
    {
        fmpzi_divrem(q, r, x, y);
        fmpzi_swap(x, y);
        fmpzi_swap(y, r);
    }

    fmpzi_swap(res, x);
    fmpzi_canonicalise_unit(res, res);

    fmpzi_clear(x);
    fmpzi_clear(y);
    fmpzi_clear(q);
    fmpzi_clear(r);
}
