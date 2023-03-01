/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpzi.h"

void
fmpzi_div_i(fmpzi_t z, const fmpzi_t x)
{
    if (z == x)
    {
        fmpz_swap(fmpzi_realref(z), fmpzi_imagref(z));
        fmpz_neg(fmpzi_imagref(z), fmpzi_imagref(z));
    }
    else
    {
        fmpz_set(fmpzi_realref(z), fmpzi_imagref(x));
        fmpz_neg(fmpzi_imagref(z), fmpzi_realref(x));
    }
}
