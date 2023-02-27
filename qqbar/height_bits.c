/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "qqbar.h"

slong
qqbar_height_bits(const qqbar_t x)
{
    slong bits;

    bits = fmpz_poly_max_bits(QQBAR_POLY(x));
    return FLINT_ABS(bits);
}
