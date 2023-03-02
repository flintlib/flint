/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

flint_bitcnt_t
fmpq_height_bits(const fmpq_t x)
{
    flint_bitcnt_t a, b;

    a = fmpz_bits(fmpq_numref(x));
    b = fmpz_bits(fmpq_denref(x));

    return FLINT_MAX(a, b);
}
