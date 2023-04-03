/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

#define ONE_OVER_LOG2_PHI (1.44042009041255648 + 1e-13)

slong
fmpq_cfrac_bound(const fmpq_t x)
{
    if (fmpz_is_one(fmpq_denref(x)))
        return 1;

    return fmpz_bits(fmpq_denref(x)) * ONE_OVER_LOG2_PHI + 2;
}
