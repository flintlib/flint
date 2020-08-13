/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

flint_bitcnt_t fmpz_bits(const fmpz_t f)
{
    fmpz d = *f;

    if (!COEFF_IS_MPZ(d)) return FLINT_BIT_COUNT(FLINT_ABS(d));
    else return mpz_sizeinbase(COEFF_TO_PTR(d), 2);
}
