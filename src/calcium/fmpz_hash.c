/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "calcium.h"

/* slower alternative: fmpz_fdiv_ui(x 1000000007) */

ulong calcium_fmpz_hash(const fmpz_t x)
{
    if (!COEFF_IS_MPZ(*x))
        return *x;
    else
    {
        mpz_ptr z = COEFF_TO_PTR(*x);
        return (z->_mp_size > 0) ? z->_mp_d[0] : -z->_mp_d[0];
    }
}
