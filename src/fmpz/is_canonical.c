/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"

int _fmpz_is_canonical(const fmpz_t x)
{
    __mpz_struct * z;
    mp_size_t n;

    if (!COEFF_IS_MPZ(*x))
        return 1;

    z = COEFF_TO_PTR(*x);
    n = FLINT_ABS(z->_mp_size);

    if (n == 0)
        return 0;

    if (n == 1)
        return z->_mp_d[0] > (mp_limb_t) COEFF_MAX;

    return z->_mp_d[n - 1] != 0;
}
