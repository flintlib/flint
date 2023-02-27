/*
    Copyright (C) 2015 Kushagra Singh

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "mpn_extras.h"
#include "fmpz.h"

int
fmpz_get_mpn(mp_ptr *n, fmpz_t n_in)
{
    mp_limb_t n_size;
    mp_ptr temp;

    n_size = fmpz_size(n_in);
    *n = flint_malloc(n_size * sizeof(mp_limb_t));

    if (n_size <= 1)
    {
        (*n)[0] = fmpz_get_ui(n_in);
        return 1;
    }
    else
    {
        temp = COEFF_TO_PTR(*n_in)->_mp_d;
        flint_mpn_copyi(*n, temp, n_size);
        return n_size;
    }
}
