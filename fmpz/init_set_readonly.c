/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"

void fmpz_init_set_readonly(fmpz_t f, const mpz_t z)
{
    if (z->_mp_size == 1 && z->_mp_d[0] <= COEFF_MAX)
    {
        *f = z->_mp_d[0];
    }
    else if (z->_mp_size == -1 && z->_mp_d[0] <= COEFF_MAX)
    {
        *f = -(z->_mp_d[0]);
    }
    else if (z->_mp_size)
    {
        _fmpz_init_readonly_mpz(f, z);
    }
    else
    {
        *f = WORD(0);
    }
}

