/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmpcompat.h"
#include "fmpz.h"

void flint_mpz_clear_readonly(mpz_t z)
{
    _fmpz_clear_readonly_mpz(z);
}

void flint_mpz_init_set_readonly(mpz_t z, const fmpz_t f)
{
    if (COEFF_IS_MPZ(*f))
    {
        *z = *COEFF_TO_PTR(*f);
    }
    else
    {
        flint_mpz_init_set_si(z, *f);
    }
}
