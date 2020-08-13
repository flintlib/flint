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

void fmpz_clear_readonly(fmpz_t f)
{
    if (COEFF_IS_MPZ(*f))
    {
        __mpz_struct *ptr = COEFF_TO_PTR(*f);

        mpz_init(ptr);
        _fmpz_clear_mpz(*f);
        *f = WORD(0);
    }
}

