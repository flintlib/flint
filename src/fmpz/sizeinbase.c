/*
    Copyright (C) 2010, 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "long_extras.h"

size_t fmpz_sizeinbase(const fmpz_t f, int b)
{
    fmpz d = *f;
    
    if (!COEFF_IS_MPZ(d))
        return z_sizeinbase(d, b);
    else
        return mpz_sizeinbase(COEFF_TO_PTR(d), b);
}

