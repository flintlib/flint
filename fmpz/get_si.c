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

slong 
fmpz_get_si(const fmpz_t f)
{
    return (!COEFF_IS_MPZ(*f) ? *f : flint_mpz_get_si(COEFF_TO_PTR(*f)));
}
