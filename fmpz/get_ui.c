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

ulong
fmpz_get_ui(const fmpz_t f)
{
    if (!COEFF_IS_MPZ(*f))      /*value is small */
        return (*f < WORD(0) ? -*f : *f);
    else                        /* value is large */
        return flint_mpz_get_ui(COEFF_TO_PTR(*f));
}
