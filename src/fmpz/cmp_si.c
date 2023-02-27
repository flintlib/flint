/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson

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

int
fmpz_cmp_si(const fmpz_t f, slong g)
{
    fmpz c = *f;

    if (!COEFF_IS_MPZ(c))    /* f is small */
        return c < g ? -1 : c > g;
    else                     /* f is large */
        return flint_mpz_cmp_si(COEFF_TO_PTR(c), g);
}
