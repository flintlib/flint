/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPQ_IMPL_H
#define FMPQ_IMPL_H

#include <stdio.h>
#include "mpfr.h"
#include "mpn_extras.h"
#include "fmpq.h"

/* Assumes that h is non-zero */
static __inline__
ulong
_fmpz_gcd_ui(const fmpz_t g, ulong h)
{
    if (!COEFF_IS_MPZ(*g))
    {
        ulong ag = FLINT_ABS(*g);
        return mpn_gcd_1(&ag, 1, h);
    }
    else
    {
        mpz_ptr mg = COEFF_TO_PTR(*g);
        return mpn_gcd_1(mg->_mp_d, FLINT_ABS(mg->_mp_size), h);
    }
}

#endif
