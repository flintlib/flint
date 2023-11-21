/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

/* The clear methods are not defined as inlines since this inflates
   the code size. We inline the content manually here to avoid function
   call overhead. */
void
arb_clear(arb_t x)
{
    ARF_DEMOTE(arb_midref(x));

    if (COEFF_IS_MPZ(ARF_EXP(arb_midref(x))))
        _fmpz_clear_mpz(ARF_EXP(arb_midref(x)));

    if (COEFF_IS_MPZ(MAG_EXP(arb_radref(x))))
        _fmpz_clear_mpz(MAG_EXP(arb_radref(x)));
}

