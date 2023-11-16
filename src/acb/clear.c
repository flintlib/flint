/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

/* The clear methods are not defined as inlines since this inflates
   the code size. We inline the content manually here to avoid function
   call overhead. Question: would calling arb_clear twice here be worth
   the call overhead, by improving branch prediction? */
void
acb_clear(acb_t x)
{
    ARF_DEMOTE(arb_midref(acb_realref(x)));

    if (COEFF_IS_MPZ(ARF_EXP(arb_midref(acb_realref(x)))))
        _fmpz_clear_mpz(ARF_EXP(arb_midref(acb_realref(x))));

    if (COEFF_IS_MPZ(MAG_EXP(arb_radref(acb_realref(x)))))
        _fmpz_clear_mpz(MAG_EXP(arb_radref(acb_realref(x))));

    ARF_DEMOTE(arb_midref(acb_imagref(x)));

    if (COEFF_IS_MPZ(ARF_EXP(arb_midref(acb_imagref(x)))))
        _fmpz_clear_mpz(ARF_EXP(arb_midref(acb_imagref(x))));

    if (COEFF_IS_MPZ(MAG_EXP(arb_radref(acb_imagref(x)))))
        _fmpz_clear_mpz(MAG_EXP(arb_radref(acb_imagref(x))));
}

