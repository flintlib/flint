/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"

int
arf_is_int(const arf_t x)
{
    mp_size_t xn;
    mp_srcptr xp;
    slong exp, c;

    exp = ARF_EXP(x);

    if (ARF_IS_SPECIAL(x))
        return exp == ARF_EXP_ZERO;

    if (COEFF_IS_MPZ(exp))
        return mpz_sgn(COEFF_TO_PTR(exp)) > 0;

    ARF_GET_MPN_READONLY(xp, xn, x);
    c = flint_ctz(xp[0]);
    return exp - xn * FLINT_BITS + c >= 0;
}

