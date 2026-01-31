/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmpcompat.h"
#include "fmpz.h"

void
fmpz_mul_tdiv_q_2exp(fmpz_t f, const fmpz_t g, const fmpz_t h, ulong exp)
{
    fmpz c1, c2;

    c1 = *g;

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
        fmpz_mul_si_tdiv_q_2exp(f, h, c1, exp);
        return;
    }

    c2 = *h;                    /* save h in case it is aliased with f */

    if (c2 == WORD(0))               /* special case, h = 0  */
    {
        fmpz_zero(f);
        return;
    }

    fmpz_mul(f, g, h);
    fmpz_tdiv_q_2exp(f, f, exp);
}
