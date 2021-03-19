/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"

/* r = a mod m  with various options for the range of r */
void _fmpz_mods(
    fmpz_t r,
    const fmpz_t a,
    const fmpz_t m,
    int sign, /* -1: don't care just get close,  0: 0<=r<m,  1: 2|r|<=m  */
    fmpz_t t) /* temp */
{
    FLINT_ASSERT(fmpz_sgn(m) > 0);

    if (sign < 0)
    {
        if (fmpz_cmpabs(m, a) > 0)
            fmpz_set(r, a);
        else
            fmpz_tdiv_qr(t, r, a, m);
    }
    else if (sign > 0)
    {
        if (fmpz_cmp2abs(m, a) >= 0)
        {
            fmpz_set(r, a);
            return;
        }

        fmpz_tdiv_qr(t, r, a, m);

        if (fmpz_cmp2abs(m, r) >= 0)
            return;

        if (fmpz_sgn(r) < 0)
            fmpz_add(r, r, m);
        else
            fmpz_sub(r, r, m);
    }
    else
    {
        fmpz_fdiv_qr(t, r, a, m);
    }
}

