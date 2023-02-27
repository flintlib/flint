/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2014 Abhinav Baid
    Copyright (C) 2021 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"


/* r = a mod m  with various options for the range of r */
void _fmpz_smod(
    fmpz_t r,
    const fmpz_t a,
    const fmpz_t m,
    int sign, /* -1: |r| < |m| & sgn(r) = sgn(a) or r == 0
                  0: |r| < |m| & sgn(r) = sgn(m) or r == 0
                  1: -|m| < 2r <= |m|                */
    fmpz_t t) /* temp not aliased with anything else */
{
    if (sign < 0)
    {
        if (fmpz_cmpabs(m, a) > 0)
            fmpz_set(r, a);
        else
            fmpz_tdiv_qr(t, r, a, m);
    }
    else if (sign > 0)
    {
        int cmp = fmpz_cmp2abs(m, a);

        if (cmp >= 0)
        {
            if (cmp == 0)
                fmpz_abs(r, a);
            else
                fmpz_set(r, a);
        }
        else if (m != r)
        {
            fmpz_fdiv_qr(t, r, a, m);  /* r is zero or has same sign as m */

            cmp = fmpz_cmp2abs(m, r);

            if (cmp == 0)
                fmpz_abs(r, r);
            else if (cmp < 0)
                fmpz_sub(r, r, m);
        }
        else
        {
            fmpz_set(t, m);

            fmpz_fdiv_r(r, a, t);

            cmp = fmpz_cmp2abs(t, r);

            if (cmp == 0)
                fmpz_abs(r, r);
            else if (cmp < 0)
                fmpz_sub(r, r, t);
        }
    }
    else
    {
        fmpz_fdiv_qr(t, r, a, m);
    }
}


void fmpz_smod(fmpz_t f, const fmpz_t g, const fmpz_t h)
{
    fmpz c2 = *h;

    if (!COEFF_IS_MPZ(c2))      /* h is small */
    {
        ulong tmp;

        tmp = FLINT_ABS(c2);

        fmpz_mod(f, g, h);
        if (fmpz_cmp_ui(f, tmp / 2) > 0)
        {
            fmpz_sub_ui(f, f, tmp);
        }
    }
    else                        /* h is large */
    {
        fmpz_t tmp;
        fmpz_init(tmp);
        _fmpz_smod(f, g, h, 1, tmp);
        fmpz_clear(tmp);
    }
}
