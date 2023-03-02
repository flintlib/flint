/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"

void fmpz_fmma(fmpz_t f, const fmpz_t a, const fmpz_t b,
                         const fmpz_t c, const fmpz_t d)
{
    fmpz s, t, u, v;

    s = *a;
    t = *b;
    u = *c;
    v = *d;

    if (!COEFF_IS_MPZ(s) && !COEFF_IS_MPZ(t) &&
        !COEFF_IS_MPZ(u) && !COEFF_IS_MPZ(v))
    {
        mp_limb_t sh, sl, th, tl;

        smul_ppmm(sh, sl, s, t);
        smul_ppmm(th, tl, u, v);
        add_ssaaaa(sh, sl, sh, sl, th, tl);

        fmpz_set_signed_uiui(f, sh, sl);
        return;
    }

    if (s == 0 || t == 0)
    {
        fmpz_mul(f, c, d);
        return;
    }

    if (u == 0 || v == 0)
    {
        fmpz_mul(f, a, b);
        return;
    }

    if (f == c || f == d)
    {
        if (f == a || f == b)
        {
            fmpz_t t;
            fmpz_init(t);
            fmpz_mul(t, a, b);
            fmpz_addmul(t, c, d);
            fmpz_swap(t, f);
            fmpz_clear(t);
        }
        else
        {
            fmpz_mul(f, c, d);
            fmpz_addmul(f, a, b);
        }
    }
    else
    {
        fmpz_mul(f, a, b);
        fmpz_addmul(f, c, d);
    }
}

