/*
    Copyright (C) 2012 William Hart

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

void
fmpz_xgcd(fmpz_t d, fmpz_t a, fmpz_t b, const fmpz_t f, const fmpz_t g)
{
    fmpz_t t1, t2;
    fmpz *f1, *g1;

    if (fmpz_is_zero(f))
    {
        int sign = fmpz_sgn(g);
        fmpz_abs(d, g);
        fmpz_set_ui(a, 0);
        if (sign == 0)
            fmpz_set_ui(b, 0);
        else if (sign > 0)
            fmpz_set_ui(b, 1);
        else
            fmpz_set_si(b, -1);
    }
    else if (fmpz_cmpabs(f, g) == 0)
    {
        if (fmpz_sgn(f) > 0)
        {
            fmpz_set(d, f);
            fmpz_set_ui(a, 1);
        }
        else
        {
            fmpz_neg(d, f);
            fmpz_set_si(a, -1);
        }
        fmpz_set_si(b, 0);
    }
    else
    {
        int sign1 = fmpz_sgn(f);
        int sign2 = fmpz_sgn(g);

        fmpz_init(t1);
        fmpz_init(t2);

        /* support aliasing */
        if (d == f || a == f || sign1 < 0)
        {
            f1 = t1;
            if (sign1 < 0)
                fmpz_neg(f1, f);
            else
                fmpz_set(f1, f);
        }
        else
            f1 = (fmpz *) f;

        if (d == g || a == g || sign2 < 0)
        {
            g1 = t2;
            if (sign2 < 0)
                fmpz_neg(g1, g);
            else
                fmpz_set(g1, g);
        }
        else
            g1 = (fmpz *) g;

        if (fmpz_cmp(f1, g1) < 0)
        {
            fmpz_gcdinv(d, a, f1, g1);
            fmpz_mul(t1, a, f1);
            fmpz_sub(t1, d, t1);
            fmpz_divexact(b, t1, g1);
        }
        else                    /* g < f */
        {
            fmpz_gcdinv(d, b, g1, f1);
            fmpz_mul(t2, b, g1);
            fmpz_sub(t2, d, t2);
            fmpz_divexact(a, t2, f1);
        }

        if (sign1 < 0)
            fmpz_neg(a, a);
        if (sign2 < 0)
            fmpz_neg(b, b);

        fmpz_clear(t1);
        fmpz_clear(t2);
    }
}
