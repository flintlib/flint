/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_modular.h"

static int
fmpz_kronecker1(const fmpz_t a, const fmpz_t b)
{
    if (fmpz_sgn(b) < 0)
    {
        int r;
        fmpz_t t;
        fmpz_init(t);
        fmpz_neg(t, b);
        r = fmpz_kronecker1(a, t);
        fmpz_clear(t);
        return r;
    }
    else if (fmpz_is_one(b))
    {
        return 1;
    }
    else
    {
        return fmpz_jacobi(a, b);
    }
}

int
acb_modular_epsilon_arg(const psl2z_t g)
{
#define a (&g->a)
#define b (&g->b)
#define c (&g->c)
#define d (&g->d)

    if (fmpz_is_zero(c))
    {
        return fmpz_fdiv_ui(b, 24);
    }
    else
    {
        int aa, bb, cc, dd;
        int u;

        aa = fmpz_fdiv_ui(a, 24);
        bb = fmpz_fdiv_ui(b, 24);
        cc = fmpz_fdiv_ui(c, 24);
        dd = fmpz_fdiv_ui(d, 24);

        if (cc % 2 == 1)
        {
            u = fmpz_kronecker1(a, c);
            aa = aa*bb + 2*aa*cc - 3*cc + cc*dd*(1-aa*aa);
        }
        else
        {
            u = fmpz_kronecker1(c, a);
            aa = aa*bb - aa*cc + 3*aa - 3 + cc*dd*(1-aa*aa);
        }

        if (u == -1)
        {
            aa += 12;
        }
        else if (u != 1)
        {
            flint_throw(FLINT_ERROR, "bad kronecker input\n");
        }

        /* mod 24 */
        if (aa < 0)
        {
            aa = 24 - ((-aa) % 24);
            if (aa == 24)
                aa = 0;
        }
        else
            aa = aa % 24;

        return aa;
    }

#undef a
#undef b
#undef c
#undef d
}

