/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_modular.h"

void
psl2z_randtest(psl2z_t g, flint_rand_t state, slong bits)
{
    bits = FLINT_MAX(bits, 1);

    fmpz_randtest(&g->a, state, bits);
    fmpz_randtest(&g->b, state, bits);

    if (fmpz_is_zero(&g->a) && fmpz_is_zero(&g->b))
    {
        psl2z_one(g);
    }
    else
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_xgcd(t, &g->d, &g->c, &g->a, &g->b);
        fmpz_divexact(&g->a, &g->a, t);
        fmpz_divexact(&g->b, &g->b, t);

        if (fmpz_sgn(&g->c) < 0)
            fmpz_neg(&g->c, &g->c);
        else
            fmpz_neg(&g->b, &g->b);

        if (fmpz_is_zero(&g->c) && fmpz_sgn(&g->d) < 0)
        {
            fmpz_neg(&g->a, &g->a);
            fmpz_neg(&g->b, &g->b);
            fmpz_neg(&g->c, &g->c);
            fmpz_neg(&g->d, &g->d);
        }

        fmpz_clear(t);
    }
}

