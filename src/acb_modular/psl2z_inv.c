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
psl2z_inv(psl2z_t h, const psl2z_t g)
{
    if (h != g)
        psl2z_set(h, g);

    if (fmpz_is_zero(&h->c) && fmpz_sgn(&h->a) > 0)
    {
        fmpz_neg(&h->b, &h->b);
        fmpz_swap(&h->d, &h->a);
    }
    else
    {
        fmpz_swap(&h->a, &h->d);
        fmpz_neg(&h->a, &h->a);
        fmpz_neg(&h->d, &h->d);
    }
}

