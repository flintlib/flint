/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "aprcl.h"

void
unity_zpq_clear(unity_zpq f)
{
    slong i;

    for (i = 0; i < f->p; i++)
    {
        fmpz_mod_poly_clear(f->polys[i]);
    }

    f->p = 0;
    f->q = 0;

    fmpz_clear(f->n);
    flint_free(f->polys);
}
