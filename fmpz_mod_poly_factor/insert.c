/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "fmpz_mod_poly.h"

void
fmpz_mod_poly_factor_insert(fmpz_mod_poly_factor_t fac,
                            const fmpz_mod_poly_t poly, slong exp)
{
    slong i;
    fmpz_t p;

    if (poly->length <= 1)
        return;

    for (i = 0; i < fac->num; i++)
    {
        if (fmpz_mod_poly_equal(poly, fac->poly + i))
        {
            fac->exp[i] += exp;
            return;
        }
    }

    if (fac->alloc == fac->num)
    {
        slong new_size = 2 * fac->alloc;

        fac->poly =
            flint_realloc(fac->poly, sizeof(fmpz_mod_poly_struct) * new_size);
        fac->exp = flint_realloc(fac->exp, sizeof(slong) * new_size);

        fmpz_init_set_ui(p, 5);
        for (i = fac->alloc; i < new_size; i++)
            fmpz_mod_poly_init(fac->poly + i, p);
        fmpz_clear(p);

        fac->alloc = new_size;
    }

    fmpz_mod_poly_set(fac->poly + fac->num, poly);
    fmpz_set(&(fac->poly + fac->num)->p, &poly->p);
    fac->exp[fac->num] = exp;
    fac->num++;
}
