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
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz_mod_poly.h"

void
fmpz_mod_poly_factor_init(fmpz_mod_poly_factor_t fac)
{
    slong i;
    fmpz_t p;

    fac->alloc = 5;
    fac->num = 0;
    fac->poly = flint_malloc(sizeof(fmpz_mod_poly_struct) * 5);
    fac->exp = flint_malloc(sizeof(slong) * 5);

    fmpz_init_set_ui(p, 5);
    for (i = 0; i < 5; i++)
        fmpz_mod_poly_init(fac->poly + i, p);
    fmpz_clear(p);
}
