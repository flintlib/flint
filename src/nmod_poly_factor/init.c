/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_poly.h"

void
nmod_poly_factor_init(nmod_poly_factor_t fac)
{
    slong i;

    fac->alloc = 5;
    fac->num   = 0;
    fac->p     = flint_malloc(sizeof(nmod_poly_struct) * 5);
    fac->exp   = flint_malloc(sizeof(slong) * 5);

    for (i = 0; i < 5; i++)
        nmod_poly_init_preinv(fac->p + i, 1, 0);
}
