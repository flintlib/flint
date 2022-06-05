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

#include "nmod_poly_factor-impl.h"


void
nmod_poly_factor_clear(nmod_poly_factor_t fac)
{
    slong i;

    for (i = 0; i < fac->alloc; i++)
        nmod_poly_clear(fac->p + i);

    flint_free(fac->p);
    flint_free(fac->exp);
}
