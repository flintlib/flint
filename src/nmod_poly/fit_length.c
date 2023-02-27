/*
    Copyright (C) 2007, David Howden.
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

void
nmod_poly_fit_length(nmod_poly_t poly, slong alloc)
{
    if (alloc > poly->alloc)
    {
        if (alloc < 2 * poly->alloc)
            alloc = 2 * poly->alloc;

        nmod_poly_realloc(poly, alloc);
    }
}
