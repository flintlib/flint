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
nmod_poly_realloc(nmod_poly_t poly, slong alloc)
{
    if (alloc == 0)
    {
        nmod_poly_clear(poly);
        poly->length = 0;
        poly->alloc = 0;
        poly->coeffs = NULL;

        return;
    }

    poly->coeffs = (mp_ptr) flint_realloc(poly->coeffs, alloc * sizeof(mp_limb_t));

    poly->alloc = alloc;

    /* truncate poly if necessary */
    if (poly->length > alloc)
    {
        poly->length = alloc;
        _nmod_poly_normalise(poly);
    }
}
