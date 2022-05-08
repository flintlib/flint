/*
    Copyright (C) 2007, David Howden.
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"

void
nmod_poly_init2_preinv(nmod_poly_t poly,
                       ulong n, ulong ninv, slong alloc)
{
    if (alloc)
        poly->coeffs = (ulong_ptr) flint_malloc(alloc * sizeof(ulong));
    else
        poly->coeffs = NULL;

    poly->mod.n = n;
    poly->mod.ninv = ninv;

    count_leading_zeros(poly->mod.norm, n);

    poly->alloc = alloc;
    poly->length = 0;
}

void
nmod_poly_init2(nmod_poly_t poly, ulong n, slong alloc)
{
    nmod_poly_init2_preinv(poly, n, n_preinvert_limb(n), alloc);
}
