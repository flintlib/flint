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
nmod_poly_init_preinv(nmod_poly_t poly, mp_limb_t n, mp_limb_t ninv)
{
    poly->coeffs = NULL;

    poly->alloc = 0;
    poly->length = 0;

    poly->mod.n = n;
    poly->mod.ninv = ninv;
    count_leading_zeros(poly->mod.norm, n);
}

void
nmod_poly_init(nmod_poly_t poly, mp_limb_t n)
{
    nmod_poly_init_preinv(poly, n, n_preinvert_limb(n));
}
