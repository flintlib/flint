/*
    Copyright (C) 2007, David Howden.
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

void
nmod_poly_neg(nmod_poly_t res, const nmod_poly_t poly1)
{
    nmod_poly_fit_length(res, poly1->length);

    _nmod_vec_neg(res->coeffs, poly1->coeffs, poly1->length, poly1->mod);

    res->length = poly1->length;
}
