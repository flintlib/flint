/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly_mini.h"
#include "flint-impl.h"

void nmod_poly_set(nmod_poly_t a, const nmod_poly_t b)
{
    if (a != b)
    {
        nmod_poly_fit_length(a, b->length);
        FLINT_MPN_COPYI(a->coeffs, b->coeffs, b->length);
        a->length = b->length;
    }
}
