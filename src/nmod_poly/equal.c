/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

int nmod_poly_equal(const nmod_poly_t a, const nmod_poly_t b)
{
    if (a->length != b->length)
        return 0;

    if (a != b)
        if (!_nmod_vec_equal(a->coeffs, b->coeffs, a->length))
            return 0;

    return 1;
}

int nmod_poly_equal_nmod(const nmod_poly_t poly, ulong cst)
{
    if (cst == 0)
        return nmod_poly_is_zero(poly);
    else
        return (poly->length == 1 && poly->coeffs[0] == cst);
}

int nmod_poly_equal_ui(const nmod_poly_t poly, ulong cst)
{
    return nmod_poly_equal_nmod(poly, nmod_set_ui(cst, poly->mod));
}
