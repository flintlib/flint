/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "n_poly.h"

void
n_poly_mod_product_roots_nmod_vec(n_poly_t A, ulong_srcptr r, slong n, nmod_t mod)
{
    n_poly_fit_length(A, n + 1);
    A->length = n + 1;
    _nmod_poly_product_roots_nmod_vec(A->coeffs, r, n, mod);
}
