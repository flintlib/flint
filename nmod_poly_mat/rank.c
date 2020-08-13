/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "nmod_poly.h"
#include "nmod_poly_mat.h"


slong
nmod_poly_mat_rank(const nmod_poly_mat_t A)
{
    nmod_poly_mat_t tmp;
    nmod_poly_t den;
    slong rank;

    if (nmod_poly_mat_is_empty(A))
        return 0;

    nmod_poly_mat_init_set(tmp, A);
    nmod_poly_init(den, nmod_poly_mat_modulus(A));
    rank = nmod_poly_mat_fflu(tmp, den, NULL, tmp, 0);
    nmod_poly_mat_clear(tmp);
    nmod_poly_clear(den);
    return rank;
}
