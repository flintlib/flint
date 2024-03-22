/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "nmod_poly.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"

void fmpz_mod_poly_get_fmpz_poly(fmpz_poly_t f, const fmpz_mod_poly_t g, const fmpz_mod_ctx_t FLINT_UNUSED(ctx))
{
    fmpz_poly_fit_length(f, g->length);
    _fmpz_poly_set_length(f, g->length);
    _fmpz_vec_set(f->coeffs, g->coeffs, g->length);
}

void fmpz_mod_poly_get_nmod_poly(nmod_poly_t f, const fmpz_mod_poly_t g)
{
    slong i;

    nmod_poly_fit_length(f, g->length);
    _nmod_poly_set_length(f, g->length);

    for (i = 0; i < g->length; i++)
       f->coeffs[i] = fmpz_get_ui(g->coeffs + i);
}
