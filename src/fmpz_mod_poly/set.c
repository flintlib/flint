/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_vec.h"
#include "fmpz_mod_poly.h"

void fmpz_mod_poly_set(fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2,
                                                      const fmpz_mod_ctx_t ctx)
{
    if (poly1 != poly2)         /* Aliasing is trivial */
    {
        slong i, len = poly2->length;

        fmpz_mod_poly_fit_length(poly1, len, ctx);

        for (i = 0; i < len; i++)
            fmpz_set(poly1->coeffs + i, poly2->coeffs + i);

        _fmpz_mod_poly_set_length(poly1, len);
    }
}

void fmpz_mod_poly_set_fmpz(fmpz_mod_poly_t poly, const fmpz_t c,
                                                      const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_poly_fit_length(poly, 1, ctx);
    fmpz_mod(poly->coeffs, c, fmpz_mod_ctx_modulus(ctx));
    _fmpz_mod_poly_set_length(poly, 1);
    _fmpz_mod_poly_normalise(poly);
}

void fmpz_mod_poly_set_fmpz_poly(fmpz_mod_poly_t f, const fmpz_poly_t g,
                                                      const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_poly_fit_length(f, g->length, ctx);
    _fmpz_mod_poly_set_length(f, g->length);
    _fmpz_mod_vec_set_fmpz_vec(f->coeffs, g->coeffs, g->length, ctx);
    _fmpz_mod_poly_normalise(f);
}

void fmpz_mod_poly_set_nmod_poly(fmpz_mod_poly_t f, const nmod_poly_t g)
{
    slong i;

    _fmpz_mod_poly_fit_length(f, g->length);
    _fmpz_mod_poly_set_length(f, g->length);

    for (i = 0; i < g->length; i++)
       fmpz_set_ui(f->coeffs + i, g->coeffs[i]);
}

void fmpz_mod_poly_set_ui(fmpz_mod_poly_t f, ulong x, const fmpz_mod_ctx_t ctx)
{
    if (x == 0)
    {
        fmpz_mod_poly_zero(f, ctx);
    }
    else
    {
        fmpz_mod_poly_fit_length(f, 1, ctx);
        _fmpz_mod_poly_set_length(f, 1);
        fmpz_set_ui(f->coeffs, x);
        fmpz_mod(f->coeffs, f->coeffs, fmpz_mod_ctx_modulus(ctx));
        _fmpz_mod_poly_normalise(f);
    }
}
