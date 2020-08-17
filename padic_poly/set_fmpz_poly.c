/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "padic_poly.h"

void padic_poly_set_fmpz_poly(padic_poly_t f, const fmpz_poly_t g, 
                              const padic_ctx_t ctx)
{
    const slong len = g->length;

    padic_poly_fit_length(f, len);
    _padic_poly_set_length(f, len);
    _fmpz_vec_set(f->coeffs, g->coeffs, len);
    f->val = 0;

    padic_poly_canonicalise(f, ctx->p);
    padic_poly_reduce(f, ctx);
}

