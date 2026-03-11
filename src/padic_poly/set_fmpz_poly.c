/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "padic_poly.h"

void padic_poly_set_fmpz_poly(padic_poly_t rop, const fmpz_poly_t op,
                              const padic_ctx_t ctx)
{
    const slong len = op->length;

    padic_poly_fit_length(rop, len);
    _padic_poly_set_length(rop, len);
    _fmpz_vec_set(rop->coeffs, op->coeffs, len);
    rop->val = 0;

    padic_poly_canonicalise(rop, ctx->p);
    padic_poly_reduce(rop, ctx);
}
