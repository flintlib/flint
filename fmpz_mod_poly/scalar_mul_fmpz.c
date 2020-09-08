/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"

void _fmpz_mod_poly_scalar_mul_fmpz(fmpz *res, const fmpz *poly, slong len, 
                                    const fmpz_t x, const fmpz_t p)
{
    _fmpz_vec_scalar_mul_fmpz(res, poly, len, x);
    _fmpz_vec_scalar_mod_fmpz(res, res, len, p);
}

void fmpz_mod_poly_scalar_mul_fmpz(fmpz_mod_poly_t res,
          const fmpz_mod_poly_t poly, const fmpz_t x, const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_poly_fit_length(res, poly->length, ctx);
    _fmpz_mod_poly_scalar_mul_fmpz(res->coeffs, poly->coeffs, poly->length,
                                                 x, fmpz_mod_ctx_modulus(ctx));

    _fmpz_mod_poly_set_length(res, poly->length);
    _fmpz_mod_poly_normalise(res);
}

