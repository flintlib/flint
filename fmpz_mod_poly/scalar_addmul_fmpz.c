/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_poly.h"

void fmpz_mod_poly_scalar_addmul_fmpz(fmpz_mod_poly_t A,
             const fmpz_mod_poly_t B, const fmpz_t x, const fmpz_mod_ctx_t ctx)
{
    slong len = FLINT_MAX(A->length, B->length);

    if (fmpz_is_zero(x) || B->length < 1)
        return;

    fmpz_mod_poly_fit_length(A, B->length, ctx);
    
    if (B->length > A->length)
        _fmpz_vec_zero(A->coeffs + A->length, B->length - A->length);

    _fmpz_vec_scalar_mod_fmpz(A->coeffs, A->coeffs, len,
                                                    fmpz_mod_ctx_modulus(ctx));
    _fmpz_mod_poly_set_length(A, len);
    _fmpz_mod_poly_normalise(A);
}

