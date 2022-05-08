/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"
#include "fq.h"

void
_fq_dense_reduce(fmpz* R, slong lenR, const fq_ctx_t ctx)
{
    fmpz  *q, *r;

    if (lenR < ctx->modulus->length)
    {
        _fmpz_vec_scalar_mod_fmpz(R, R, lenR, fq_ctx_prime(ctx));
        return;
    }
    
    q = _fmpz_vec_init(lenR - ctx->modulus->length + 1);
    r = _fmpz_vec_init(ctx->modulus->length - 1);

    _fmpz_mod_poly_divrem_newton_n_preinv(q, r, R, lenR, 
                                        ctx->modulus->coeffs, ctx->modulus->length,
                                        ctx->inv->coeffs, ctx->inv->length,
                                        fq_ctx_prime(ctx));

    _fmpz_vec_set(R, r, ctx->modulus->length - 1);
    _fmpz_vec_clear(q, lenR - ctx->modulus->length + 1);
    _fmpz_vec_clear(r, ctx->modulus->length - 1);
}
