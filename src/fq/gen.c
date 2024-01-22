/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_poly.h"
#include "fq.h"

void fq_gen(fq_t rop, const fq_ctx_t ctx)
{
    if (ctx->modulus->length == 2)
    {
        fmpz_poly_fit_length(rop, 1);
        fmpz_invmod(rop->coeffs, ctx->modulus->coeffs + 1, fq_ctx_prime(ctx));
        fmpz_inplace_neg(rop->coeffs);
        fmpz_mul(rop->coeffs, rop->coeffs, ctx->modulus->coeffs);
        fmpz_mod(rop->coeffs, rop->coeffs, fq_ctx_prime(ctx));
        _fmpz_poly_set_length(rop, !fmpz_is_zero(rop->coeffs));
    }
    else
    {
        fmpz_poly_zero(rop);
        fmpz_poly_set_coeff_ui(rop, 0, 0);
        fmpz_poly_set_coeff_ui(rop, 1, 1);
    }
}
