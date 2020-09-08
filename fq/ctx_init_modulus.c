/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <string.h>

#include "fq.h"
#include "fq_poly.h"

void fq_ctx_init_modulus(fq_ctx_t ctx, const fmpz_mod_poly_t modulus,
                                    const fmpz_mod_ctx_t ctxp, const char *var)
{
    slong nz;
    int i, j;
    const fmpz * p = fmpz_mod_ctx_modulus(ctxp);
    fmpz_t inv;

    /* Count number of nonzero coefficients */
    nz = 0;
    for (i = 0; i < modulus->length; i++)
        nz += !fmpz_is_zero(modulus->coeffs + i);

    ctx->len = nz;
    ctx->a = _fmpz_vec_init(ctx->len);
    ctx->j = flint_malloc(ctx->len * sizeof(slong));

    fmpz_init(inv);
    fmpz_invmod(inv, modulus->coeffs + modulus->length - 1, p);

    /* Copy the polynomial */
    j = 0;
    for (i = 0; i < modulus->length; i++)
    {
        if (!fmpz_is_zero(modulus->coeffs + i))
        {
            fmpz_mul(ctx->a + j, inv, modulus->coeffs + i);
            fmpz_mod(ctx->a + j, ctx->a + j, p);
            ctx->j[j] = i;
            j++;
        }
    }

    fmpz_clear(inv);

    if (ctx->len < 6)
        ctx->sparse_modulus = 1;
    else
        ctx->sparse_modulus = 0;

    fmpz_mod_ctx_init(ctx->ctxp, p);

    ctx->var = flint_malloc(strlen(var) + 1);
    strcpy(ctx->var, var);

    /* Set the modulus */
    fmpz_mod_poly_init(ctx->modulus, ctx->ctxp);
    fmpz_mod_poly_set(ctx->modulus, modulus, ctx->ctxp);

    /* Precompute the inverse of the modulus */
    fmpz_mod_poly_init(ctx->inv, ctx->ctxp);
    fmpz_mod_poly_reverse(ctx->inv, ctx->modulus, ctx->modulus->length, ctx->ctxp);
    fmpz_mod_poly_inv_series_newton(ctx->inv, ctx->inv, ctx->modulus->length, ctx->ctxp);

    ctx->is_conway = 0;
}
