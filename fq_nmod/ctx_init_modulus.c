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

#include "fq_nmod.h"

void
fq_nmod_ctx_init_modulus(fq_nmod_ctx_t ctx, const nmod_poly_t modulus,
                         const char *var)
{
    slong nz;
    int i, j;
    mp_limb_t inv;

    fmpz_init(fq_nmod_ctx_prime(ctx));
    fmpz_set_ui(fq_nmod_ctx_prime(ctx), modulus->mod.n);
    
    ctx->mod.n = modulus->mod.n;
    ctx->mod.ninv = modulus->mod.ninv;
    ctx->mod.norm = modulus->mod.norm;

    /* Count number of nonzero coefficients */
    nz = 0;
    for (i = 0; i < modulus->length; i++)
    {
        if (modulus->coeffs[i] != 0)
        {
            nz += 1;
        }
    }

    ctx->len = nz;
    ctx->a = _nmod_vec_init(ctx->len);
    ctx->j = flint_malloc(ctx->len * sizeof(mp_limb_t));

    inv = n_invmod(modulus->coeffs[modulus->length - 1], ctx->mod.n);

    /* Copy the polynomial */
    j = 0;
    for (i = 0; i < modulus->length; i++)
    {
        if (modulus->coeffs[i] != 0)
        {
            ctx->a[j] = n_mulmod2_preinv(inv, modulus->coeffs[i],
                                         ctx->mod.n, ctx->mod.ninv);
            ctx->j[j] = i;
            j++;
        }
    }

    if (ctx->len < 6)
        ctx->sparse_modulus = 1;
    else
        ctx->sparse_modulus = 0;

    ctx->var = flint_malloc(strlen(var) + 1);
    strcpy(ctx->var, var);

    /* Set the modulus */
    nmod_poly_init(ctx->modulus, ctx->mod.n);
    nmod_poly_set(ctx->modulus, modulus);

    /* Precompute the inverse of the modulus */
    nmod_poly_init(ctx->inv, ctx->mod.n);
    nmod_poly_reverse(ctx->inv, ctx->modulus, ctx->modulus->length);
    nmod_poly_inv_series_newton(ctx->inv, ctx->inv, ctx->modulus->length);

    ctx->is_conway = 0;
}

