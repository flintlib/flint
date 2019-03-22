/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <string.h>

#include "fq_nmod.h"

/*
void fq_nmod_ctx_init_modulusdr(fq_nmod_ctx_t ctx, const nmod_polydr_t modulus, mp_limb_t p,
                         const char *var)
{
    nmod_ctx_t fpctx;
    nmod_poly_t m;

    nmod_ctx_init(fpctx, p);
    m->mod = fpctx->mod;
    m->coeffs = modulus->coeffs;
    m->length = modulus->length;
    m->alloc = modulus->alloc;

    fq_nmod_ctx_init_modulus(ctx, m, var);

    nmod_ctx_clear(fpctx);
}
*/
void
fq_nmod_ctx_init_modulus(fq_nmod_ctx_t ctx, const nmod_poly_t modulus,
                         const char *var)
{
    slong nz;
    int i, j;
    mp_limb_t inv;

    fmpz_init_set_ui(fq_nmod_ctx_prime(ctx), modulus->mod.n);
    nmod_ctx_init_mod(ctx->fpctx, modulus->mod);

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

    inv = n_invmod(modulus->coeffs[modulus->length - 1], ctx->fpctx->mod.n);

    /* Copy the polynomial */
    j = 0;
    for (i = 0; i < modulus->length; i++)
    {
        if (modulus->coeffs[i] != 0)
        {
            ctx->a[j] = n_mulmod2_preinv(inv, modulus->coeffs[i],
                                      ctx->fpctx->mod.n, ctx->fpctx->mod.ninv);
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
    nmod_poly_init(ctx->modulus, modulus->mod.n);
    nmod_poly_set(ctx->modulus, modulus);

    /* Precompute the inverse of the modulus */
    nmod_polydr_init(ctx->inv, ctx->fpctx);
    nmod_polydr_reverse(ctx->inv, ctx->modulus, ctx->modulus->length, ctx->fpctx);
    nmod_polydr_inv_series_newton(ctx->inv, ctx->inv, ctx->modulus->length, ctx->fpctx);
}
