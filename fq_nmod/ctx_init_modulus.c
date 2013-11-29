/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2013 Mike Hansen

******************************************************************************/

#include <stdio.h>
#include <string.h>

#include "fq_nmod.h"

void
fq_nmod_ctx_init_modulus(fq_nmod_ctx_t ctx, const fmpz_t p, slong d,
                         const nmod_poly_t modulus, const char *var)
{
    slong nz;
    int i, j;
    mp_limb_t inv;

    fmpz_init_set(fq_nmod_ctx_prime(ctx), p);
    ctx->mod.n = fmpz_get_ui(p);
    ctx->mod.ninv = n_preinvert_limb(ctx->mod.n);
    count_leading_zeros(ctx->mod.norm, ctx->mod.n);

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
    nmod_poly_inv_series_newton(ctx->inv, ctx->inv, ctx->inv->length);
}
