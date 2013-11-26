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

#include "fq_nmod.h"

void
fq_nmod_pow_pn_init_precomp_matrix(nmod_mat_t A, const fq_nmod_ctx_t ctx)
{
    nmod_mat_init(A, n_sqrt(ctx->modulus->length - 1) + 1,
                  ctx->modulus->length - 1, ctx->mod.n);
}

void
fq_nmod_pow_pn_clear_precomp_matrix(nmod_mat_t A, const fq_nmod_ctx_t ctx)
{
    nmod_mat_clear(A);
}

void
fq_nmod_pow_pn_precompute_matrix_fq_nmod(nmod_mat_t A, const fq_nmod_t op,
                                         const fq_nmod_ctx_t ctx)
{
    nmod_poly_precompute_matrix(A, op, ctx->modulus, ctx->inv);
}

void
fq_nmod_pow_pn_precompute_matrix_ui(nmod_mat_t A, ulong n,
                                    const fq_nmod_ctx_t ctx)
{
    slong j;
    fq_nmod_t xp, xpi;
    fq_nmod_init(xp, ctx);
    fq_nmod_gen(xp, ctx);
    fq_nmod_pow(xp, xp, fq_nmod_ctx_prime(ctx), ctx);
    fq_nmod_pow_pn_precompute_matrix_fq_nmod(A, xp, ctx);
    if (n == 1)
    {
        fq_nmod_clear(xp, ctx);
        return;
    }
        
    fq_nmod_init(xpi, ctx);
    fq_nmod_set(xpi, xp, ctx);
    for (j = 1; j < n; j++)
        fq_nmod_pow_pn_precomp(xpi, xpi, A, ctx);

    fq_nmod_pow_pn_precompute_matrix_fq_nmod(A, xpi, ctx);
    fq_nmod_clear(xpi, ctx);
    fq_nmod_clear(xp, ctx);
}

void
fq_nmod_pow_pn_precomp(fq_nmod_t rop, const fq_nmod_t op1, const nmod_mat_t A,
                       const fq_nmod_ctx_t ctx)
{
    nmod_poly_compose_mod_brent_kung_precomp_preinv(rop, op1, A, ctx->modulus,
                                                    ctx->inv);
}

void
fq_nmod_pow_pn(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_t op2,
                       const fq_nmod_ctx_t ctx)
{
    nmod_poly_compose_mod_brent_kung_preinv(rop, op1, op2, ctx->modulus, ctx->inv);
}
