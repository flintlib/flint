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

#include "fq_zech.h"

void
fq_zech_pow_pn_init_precomp_matrix(nmod_mat_t A, const fq_zech_ctx_t ctx)
{
    nmod_mat_init(A, 1, 1, ctx->qm1);
}

void
fq_zech_pow_pn_clear_precomp_matrix(nmod_mat_t A, const fq_zech_ctx_t ctx)
{
    nmod_mat_clear(A);
}

void
fq_zech_pow_pn_precompute_matrix_fq_zech(nmod_mat_t A, const fq_zech_t op,
                                         const fq_zech_ctx_t ctx)
{
    nmod_mat_entry(A, 0, 0) = op->value;
}

void
fq_zech_pow_pn_precompute_matrix_ui(nmod_mat_t A, ulong n,
                                    const fq_zech_ctx_t ctx)
{
    slong j;
    fq_zech_t xp, xpi;
    fq_zech_init(xp, ctx);
    fq_zech_gen(xp, ctx);
    fq_zech_pow(xp, xp, fq_zech_ctx_prime(ctx), ctx);
    fq_zech_pow_pn_precompute_matrix_fq_zech(A, xp, ctx);
    if (n == 1)
    {
        fq_zech_clear(xp, ctx);
        return;
    }
        
    fq_zech_init(xpi, ctx);
    fq_zech_set(xpi, xp, ctx);
    for (j = 1; j < n; j++)
        fq_zech_pow_pn_precomp(xpi, xpi, A, ctx);

    fq_zech_pow_pn_precompute_matrix_fq_zech(A, xpi, ctx);
    fq_zech_clear(xpi, ctx);
    fq_zech_clear(xp, ctx);
}

void
fq_zech_pow_pn_precomp(fq_zech_t rop, const fq_zech_t op1, const nmod_mat_t A,
                       const fq_zech_ctx_t ctx)
{
    if (op1->value == ctx->qm1)
    {
        fq_zech_zero(rop, ctx);
        return;        
    }
    rop->value = n_mulmod_precomp(op1->value, nmod_mat_entry(A, 0, 0),
                                  ctx->qm1, ctx->qm1inv);
}
