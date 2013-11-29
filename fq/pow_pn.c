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

#include "fq.h"

void
fq_pow_pn_init_precomp_matrix(fmpz_mat_t A, const fq_ctx_t ctx)
{
    fmpz_mat_init(A, n_sqrt(ctx->modulus->length - 1) + 1,
                  ctx->modulus->length - 1);
}

void
fq_pow_pn_clear_precomp_matrix(fmpz_mat_t A, const fq_ctx_t ctx)
{
    fmpz_mat_clear(A);
}

void
fq_pow_pn_precompute_matrix_fq(fmpz_mat_t A, fq_t op, const fq_ctx_t ctx)
{
    fmpz_poly_fit_length(op, fq_ctx_degree(ctx));
    _fmpz_mod_poly_precompute_matrix(A, op->coeffs,
                                     ctx->modulus->coeffs, ctx->modulus->length,
                                     ctx->inv->coeffs, ctx->inv->length,
                                     fq_ctx_prime(ctx));
}

void
fq_pow_pn_precompute_matrix_ui(fmpz_mat_t A, ulong n, const fq_ctx_t ctx)
{
    slong j;
    fq_t xp, xpi;
    fq_init(xp, ctx);
    fq_gen(xp, ctx);
    fq_pow(xp, xp, fq_ctx_prime(ctx), ctx);
    fq_pow_pn_precompute_matrix_fq(A, xp, ctx);
    if (n == 1)
    {
        fq_clear(xp, ctx);
        return;
    }
        
    fq_init(xpi, ctx);
    fq_set(xpi, xp, ctx);
    for (j = 1; j < n; j++)
        fq_pow_pn_precomp(xpi, xpi, A, ctx);

    fq_pow_pn_precompute_matrix_fq(A, xpi, ctx);
    fq_clear(xpi, ctx);
    fq_clear(xp, ctx);
}

void
fq_pow_pn_precomp(fq_t rop, const fq_t op1, const fmpz_mat_t A,
                  const fq_ctx_t ctx)
{
    slong len1 = op1->length;
    slong len3 = ctx->modulus->length;

    if (len1 == 0 || len3 == 1)
    {
        fq_zero(rop, ctx);
        return;
    }

    if (len1 == 1)
    {
        fq_set(rop, op1, ctx);
        return;
    }

    /* Handle aliasing */
    if (rop == op1)
    {
        fq_t tmp;
        fq_init(tmp, ctx);
        fq_pow_pn_precomp(tmp, op1, A, ctx);
        fq_swap(tmp, rop, ctx);
        fq_clear(tmp, ctx);
        return;
    }

    fmpz_poly_fit_length(rop, fq_ctx_degree(ctx));
    _fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(rop->coeffs,
                                                         op1->coeffs, op1->length, A,
                                                         ctx->modulus->coeffs, ctx->modulus->length,
                                                         ctx->inv->coeffs, ctx->inv->length,
                                                         fq_ctx_prime(ctx));
    _fmpz_poly_set_length(rop, fq_ctx_degree(ctx));
    _fmpz_poly_normalise(rop);
}
