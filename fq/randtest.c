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

    Copyright (C) 2011 Sebastian Pancratz, 2012 Andres Goens

******************************************************************************/

#include "fq.h"

void
fq_randtest(fq_t x, flint_rand_t state, const fq_ctx_t ctx)
{
    padic_poly_randtest(x, state, fq_ctx_dim(ctx), &ctx->pctx);
}

void
fq_randtest_not_zero(fq_t x, flint_rand_t state, const fq_ctx_t ctx)
{
    padic_poly_randtest_not_zero(x, state, fq_ctx_dim(ctx), &ctx->pctx);
}

void
fq_randtest_val(fq_t x, flint_rand_t state, long val, const fq_ctx_t ctx)
{
    padic_poly_randtest_val(x, state, val, fq_ctx_dim(ctx), &ctx->pctx);
}

void
fq_randtest_int(fq_t x, flint_rand_t state, const fq_ctx_t ctx)
{
    const long N = (&ctx->pctx)->N;

    if (N <= 0)
    {
        padic_poly_zero(x);
    }
    else
    {
        padic_poly_randtest_val(x, state, n_randint(state, N),
                                fq_ctx_dim(ctx), &ctx->pctx);
    }
}
