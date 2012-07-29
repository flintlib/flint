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

    Copyright (C) 2011 Sebastian Pancratz 
    Copyright (C) 2012 Andres Goens

******************************************************************************/

#include "fq.h"

void
fq_randtest(fq_t x, flint_rand_t state, const fq_ctx_t ctx)
{
    padic_poly_randtest_val(x, state, 0, fq_ctx_dim(ctx), &ctx->pctx);
}

void
fq_randtest_not_zero(fq_t f, flint_rand_t state,
                             const fq_ctx_t ctx)
{
    long i;

    if (fq_ctx_dim(ctx) == 0)
    {
        printf("Exception (fq_randtest_not_zero).  dim == 0.\n");
        abort();
    }

    fq_randtest(f, state, ctx);
    for (i = 0; !fq_is_zero(f) && (i < 10); i++)
        fq_randtest(f, state, ctx);

    if (fq_is_zero(f))
    {
        padic_poly_fit_length(f, 1);
        _padic_poly_set_length(f, 1);
        fmpz_set_ui(f->coeffs + 0, 1);
        f->val = 0;
    }
}

