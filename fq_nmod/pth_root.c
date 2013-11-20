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
fq_nmod_pth_root(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_ctx_t ctx)
{
    slong i, d;
    if (fq_nmod_is_zero(op1, ctx) || fq_nmod_is_one(op1, ctx))
    {
        fq_nmod_set(rop, op1, ctx);
        return;
    }

    /* TODO: Use modular composition */
    d = fq_nmod_ctx_degree(ctx) - 1;
    fq_nmod_set(rop, op1, ctx);
    for (i = 0; i < d; i++)
    {
        fq_nmod_pow(rop, rop, fq_nmod_ctx_prime(ctx), ctx);
    }
}
