/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"
#include "fq_zech.h"

void
fq_zech_mul(fq_zech_t rop, const fq_zech_t op1, const fq_zech_t op2,
            const fq_zech_ctx_t ctx)
{
    if (op1->value == ctx->qm1 || op2->value == ctx->qm1)
    {
        rop->value = ctx->qm1;
        return;
    }
    rop->value = n_addmod(op1->value, op2->value, ctx->qm1);
}

void
fq_zech_mul_fmpz(fq_zech_t rop, const fq_zech_t op, const fmpz_t x,
                 const fq_zech_ctx_t ctx)
{
    ulong ux;
    fmpz_t y;

    fmpz_init(y);
    fmpz_mod_ui(y, x, ctx->p);

    ux = fmpz_get_ui(y);

    fmpz_clear(y);

    fq_zech_mul_ui(rop, op, ux, ctx);
}

void
fq_zech_mul_si(fq_zech_t rop, const fq_zech_t op, slong x,
               const fq_zech_ctx_t ctx)
{
    ulong y;
    if (x == 0 || fq_zech_is_zero(op, ctx))
    {
        fq_zech_zero(rop, ctx);
        return;
    }
    if (x < 0)
    {
        y = -x;
        fq_zech_mul_ui(rop, op, y, ctx);
        fq_zech_neg(rop, rop, ctx);
    }
    else
    {
        y = x;
        fq_zech_mul_ui(rop, op, y, ctx);
    }
}

void
fq_zech_mul_ui(fq_zech_t rop, const fq_zech_t op, ulong x,
               const fq_zech_ctx_t ctx)
{
    ulong b;

    if (x == 0 || fq_zech_is_zero(op, ctx))
    {
        fq_zech_zero(rop, ctx);
        return;
    }

    b = x;
    if (x >= ctx->p)
        b = n_mod2_precomp(x, ctx->p, ctx->ppre);

    if (b == 0)
        fq_zech_zero(rop, ctx);
    else
        rop->value = n_addmod(op->value, ctx->prime_field_table[b], ctx->qm1);
}
