/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca.h"

void
ca_field_init_qq(ca_field_t K)
{
    K->type = CA_FIELD_TYPE_QQ;
}

void
ca_field_init_nf(ca_field_t K, const qqbar_t x)
{
    fmpq_poly_t t;

    K->type = CA_FIELD_TYPE_NF;

    qqbar_init(&K->data.nf.x);
    qqbar_set(&K->data.nf.x, x);

    /* nf_init wants an fmpq_poly_t, so mock up one */
    t->coeffs = QQBAR_POLY(x)->coeffs;
    t->den[0] = 1;
    t->length = QQBAR_POLY(x)->length;
    t->alloc = QQBAR_POLY(x)->alloc;

    nf_init(&K->data.nf.nf, t);
}

void
ca_field_init_const(ca_field_t K, ulong func)
{
    K->type = CA_FIELD_TYPE_FUNC;
    K->data.func.func = func;
    K->data.func.args_len = 0;
    K->data.func.args = NULL;

    acb_init(&K->data.func.enclosure);
    acb_indeterminate(&K->data.func.enclosure);

/*
    if (func == CA_Pi)
        acb_const_pi(&K->data.func.enclosure, 128);
    else
        flint_abort();
*/
}

void ca_field_init_fx(ca_field_t K, ulong func, const ca_t x, ca_ctx_t ctx)
{
    K->type = CA_FIELD_TYPE_FUNC;
    K->data.func.func = func;
    K->data.func.args_len = 1;
    K->data.func.args = ca_vec_init(1, ctx);
    ca_set(K->data.func.args, x, ctx);

    acb_init(&K->data.func.enclosure);
    acb_indeterminate(&K->data.func.enclosure);
}

void ca_field_init_fxy(ca_field_t K, ulong func, const ca_t x, const ca_t y, ca_ctx_t ctx)
{
    K->type = CA_FIELD_TYPE_FUNC;
    K->data.func.func = func;
    K->data.func.args_len = 2;
    K->data.func.args = ca_vec_init(2, ctx);
    ca_set(K->data.func.args, x, ctx);
    ca_set(K->data.func.args + 1, y, ctx);

    acb_init(&K->data.func.enclosure);
    acb_indeterminate(&K->data.func.enclosure);
}

void
ca_field_init_multi(ca_field_t K, slong len, ca_ctx_t ctx)
{
    K->type = CA_FIELD_TYPE_MULTI;
    K->data.multi.len = len;
    K->data.multi.ext = flint_malloc(len * sizeof(slong));
    K->data.multi.ideal = NULL;
    K->data.multi.ideal_len = 0;

    while (ctx->mctx_len < len)
    {
        slong i;
        ctx->mctx = flint_realloc(ctx->mctx, 2 * ctx->mctx_len * sizeof(fmpz_mpoly_ctx_struct));
        for (i = ctx->mctx_len; i < 2 * ctx->mctx_len; i++)
            fmpz_mpoly_ctx_init(ctx->mctx + i, i + 1, ORD_LEX);
        ctx->mctx_len *= 2;
    }
}

