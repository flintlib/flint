/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca.h"
#include "ca_ext.h"
#include "ca_field.h"

#define CA_NVARS_DEFAULT 8

void
ca_ctx_init(ca_ctx_t ctx)
{
    slong i;
    qqbar_t onei;

    ctx->mctx = flint_malloc(CA_NVARS_DEFAULT * sizeof(fmpz_mpoly_ctx_struct));
    for (i = 0; i < CA_NVARS_DEFAULT; i++)
        fmpz_mpoly_ctx_init(ctx->mctx + i, i + 1, CA_MPOLY_ORD);
    ctx->mctx_len = CA_NVARS_DEFAULT;

    ca_ext_cache_init(CA_CTX_EXT_CACHE(ctx), ctx);

    /* Always create QQ, QQ(i) */

    ctx->fields = (ca_field_struct *) flint_malloc(2 * sizeof(ca_field_struct));
    ctx->fields_len = 2;
    ctx->fields_alloc = 2;

    ca_field_init_qq(ctx->fields, ctx);

    qqbar_init(onei);
    qqbar_i(onei);
    ca_field_init_nf(ctx->fields + 1, onei, ctx);
    qqbar_clear(onei);

    ctx->options = flint_calloc(CA_OPT_NUM_OPTIONS, sizeof(slong));

    ctx->options[CA_OPT_PREC_LIMIT] = 4096;
    ctx->options[CA_OPT_QQBAR_DEG_LIMIT] = 120;
    ctx->options[CA_OPT_LOW_PREC] = 64;
    ctx->options[CA_OPT_SMOOTH_LIMIT] = 32;
}

