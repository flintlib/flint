/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca.h"
#include "ca_ext.h"
#include "ca_field.h"

void
ca_ctx_init(ca_ctx_t ctx)
{
    qqbar_t onei;
    ca_ext_t ext;
    ca_ext_struct * ext_ptr[1];

    ctx->options = flint_calloc(CA_OPT_NUM_OPTIONS, sizeof(slong));

    ctx->options[CA_OPT_PREC_LIMIT] = 4096;
    ctx->options[CA_OPT_QQBAR_DEG_LIMIT] = 120;
    ctx->options[CA_OPT_LOW_PREC] = 64;
    ctx->options[CA_OPT_SMOOTH_LIMIT] = 32;
    ctx->options[CA_OPT_LLL_PREC] = 128;
    ctx->options[CA_OPT_POW_LIMIT] = 20;
    ctx->options[CA_OPT_USE_GROEBNER] = 1;
    ctx->options[CA_OPT_GROEBNER_LENGTH_LIMIT] = 100;
    ctx->options[CA_OPT_GROEBNER_POLY_LENGTH_LIMIT] = 1000;
    ctx->options[CA_OPT_GROEBNER_POLY_BITS_LIMIT] = 10000;
    ctx->options[CA_OPT_VIETA_LIMIT] = 6;
    ctx->options[CA_OPT_PRINT_FLAGS] = CA_PRINT_DEFAULT;
    ctx->options[CA_OPT_MPOLY_ORD] = ORD_LEX;
    ctx->options[CA_OPT_TRIG_FORM] = CA_TRIG_EXPONENTIAL;

    ctx->mctx = NULL;
    ctx->mctx_len = 0;

    ca_ext_cache_init(CA_CTX_EXT_CACHE(ctx), ctx);
    ca_field_cache_init(CA_CTX_FIELD_CACHE(ctx), ctx);

    /* Always create QQ */
    ctx->field_qq = ca_field_cache_insert_ext(CA_CTX_FIELD_CACHE(ctx), NULL, 0, ctx);

    /* Always create QQ(i) */
    qqbar_init(onei);
    qqbar_i(onei);
    ca_ext_init_qqbar(ext, onei, ctx);
    ext_ptr[0] = ca_ext_cache_insert(CA_CTX_EXT_CACHE(ctx), ext, ctx);
    ctx->field_qq_i = ca_field_cache_insert_ext(CA_CTX_FIELD_CACHE(ctx), ext_ptr, 1, ctx);
    qqbar_clear(onei);
    ca_ext_clear(ext, ctx);
}
