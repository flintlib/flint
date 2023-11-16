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
_ca_ctx_init_mctx(ca_ctx_t ctx, slong len)
{
    while (ctx->mctx_len < len)
    {
        slong i, alloc;
        alloc = FLINT_MAX(1, 2 * ctx->mctx_len);
        ctx->mctx = flint_realloc(ctx->mctx, alloc * sizeof(fmpz_mpoly_ctx_struct *));

        for (i = ctx->mctx_len; i < alloc; i++)
        {
            ctx->mctx[i] = flint_malloc(sizeof(fmpz_mpoly_ctx_struct));
            fmpz_mpoly_ctx_init(ctx->mctx[i], i + 1, ctx->options[CA_OPT_MPOLY_ORD]);
        }
        ctx->mctx_len = alloc;
    }
}


void
ca_field_init_qq(ca_field_t K, ca_ctx_t ctx)
{
    CA_FIELD_LENGTH(K) = 0;
    CA_FIELD_EXT(K) = NULL;
    CA_FIELD_IDEAL_P(K) = NULL;
    CA_FIELD_IDEAL_LENGTH(K) = 0;
    CA_FIELD_IDEAL_ALLOC(K) = 0;
    CA_FIELD_HASH(K) = 0;
}

void
ca_field_init_nf(ca_field_t K, const qqbar_t x, ca_ctx_t ctx)
{
    ca_ext_ptr ext;
    ca_ext_t tmp;

    /* todo: avoid unnecessary work (especially the nf_t init) */
    ca_ext_init_qqbar(tmp, x, ctx);
    ext = ca_ext_cache_insert(CA_CTX_EXT_CACHE(ctx), tmp, ctx);
    ca_ext_clear(tmp, ctx);

    CA_FIELD_LENGTH(K) = 1;
    CA_FIELD_EXT(K) = flint_malloc(sizeof(ca_ext_ptr));
    CA_FIELD_EXT_ELEM(K, 0) = ext;
    CA_FIELD_IDEAL_P(K) = NULL;
    CA_FIELD_IDEAL_LENGTH(K) = -1;
    CA_FIELD_IDEAL_ALLOC(K) = 0;
    CA_FIELD_HASH(K) = CA_EXT_HASH(ext);
}

void
ca_field_init_const(ca_field_t K, calcium_func_code func, ca_ctx_t ctx)
{
    ca_ext_ptr ext;
    ca_ext_t tmp;

    /* todo: avoid unnecessary work */
    ca_ext_init_const(tmp, func, ctx);
    ext = ca_ext_cache_insert(CA_CTX_EXT_CACHE(ctx), tmp, ctx);
    ca_ext_clear(tmp, ctx);

    CA_FIELD_LENGTH(K) = 1;
    CA_FIELD_EXT(K) = flint_malloc(sizeof(ca_ext_ptr));
    CA_FIELD_EXT_ELEM(K, 0) = ext;
    CA_FIELD_IDEAL_P(K) = NULL;
    CA_FIELD_IDEAL_LENGTH(K) = 0;
    CA_FIELD_IDEAL_ALLOC(K) = 0;
    CA_FIELD_HASH(K) = CA_EXT_HASH(ext);

    _ca_ctx_init_mctx(ctx, 1);
}

void ca_field_init_fx(ca_field_t K, calcium_func_code func, const ca_t x, ca_ctx_t ctx)
{
    ca_ext_ptr ext;
    ca_ext_t tmp;

    /* todo: avoid unnecessary work */
    ca_ext_init_fx(tmp, func, x, ctx);
    ext = ca_ext_cache_insert(CA_CTX_EXT_CACHE(ctx), tmp, ctx);
    ca_ext_clear(tmp, ctx);

    CA_FIELD_LENGTH(K) = 1;
    CA_FIELD_EXT(K) = flint_malloc(sizeof(ca_ext_ptr));
    CA_FIELD_EXT_ELEM(K, 0) = ext;
    CA_FIELD_IDEAL_P(K) = NULL;
    CA_FIELD_IDEAL_LENGTH(K) = 0;
    CA_FIELD_IDEAL_ALLOC(K) = 0;
    CA_FIELD_HASH(K) = CA_EXT_HASH(ext);

    _ca_ctx_init_mctx(ctx, 1);
}

void ca_field_init_fxy(ca_field_t K, calcium_func_code func, const ca_t x, const ca_t y, ca_ctx_t ctx)
{
    ca_ext_ptr ext;
    ca_ext_t tmp;

    /* todo: avoid unnecessary work */
    ca_ext_init_fxy(tmp, func, x, y, ctx);
    ext = ca_ext_cache_insert(CA_CTX_EXT_CACHE(ctx), tmp, ctx);
    ca_ext_clear(tmp, ctx);

    CA_FIELD_LENGTH(K) = 1;
    CA_FIELD_EXT(K) = flint_malloc(sizeof(ca_ext_ptr));
    CA_FIELD_EXT_ELEM(K, 0) = ext;
    CA_FIELD_IDEAL_P(K) = NULL;
    CA_FIELD_IDEAL_LENGTH(K) = 0;
    CA_FIELD_IDEAL_ALLOC(K) = 0;
    CA_FIELD_HASH(K) = CA_EXT_HASH(ext);

    _ca_ctx_init_mctx(ctx, 2);
}

void
ca_field_init_multi(ca_field_t K, slong len, ca_ctx_t ctx)
{
    CA_FIELD_LENGTH(K) = len;
    CA_FIELD_EXT(K) = flint_malloc(len * sizeof(ca_ext_ptr));
    CA_FIELD_IDEAL_P(K) = NULL;
    CA_FIELD_IDEAL_LENGTH(K) = 0;
    CA_FIELD_IDEAL_ALLOC(K) = 0;
    CA_FIELD_HASH(K) = 0;

    _ca_ctx_init_mctx(ctx, len);
}
