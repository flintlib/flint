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
ca_ctx_clear(ca_ctx_t ctx)
{
    slong i;

    CA_INFO(ctx, ("%wd extension numbers cached at time of destruction\n", CA_CTX_EXT_CACHE(ctx)->length));
    CA_INFO(ctx, ("%wd fields cached at time of destruction\n", CA_CTX_FIELD_CACHE(ctx)->length));

    ca_ext_cache_clear(CA_CTX_EXT_CACHE(ctx), ctx);
    ca_field_cache_clear(CA_CTX_FIELD_CACHE(ctx), ctx);

    for (i = 0; i < ctx->mctx_len; i++)
        flint_free(ctx->mctx[i]);

    flint_free(ctx->mctx);
    flint_free(ctx->options);
}

