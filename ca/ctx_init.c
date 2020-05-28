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
ca_ctx_init(ca_ctx_t ctx)
{
    qqbar_t onei;

    /* Always create QQ, QQ(i) */

    ctx->fields = (ca_field_struct *) flint_malloc(2 * sizeof(ca_field_struct));
    ctx->fields_len = 2;
    ctx->fields_alloc = 2;

    ctx->extensions = (ca_extension_struct *) flint_malloc(1 * sizeof(ca_extension_struct));
    ctx->extensions_len = 1;
    ctx->extensions_alloc = 1;

    qqbar_init(onei);
    qqbar_i(onei);
    ca_extension_init_qqbar(ctx->extensions, onei);
    qqbar_clear(onei);

    ca_field_init_qq(ctx->fields);
    ca_field_init_nf(ctx->fields + 1, ctx->extensions);
}

