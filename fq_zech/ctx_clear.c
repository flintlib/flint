/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <string.h>

#include "fq_zech.h"

void
fq_zech_ctx_clear(fq_zech_ctx_t ctx)
{
    flint_free(ctx->zech_log_table);
    flint_free(ctx->prime_field_table);
    flint_free(ctx->eval_table);

    if (ctx->owns_fq_nmod_ctx)
    {
        fq_nmod_ctx_clear(ctx->fq_nmod_ctx);
        flint_free(ctx->fq_nmod_ctx);
    }
}
