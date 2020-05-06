/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_zech.h"
#include <string.h>

char *
fq_zech_get_str_pretty(const fq_zech_t op, const fq_zech_ctx_t ctx)
{
    slong value = (slong)op->value;
    size_t len;
    char *s;

    if (value == 0)
        len = 2U;
    else if (value < LONG_MAX)
        len = n_clog(op->value + 1, 10) + 1U;
    else if (-value == value)
        len = n_clog(op->value, 10) + 2U;
    else
        len = n_clog(-value + 1, 10) + 2U;
    s = flint_malloc((len + strlen(ctx->fq_nmod_ctx->var) + 1U) * sizeof(char));
    flint_sprintf(s, "%s^%wd", ctx->fq_nmod_ctx->var, op->value);
    return s;
}
