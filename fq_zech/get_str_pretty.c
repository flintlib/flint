/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef strlen
# ifdef __GNUC__
#  define strlen __builtin_strlen
# else
#  include <string.h>
# endif
#endif
#include <stdio.h>
#include "flint-impl.h"
#include "ulong_extras.h"
#include "fq_zech.h"

char *
fq_zech_get_str_pretty(const fq_zech_t op, const fq_zech_ctx_t ctx)
{
    char * s;
    
    if (fq_zech_is_zero(op, ctx))
    {
        s = flint_malloc(2*sizeof(char));
        sprintf(s, "0");
    }
    else
    {
       slong num_chars = op->value == 0 ? 1 : n_clog(op->value + 1, 10);
       s = flint_malloc((num_chars + strlen(ctx->fq_nmod_ctx->var) + 2) * sizeof(char));
       sprintf(s, "%s^" WORD_FMT "d", ctx->fq_nmod_ctx->var, op->value);
    }
    
    return s;
}
