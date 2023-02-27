/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include "fmpz_mod_mpoly.h"

int fmpz_mod_mpoly_set_str_pretty(
    fmpz_mod_mpoly_t poly,
    const char * str,
    const char** x,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int ret;
    slong i;
    fmpz_mod_mpoly_t val;
    mpoly_parse_t E;
    char dummy[FLINT_BITS/2];

    mpoly_void_ring_init_fmpz_mod_mpoly_ctx(E->R, ctx);
    mpoly_parse_init(E);

    fmpz_mod_mpoly_init(val, ctx);
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        fmpz_mod_mpoly_gen(val, i, ctx);
        if (x == NULL)
        {
            flint_sprintf(dummy, "x%wd", i + 1);
            mpoly_parse_add_terminal(E, dummy, (const void *)val);
        }
        else
        {
            mpoly_parse_add_terminal(E, x[i], (const void *)val);
        }
    }
    fmpz_mod_mpoly_clear(val, ctx);

    ret = mpoly_parse_parse(E, poly, str, strlen(str));

    mpoly_parse_clear(E);

    return ret;
}
