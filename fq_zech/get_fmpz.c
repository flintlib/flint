/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech.h"

int fq_zech_get_fmpz(fmpz_t a, const fq_zech_t b, const fq_zech_ctx_t ctx)
{
    ulong q = ctx->eval_table[b->value];
    if (q >= ctx->p)
        return 0;

    fmpz_set_ui(a, q);
    return 1;
}

