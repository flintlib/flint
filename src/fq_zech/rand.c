/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fq_zech.h"

void
fq_zech_rand(fq_zech_t rop, flint_rand_t state, const fq_zech_ctx_t ctx)
{
    rop->value = n_urandint(state, ctx->qm1 + 1);
}

void
fq_zech_rand_not_zero(fq_zech_t rop, flint_rand_t state, const fq_zech_ctx_t ctx)
{
    rop->value = n_urandint(state, ctx->qm1);
}
