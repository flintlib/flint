/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech.h"

void
fq_zech_randtest(fq_zech_t rop, flint_rand_t state, const fq_zech_ctx_t ctx)
{
    rop->value = n_randint(state, ctx->qm1 + 1);
}

void
fq_zech_randtest_not_zero(fq_zech_t rop, flint_rand_t state,
                          const fq_zech_ctx_t ctx)
{
    rop->value = n_randint(state, ctx->qm1);
}
