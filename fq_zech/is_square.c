/*
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech.h"

int
fq_zech_is_square(const fq_zech_t op1, const fq_zech_ctx_t ctx)
{
    if (fq_zech_is_zero(op1, ctx) || fq_zech_is_one(op1, ctx) || ctx->p == 2)
        return 1;

    return (op1->value & 1) == 0;
}
