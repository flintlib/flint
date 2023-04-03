/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq.h"

void
fq_sub_one(fq_t rop, const fq_t op1, const fq_ctx_t ctx)
{
    fq_t one;

    fq_init(one, ctx);
    fq_one(one, ctx);
    fq_sub(rop, op1, one, ctx);
    fq_clear(one, ctx);
}
