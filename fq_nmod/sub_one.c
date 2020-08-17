/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod.h"

void fq_nmod_sub_one(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_ctx_t ctx)
{
    fq_nmod_t one;

    fq_nmod_init(one, ctx);
    fq_nmod_one(one, ctx);
    fq_nmod_sub(rop, op1, one, ctx);
    fq_nmod_clear(one, ctx);
}
