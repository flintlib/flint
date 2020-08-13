/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod.h"

void
fq_nmod_gcdinv(fq_nmod_t rop, fq_nmod_t inv, const fq_nmod_t op,
               const fq_nmod_ctx_t ctx)
{
    nmod_poly_gcdinv(rop, inv, op, ctx->modulus);
}
