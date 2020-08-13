/*
    Copyright (C) 2011, 2012 Sebastian Pancratz 
    Copyright (C) 2012 Andres Goens
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod.h"

void fq_nmod_mul(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_t op2, const fq_nmod_ctx_t ctx)
{
    nmod_poly_mul(rop, op1, op2);

    fq_nmod_reduce(rop, ctx);
}
