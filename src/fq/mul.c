/*
    Copyright (C) 2011, 2012 Sebastian Pancratz 
    Copyright (C) 2012 Andres Goens

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq.h"

void
fq_mul(fq_t rop, const fq_t op1, const fq_t op2, const fq_ctx_t ctx)
{
    fmpz_poly_mul(rop, op1, op2);

    fq_reduce(rop, ctx);
}
