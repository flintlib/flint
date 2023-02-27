/*
    Copyright (C) 2012 Sebastian Pancratz 
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq.h"

void
fq_mul_ui(fq_t rop, const fq_t op, ulong x, const fq_ctx_t ctx)
{
    fmpz_poly_scalar_mul_ui(rop, op, x);

    fq_reduce(rop, ctx);
}
