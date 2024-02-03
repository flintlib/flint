/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "fmpz.h"
#include "fq_nmod.h"

void fq_nmod_mul_fmpz(fq_nmod_t rop, const fq_nmod_t op, const fmpz_t x, const fq_nmod_ctx_t ctx)
{
    ulong rx;

    rx = fmpz_fdiv_ui(x, fq_nmod_ctx_prime(ctx));
    nmod_poly_scalar_mul_nmod(rop, op, rx);
}
