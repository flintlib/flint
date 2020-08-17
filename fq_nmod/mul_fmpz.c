/*
    Copyright (C) 2012 Sebastian Pancratz 
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod.h"

void fq_nmod_mul_fmpz(fq_nmod_t rop, const fq_nmod_t op, const fmpz_t x, const fq_nmod_ctx_t ctx)
{
    fmpz_t rx;
    fmpz_init(rx);

    fmpz_mod(rx, x, fq_nmod_ctx_prime(ctx));
    nmod_poly_scalar_mul_nmod(rop, op, fmpz_get_ui(rx));

    fmpz_clear(rx);
}
