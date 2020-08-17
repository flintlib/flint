/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod.h"

void fq_nmod_mul_si(fq_nmod_t rop, const fq_nmod_t op, slong x, const fq_nmod_ctx_t ctx)
{
    mp_limb_t rx = x < 0 ? -x : x;
    
    nmod_poly_scalar_mul_nmod(rop, op, rx);
    
    if (x < 0)
        nmod_poly_neg(rop, rop);
}
