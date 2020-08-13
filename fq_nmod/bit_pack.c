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
fq_nmod_bit_pack(fmpz_t f, const fq_nmod_t op, flint_bitcnt_t bit_size,
                 const fq_nmod_ctx_t ctx)
{
    nmod_poly_bit_pack(f, op, bit_size);
}

