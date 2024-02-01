/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "fq.h"

void
fq_bit_pack(fmpz_t f, const fq_t op, flint_bitcnt_t bit_size,
            const fq_ctx_t ctx)
{
    fmpz_poly_bit_pack(f, op, bit_size);
}
