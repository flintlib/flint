/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fq_zech.h"
#include "fq_zech_poly_factor.h"

int FQ_ZECH_POLY_ITERATED_FROBENIUS_CUTOFF(const fq_zech_ctx_t ctx, slong length)
{
    int result;
    ulong q;

    q = fq_zech_ctx_order_ui(ctx);

    if (2 * FLINT_BIT_COUNT(q) < 3 * (n_sqrt(length) + 1))
        result = 1;
    else
        result = 0;

    return result;
}
