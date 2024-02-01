/*
    Copyright (C) 2023 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"
#include "fq_nmod.h"
#include "fq_nmod_poly_factor.h"

int FQ_NMOD_POLY_ITERATED_FROBENIUS_CUTOFF(const fq_nmod_ctx_t ctx, slong length)
{
    int result;
    fmpz_t q;
    fmpz_init(q);
    fq_nmod_ctx_order(q, ctx);
    if ( 2 * fmpz_sizeinbase(q, 2) < 3 * (n_sqrt(length) + 1))
        result = 1;
    else
        result = 0;
    fmpz_clear(q);
    return result;
}
