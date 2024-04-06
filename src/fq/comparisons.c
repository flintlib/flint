/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "fq.h"

int fq_equal(const fq_t op1, const fq_t op2, const fq_ctx_t FLINT_UNUSED(ctx))
{
    return fmpz_poly_equal(op1, op2);
}

int fq_is_zero(const fq_t op, const fq_ctx_t FLINT_UNUSED(ctx))
{
    return fmpz_poly_is_zero(op);
}

int fq_is_one(const fq_t op, const fq_ctx_t FLINT_UNUSED(ctx))
{
    return fmpz_poly_is_one(op);
}
