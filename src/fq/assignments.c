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

void fq_set(fq_t rop, const fq_t op, const fq_ctx_t ctx)
{
    fmpz_poly_set(rop, op);
}

void fq_set_fmpz(fq_t rop, const fmpz_t x, const fq_ctx_t ctx)
{
    fmpz_poly_set_fmpz(rop, x);
    fq_reduce(rop, ctx);
}

void fq_set_ui(fq_t rop, const ulong x, const fq_ctx_t ctx)
{
    fmpz_poly_set_ui(rop, x);
    fq_reduce(rop, ctx);
}

void fq_set_si(fq_t rop, const slong x, const fq_ctx_t ctx)
{
    fmpz_poly_set_si(rop, x);
    fq_reduce(rop, ctx);
}

void fq_zero(fq_t rop,  const fq_ctx_t ctx)
{
    fmpz_poly_zero(rop);
}

void fq_one(fq_t rop,  const fq_ctx_t ctx)
{
    fmpz_poly_one(rop);
}
