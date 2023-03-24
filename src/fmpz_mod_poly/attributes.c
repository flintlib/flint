/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"

slong fmpz_mod_poly_degree(const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)
{
    return poly->length - 1;
}

slong fmpz_mod_poly_length(const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)
{
    return poly->length;
}

fmpz * fmpz_mod_poly_lead(const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)
{
    if (poly->length)
        return poly->coeffs + (poly->length - 1);
    else
        return NULL;
}

int fmpz_mod_poly_is_monic(const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)
{
    return f->length > 0 && fmpz_is_one(f->coeffs + f->length - 1);
}

int fmpz_mod_poly_is_one(const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)
{
   return poly->length == 1 && fmpz_is_one(poly->coeffs + 0);
}

int fmpz_mod_poly_is_gen(const fmpz_mod_poly_t op, const fmpz_mod_ctx_t ctx)
{
    return (op->length) == 2 && (*(op->coeffs + 1) == WORD(1)) && (*(op->coeffs + 0) == WORD(0));
}

/* bogus for non-prime modulus */
int fmpz_mod_poly_is_unit(const fmpz_mod_poly_t op, const fmpz_mod_ctx_t ctx)
{
    return (op->length == 1) && fmpz_mod_is_invertible(op->coeffs + 0, ctx);
}
