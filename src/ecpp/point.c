/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5.1

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mod.h"
#include "ecpp.h"

/*
    Points in Jacobian coordinates (X : Y : Z), Z = 0 for the point at
    infinity. The arithmetic is in point_gr.c, over a generic ring.
*/

void
ecpp_point_init(ecpp_point_t P)
{
    fmpz_init(P->X);
    fmpz_init(P->Y);
    fmpz_init(P->Z);
}

void
ecpp_point_clear(ecpp_point_t P)
{
    fmpz_clear(P->X);
    fmpz_clear(P->Y);
    fmpz_clear(P->Z);
}

void
ecpp_point_set_affine(ecpp_point_t P, const fmpz_t x, const fmpz_t y)
{
    fmpz_set(P->X, x);
    fmpz_set(P->Y, y);
    fmpz_one(P->Z);
}

int
ecpp_point_is_zero(const ecpp_point_t P)
{
    return fmpz_is_zero(P->Z);
}
