/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

void fmpq_mpoly_make_monic_inplace(fmpq_mpoly_t poly1, const fmpq_mpoly_ctx_t ctx)
{
    if (poly1->zpoly->length == 0)
        flint_throw(FLINT_ERROR, "zero polynomial in fmpq_mpoly_make_monic_inplace");

    fmpz_one(fmpq_numref(poly1->content));
    fmpz_set(fmpq_denref(poly1->content), poly1->zpoly->coeffs + 0);
}

