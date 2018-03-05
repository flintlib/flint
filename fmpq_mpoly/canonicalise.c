/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"


void fmpq_mpoly_canonicalise(fmpq_mpoly_t poly, const fmpq_mpoly_ctx_t ctx)
{
    slong i;
    fmpz_t c;

    if (poly->zpoly->length == 0) {
        fmpq_zero(poly->content);
        return;
    }

    fmpz_init(c);
    _fmpz_vec_content(c, poly->zpoly->coeffs, poly->zpoly->length);
    if (fmpz_sgn(poly->zpoly->coeffs + 0) < 0)
        fmpz_neg(c, c);

    fmpq_mul_fmpz(poly->content, poly->content, c);
    for (i = 0; i < poly->zpoly->length; i++)
        fmpz_divexact(poly->zpoly->coeffs + i, poly->zpoly->coeffs + i, c);

    fmpz_clear(c);
}
