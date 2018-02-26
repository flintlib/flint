/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

void
fmpq_mpoly_get_coeff_fmpq(fmpq_t c, const fmpq_mpoly_t poly,
                                           slong n, const fmpq_mpoly_ctx_t ctx)
{
    if (n < poly->zpoly->length)
    {
        fmpq_mul_fmpz(c, poly->content, poly->zpoly->coeffs + n);
    } else
    {
        fmpq_zero(c);
    }
}
