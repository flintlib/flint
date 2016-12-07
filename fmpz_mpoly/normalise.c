/*
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"

void
_fmpz_mpoly_normalise(fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx)
{
    slong i;

    for (i = poly->length - 1; i >= 0 && poly->coeffs[i] == 0; i--) ;

    poly->length = i + 1;
}
