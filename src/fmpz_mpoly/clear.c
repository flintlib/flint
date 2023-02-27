/*
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"

void fmpz_mpoly_clear(fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx)
{
    slong i;

    for (i = 0; i < poly->alloc; i++)
        _fmpz_demote(poly->coeffs + i);

    flint_free(poly->coeffs);
    flint_free(poly->exps);
}
