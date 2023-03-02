/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

void fmpz_mod_mpoly_clear(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;

    for (i = 0; i < A->coeffs_alloc; i++)
        fmpz_clear(A->coeffs + i);

    if (A->coeffs_alloc > 0)
        flint_free(A->coeffs);

    if (A->exps_alloc > 0)
        flint_free(A->exps);
}
