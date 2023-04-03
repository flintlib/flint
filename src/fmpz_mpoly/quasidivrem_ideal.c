/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void fmpz_mpoly_quasidivrem_ideal(fmpz_t scale, fmpz_mpoly_struct ** Q,
     fmpz_mpoly_t R, const fmpz_mpoly_t A, fmpz_mpoly_struct * const * B,
                                        slong len, const fmpz_mpoly_ctx_t ctx)
{
    /* TODO !!! */
    fmpz_mpoly_quasidivrem_ideal_heap(scale, Q, R, A, B, len, ctx);
}
