/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly_factor.h"

void fmpq_mpoly_factor_fit_length(fmpq_mpoly_factor_t f, slong len,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    if (len > f->alloc)
    {
        len = FLINT_MAX(len, f->alloc + f->alloc/2);
        fmpq_mpoly_factor_realloc(f, len, ctx);
    }
}
