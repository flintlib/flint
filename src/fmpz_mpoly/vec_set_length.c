/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void
fmpz_mpoly_vec_set_length(fmpz_mpoly_vec_t vec, slong len, const fmpz_mpoly_ctx_t ctx)
{
    slong i;

    if (len > vec->length)
    {
        fmpz_mpoly_vec_fit_length(vec, len, ctx);
        for (i = vec->length; i < len; i++)
            fmpz_mpoly_zero(vec->p + i, ctx);
    }
    else if (len < vec->length)
    {
        for (i = len; i < vec->length; i++)
           fmpz_mpoly_zero(vec->p + i, ctx);
    }

    vec->length = len;
}
