/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod.h"

int fq_nmod_get_fmpz(fmpz_t a, const fq_nmod_t b, const fq_nmod_ctx_t ctx)
{
    if (b->length > 1)
        return 0;

    if (b->length == 1)
        fmpz_set_ui(a, b->coeffs[0]);
    else
        fmpz_zero(a);

    return 1;
}
