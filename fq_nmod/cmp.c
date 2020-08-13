/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod.h"

int fq_nmod_cmp(const fq_nmod_t a, const fq_nmod_t b, const fq_nmod_ctx_t ctx)
{
    slong i;

    if (a->length != b->length)
        return a->length < b->length ? -1 : 1;

    for (i = 0; i < a->length; i++)
    {
        if (a->coeffs[i] != b->coeffs[i])
            return a->coeffs[i] < b->coeffs[i] ? -1 : 1;
    }

    return 0;
}
