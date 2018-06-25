/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

void nmod_mpoly_clear(nmod_mpoly_t poly, const nmod_mpoly_ctx_t ctx)
{
    if (poly->alloc != 0)
    {
        flint_free(poly->coeffs);
        flint_free(poly->exps);
    }
}
