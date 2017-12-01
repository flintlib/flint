/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

void nmod_mpoly_pow_rmul(nmod_mpoly_t poly1, const nmod_mpoly_t poly2,
                                         slong pow, const nmod_mpoly_ctx_t ctx)
{
    nmod_mpoly_t t;

    nmod_mpoly_init(t, ctx);

    if (poly1 == poly2)
    {
        nmod_mpoly_pow_rmul(t, poly2, pow, ctx);
        nmod_mpoly_swap(t, poly1, ctx);
        goto done;
    }

    nmod_mpoly_set_ui(poly1, 1, ctx);
    while (pow >= 1)
    { 
        nmod_mpoly_mul_johnson(t, poly1, poly2, ctx);
        if (pow == 1)
        {
            nmod_mpoly_swap(poly1, t, ctx);
            break;
        }
        nmod_mpoly_mul_johnson(poly1, t, poly2, ctx);
        pow -= 2;
    }

done:
    nmod_mpoly_clear(t, ctx);
}
