/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

void _nmod_mpoly_realloc(mp_limb_t ** coeff, ulong ** exp,
                                             slong * alloc, slong len, slong N)
{
    (* coeff) = (mp_limb_t *) flint_realloc(* coeff, len*sizeof(mp_limb_t));
    (* exp) = (ulong *) flint_realloc(* exp, len*N*sizeof(ulong));    
    (* alloc) = len;
}

void nmod_mpoly_realloc(nmod_mpoly_t poly,
                                       slong alloc, const nmod_mpoly_ctx_t ctx)
{
    slong N;

    if (alloc == 0)             /* Clear up, reinitialise */
    {
        nmod_mpoly_clear(poly, ctx);
        nmod_mpoly_init(poly, ctx);
        return;
    }

    N = mpoly_words_per_exp(poly->bits, ctx->minfo);

    if (poly->alloc != 0)            /* Realloc */
    {
        nmod_mpoly_truncate(poly, alloc, ctx);

        poly->coeffs = (mp_limb_t *) flint_realloc(poly->coeffs, alloc*sizeof(mp_limb_t));
        poly->exps   = (ulong *) flint_realloc(poly->exps, alloc*N*sizeof(ulong));

    }
    else                        /* Nothing allocated already so do it now */
    {
        poly->coeffs = (mp_limb_t *) flint_malloc(alloc*sizeof(mp_limb_t));
        poly->exps   = (ulong *) flint_malloc(alloc*N*sizeof(ulong));
    }

    poly->alloc = alloc;
}

