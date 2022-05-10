/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "fq_nmod.h"

int
fq_nmod_next(fq_nmod_t alpha, const fq_nmod_ctx_t fqctx)
{
    slong i;
    slong deg = nmod_poly_degree(fqctx->modulus);

    for (i = 0; i < deg; i++)
    {
        ulong c = nmod_poly_get_coeff_ui(alpha, i);
        c += UWORD(1);
        if (c < fqctx->mod.n)
        {
            nmod_poly_set_coeff_ui(alpha, i, c);
            return 1;
        }
        nmod_poly_set_coeff_ui(alpha, i, 0);
    }

    return 0;
}

void
fq_nmod_next_not_zero(fq_nmod_t alpha, const fq_nmod_ctx_t fqctx)
{
    slong i;
    slong deg = fqctx->modulus->length - 1;

    for (i = 0; i < deg; i++) {
        ulong c = nmod_poly_get_coeff_ui(alpha, i);
        c += UWORD(1);
        if (c < fqctx->mod.n) {
            nmod_poly_set_coeff_ui(alpha, i, c);
            return;
        }
        nmod_poly_set_coeff_ui(alpha, i, UWORD(0));
    }

    /* we hit zero, so skip to 1 */
    nmod_poly_set_coeff_ui(alpha, 0, UWORD(1));
}
