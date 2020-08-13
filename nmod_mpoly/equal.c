/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

int _nmod_mpoly_equal(const mp_limb_t * coeff1, const ulong * exp1,
                      const mp_limb_t * coeff2, const ulong * exp2,
                                                            slong len, slong N)
{
    slong i;

    if (coeff1 != coeff2)
    {
        for (i = 0; i < len; i++)
            if (coeff1[i] != coeff2[i])
                return 0;
    }

    if (exp1 != exp2)
    {
        for (i = 0; i < len; i++)
            if (!mpoly_monomial_equal(exp1 + N*i, exp2 + N*i, N))
                return 0;
    }

    return 1;
}

int nmod_mpoly_equal(const nmod_mpoly_t poly1, const nmod_mpoly_t poly2,
                                                    const nmod_mpoly_ctx_t ctx)
{
    ulong * ptr1 = poly1->exps, * ptr2 = poly2->exps;
    slong max_bits, N;
    int r, free1 = 0, free2 = 0;

    if (poly1 == poly2)
        return 1;

    if (poly1->length != poly2->length)
        return 0;

    max_bits = FLINT_MAX(poly1->bits, poly2->bits);
    N = mpoly_words_per_exp(max_bits, ctx->minfo);

    if (max_bits > poly1->bits)
    {
        free1 = 1;
        ptr1 = (ulong *) flint_malloc(N*poly1->length*sizeof(ulong));
        mpoly_repack_monomials(ptr1, max_bits, poly1->exps, poly1->bits,
                                                    poly1->length, ctx->minfo);
    }

    if (max_bits > poly2->bits)
    {
        free2 = 1;
        ptr2 = (ulong *) flint_malloc(N*poly2->length*sizeof(ulong));
        mpoly_repack_monomials(ptr2, max_bits, poly2->exps, poly2->bits,
                                                    poly2->length, ctx->minfo);
    }

    r = _nmod_mpoly_equal(poly1->coeffs, ptr1,
                                        poly2->coeffs, ptr2, poly2->length, N);

    if (free1)
        flint_free(ptr1);

    if (free2)
        flint_free(ptr2);

    return r;
}
