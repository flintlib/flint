/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_poly.h"
#include "mpn_extras.h"
#include "nmod_vec.h"
#include "mpoly.h"

void n_polyu_clear(n_polyu_t A)
{
    if (A->alloc > 0)
    {
        flint_free(A->exps);
        flint_free(A->coeffs);
    }
    else
    {
        FLINT_ASSERT(A->exps == NULL);
        FLINT_ASSERT(A->coeffs == NULL);
    }
}

void n_polyu_realloc(n_polyu_t A, slong len)
{
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(len, old_alloc + 1 + old_alloc/2);

    FLINT_ASSERT(A->alloc >= 0);
    if (len <= A->alloc)
        return;

    if (old_alloc > 0)
    {
        A->exps = (ulong *) flint_realloc(A->exps, new_alloc*sizeof(ulong));
        A->coeffs = (mp_limb_t *) flint_realloc(A->coeffs, new_alloc*sizeof(mp_limb_t));
    }
    else
    {
        FLINT_ASSERT(A->exps == NULL);
        FLINT_ASSERT(A->coeffs == NULL);
        A->exps = (ulong *) flint_malloc(new_alloc*sizeof(ulong));
        A->coeffs = (mp_limb_t *) flint_malloc(new_alloc*sizeof(mp_limb_t));
    }

    A->alloc = new_alloc;
}

void n_polyu3_print_pretty(
    const n_polyu_t A,
    const char * var0,
    const char * var1,
    const char * var2)
{
    slong i;
    int first = 1;

    for (i = 0; i < A->length; i++)
    {
        if (!first)
            printf(" + ");
        first = 0;
        flint_printf("%wu*%s^%wu*%s^%wu*%s^%wu",
            A->coeffs[i],
            var0, extract_exp(A->exps[i], 2, 3),
            var1, extract_exp(A->exps[i], 1, 3),
            var2, extract_exp(A->exps[i], 0, 3));
    }

    if (first)
        flint_printf("0");
}

void n_polyu3_degrees(
    slong * deg0,
    slong * deg1,
    slong * deg2,
    const n_polyu_t A)
{
    slong i;
    ulong m;
    ulong mask = mpoly_overflow_mask_sp(FLINT_BITS/3);

    if (A->length <= 0)
    {
        *deg0 = *deg1 = *deg2 = -1;
        return;
    }

    m = A->exps[0];
    for (i = 1; i < A->length; i++)
        m = mpoly_monomial_max1(m, A->exps[i], FLINT_BITS/3, mask);

    *deg0 = extract_exp(m, 2, 3);
    *deg1 = extract_exp(m, 1, 3);
    *deg2 = extract_exp(m, 0, 3);
}

