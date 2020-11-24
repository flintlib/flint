/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

/* this file does not need to change with new orderings */

/* try to prove A is not a square */
int mpoly_is_proved_not_square(
    const ulong * Aexps,
    slong Alen,
    flint_bitcnt_t Abits,
    slong N,
    ulong * t)  /* temp of size N */
{
    slong i;

    if (Alen < 1)
        return 0;

    /* check for odd degrees & check total degree too in degree orderings */
    mpoly_monomial_set(t, Aexps + N*0, N);
    if (Abits <= FLINT_BITS)
    {
        ulong mask = mpoly_overflow_mask_sp(Abits);

        for (i = 1; i < Alen; i++)
            mpoly_monomial_max(t, t, Aexps + N*i, Abits, N, mask);

        return !mpoly_monomial_halves(t, t, N, mask);
    }
    else
    {
        for (i = 1; i < Alen; i++)
            mpoly_monomial_max_mp(t, t, Aexps + N*i, Abits, N);

        return !mpoly_monomial_halves_mp(t, t, N, Abits);
    }
}

