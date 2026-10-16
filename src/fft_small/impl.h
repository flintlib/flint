/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FFT_SMALL_IMPL_H
#define FFT_SMALL_IMPL_H

#include "fft_small.h"

#ifdef __cplusplus
extern "C" {
#endif

/* The index i such that n = R->ffts[i].mod.n, one of the primes of the
   context R, or -1 if n is none of them. Modulo such a prime, fft_small
   multiplies with a single prime and no chinese remaindering. The primes
   of a context all have 50 bits. */
FLINT_FORCE_INLINE slong
_fft_small_mpn_ctx_prime_index(const mpn_ctx_struct * R, ulong n)
{
    slong i;

    if (FLINT_BIT_COUNT(n) != 50)
        return -1;

    for (i = 0; i < MPN_CTX_NCRTS; i++)
        if (R->ffts[i].mod.n == n)
            return i;

    return -1;
}

#ifdef __cplusplus
}
#endif

#endif
