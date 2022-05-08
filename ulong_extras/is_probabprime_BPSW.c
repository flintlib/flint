/*
    Copyright (C) 2008 Peter Shrimpton
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"

int
n_is_probabprime_BPSW(ulong n)
{
    if (n <= UWORD(1))
        return 0;

    if ((n & UWORD(1)) == UWORD(0))
    {
        if (n == UWORD(2))
            return 1;
        return 0;
    }

    if (((n % 10) == 3) || ((n % 10) == 7))
    {
        if (n_is_probabprime_fermat(n, 2) == 0)
            return 0;

        return n_is_probabprime_fibonacci(n);
    }
    else
    {
        ulong d;

        d = n - UWORD(1);
        while ((d & UWORD(1)) == UWORD(0))
            d >>= 1;

        if (FLINT_BIT_COUNT(n) <= FLINT_D_BITS)
        {
            double npre = n_precompute_inverse(n);
            if (n_is_strong_probabprime_precomp(n, npre, WORD(2), d) == 0)
                return 0;
        }
        else
        {
            ulong ninv = n_preinvert_limb(n);
            if (n_is_strong_probabprime2_preinv(n, ninv, WORD(2), d) == 0)
                return 0;
        }

        return (n_is_probabprime_lucas(n) == 1);
    }
}
