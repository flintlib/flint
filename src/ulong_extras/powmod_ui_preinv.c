/*
    Copyright (C) 2013, 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

ulong
n_powmod_ui_preinv(ulong a, ulong exp, ulong n, ulong ninv, ulong norm)
{
    ulong x;

    /* check normalisation */
    FLINT_ASSERT((n & (UWORD(1) << (FLINT_BITS - 1))) != 0);
    FLINT_ASSERT(a < n);
    FLINT_ASSERT(n != 0);

    if (exp == 0)
    {
        x = UWORD(1) << norm;

        /* anything modulo 1 is 0 */
        return x == n ? 0 : x;
    }

    if (a == 0)
        return 0;

    /* find first 1 in binary repn of exp */
    while ((exp & 1) == 0)
    {
        a = n_mulmod_preinv(a, a, n, ninv, norm);
        exp >>= 1;
    }

    x = a;

    /* binary exponentiation */
    while (exp >>= 1)
    {
        a = n_mulmod_preinv(a, a, n, ninv, norm);

        if (exp & 1)
            x = n_mulmod_preinv(x, a, n, ninv, norm);
    }

    return x;
}
