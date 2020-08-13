/*
    Copyright (C) 2018, Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

ulong
n_powmod2_fmpz_preinv(ulong a, const fmpz_t exp, ulong n, ulong ninv)
{
    flint_bitcnt_t i, bits;
    ulong x, norm;

    FLINT_ASSERT(n != 0);
    FLINT_ASSERT(a < n);
    FLINT_ASSERT(fmpz_sgn(exp) >= 0);

    if (fmpz_is_zero(exp))
    {
        /* anything modulo 1 is 0 */
        return n == 1 ? 0 : 1;
    }

    if (a == 0)
        return 0;

    count_leading_zeros(norm, n);
    a <<= norm;
    n <<= norm;

    bits = fmpz_bits(exp);
    i = 0;

    while (i < bits && fmpz_tstbit(exp, i) == 0)
    {
        a = n_mulmod_preinv(a, a, n, ninv, norm);
        i++;
    }

    x = a;

    i++;
    while (i < bits)
    {
        a = n_mulmod_preinv(a, a, n, ninv, norm);
        if (fmpz_tstbit(exp, i) != 0) {
            x = n_mulmod_preinv(x, a, n, ninv, norm);
        }
        i++;
    }

    return x >> norm;
}
