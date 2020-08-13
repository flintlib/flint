/*
    Copyright (C) 2009, 2010, 2015 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

/* 
   Method of Niels Moller and Torbjorn Granlund see paper:
   Improved Division by Invariant Integers: (algorithm 4)
   https://gmplib.org/~tege/division-paper.pdf 
*/

ulong
n_divrem2_preinv(ulong * q, ulong a, ulong n, ulong ninv)
{
    ulong norm, q1, q0, r;

    FLINT_ASSERT(n != 0);

    count_leading_zeros(norm, n);
    n <<= norm;

    {
        const ulong u1 = r_shift(a, FLINT_BITS - norm);
        const ulong u0 = (a << norm);

        umul_ppmm(q1, q0, ninv, u1);
        add_ssaaaa(q1, q0, q1, q0, u1, u0);

        (*q) = q1 + 1;
        r = u0 - (*q) * n;

        if (r > q0)
        {
            r += n;
            (*q)--;
        }

        if (r >= n)
        {
            (*q)++;
            r -= n;
        }

        return r >> norm;
    }
}
