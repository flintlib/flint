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
   Method of Niels Moller and Torbjorn Granlund, see "Improved Division by
   Invariant Integers: (Algorithm 4)
   https://gmplib.org/~tege/division-paper.pdf 
*/

ulong
n_ll_mod_preinv(ulong a_hi, ulong a_lo, ulong n, ulong ninv)
{
    ulong q0, q1, r, norm;

    FLINT_ASSERT(n != 0);

    count_leading_zeros(norm, n);

    /* reduce a_hi modulo n */
    if (a_hi >= n)
    {
        const ulong u1 = r_shift(a_hi, FLINT_BITS - norm);
        const ulong u0 = (a_hi << norm);

        n <<= norm;

        umul_ppmm(q1, q0, ninv, u1);
        add_ssaaaa(q1, q0, q1, q0, u1, u0);

        a_hi = (u0 - (q1 + 1) * n);

        if (a_hi > q0)
            a_hi += n;

        if (a_hi >= n)
            a_hi -= n;
    }
    else
    {
        n <<= norm;
        a_hi <<= norm;
    }

    /* now reduce the rest of the way */
    {
        const ulong u1 = a_hi + r_shift(a_lo, FLINT_BITS - norm);
        const ulong u0 = (a_lo << norm);

        umul_ppmm(q1, q0, ninv, u1);
        add_ssaaaa(q1, q0, q1, q0, u1, u0);

        r = (u0 - (q1 + 1) * n);

        if (r > q0)
            r += n;

        return (r < n) ? (r >> norm) : ((r - n) >> norm);
    }
}
