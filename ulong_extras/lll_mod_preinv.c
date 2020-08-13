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
n_lll_mod_preinv(ulong a_hi, ulong a_mi, ulong a_lo, ulong n, ulong ninv)
{
    ulong q0, q1, r, norm;

    count_leading_zeros(norm, n);
    n <<= norm;

    /* 
       a_hi is already reduced, so first reduce a_hi, a_mi mod n 
     */
    {
        const ulong u1 = (a_hi << norm) + r_shift(a_mi, FLINT_BITS - norm);
        const ulong u0 = (a_mi << norm);

        umul_ppmm(q1, q0, ninv, u1);
        add_ssaaaa(q1, q0, q1, q0, u1, u0);

        a_mi = (u0 - (q1 + 1) * n);

        if (a_mi > q0)
            a_mi += n;

        if (a_mi >= n)
            a_mi -= n;
    }

    /* 
       a_mi is now reduced mod n, so reduce a_mi, a_lo mod n 
     */
    {
        const ulong u1 = a_mi + r_shift(a_lo, FLINT_BITS - norm);
        const ulong u0 = (a_lo << norm);

        umul_ppmm(q1, q0, ninv, u1);
        add_ssaaaa(q1, q0, q1, q0, u1, u0);

        r = (u0 - (q1 + 1) * n);

        if (r > q0)
            r += n;

        return (r < n) ? (r >> norm) : ((r - n) >> norm);
    }
}
