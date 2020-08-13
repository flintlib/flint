/*
    Copyright (C) 2009, 2010, 2012, 2016 William Hart

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
n_mulmod_preinv(ulong a, ulong b, ulong n, ulong ninv, ulong norm)
{
    ulong q0, q1, r, p_hi, p_lo;

    /* check normalisation */
    FLINT_ASSERT((n & (UWORD(1) << (FLINT_BITS - 1))) != 0);
    FLINT_ASSERT(a < n);
    FLINT_ASSERT(b < n);
    FLINT_ASSERT(n != 0);

    /* renormalise product */
    a >>= norm;

    /* multiply */
    umul_ppmm(p_hi, p_lo, a, b);

    /* reduce mod n */
    {
        umul_ppmm(q1, q0, ninv, p_hi);
        add_ssaaaa(q1, q0, q1, q0, p_hi, p_lo);

        r = (p_lo - (q1 + 1) * n);

        if (r > q0)
            r += n;

        return (r < n ? r : r - n);
    }
}
