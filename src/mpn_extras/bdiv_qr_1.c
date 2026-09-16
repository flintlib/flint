/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

/*
    Single-limb Hensel division (port of radix_divmod_bn_1 with B = 2^64):
    develops n limbs of q with q * b == a (mod B^n), least significant limb
    first, for odd b. With binv = b^(-1) mod B, a single-limb carry c
    (initially 0) is maintained such that after limb i,

        q_partial * b == (a mod B^(i+1)) + c B^(i+1).

    At limb i, w = a_i - c, q_i = w binv, and the new carry is the high limb
    of c + q_i b (whose low limb is a_i).

    q may alias a. If r != NULL, r[0] receives (a - q b) / B^n mod B, so that
    a == q b + B^n r (mod B^(n+1)).
*/
void
flint_mpn_bdiv_qr_1(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t an,
    mp_limb_t b, mp_size_t n)
{
    mp_limb_t binv, cy, hi, lo, w, qi;
    mp_size_t i;

    FLINT_ASSERT(an >= 1);
    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(b & 1);

    binv = n_binvert(b);
    cy = 0;

    for (i = 0; i < n; i++)
    {
        w = ((i < an) ? a[i] : 0) - cy;
        qi = w * binv;
        q[i] = qi;
        umul_ppmm(hi, lo, qi, b);
        add_ssaaaa(hi, lo, hi, lo, 0, cy);
        FLINT_ASSERT(lo == ((i < an) ? a[i] : 0));
        cy = hi;
    }

    if (r != NULL)
        r[0] = ((n < an) ? a[n] : 0) - cy;
}
