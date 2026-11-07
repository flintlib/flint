/*
    Copyright (C) 2026 Fredrik Johansson
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Internal helpers of the radix module, included as "impl.h". */

#ifndef RADIX_IMPL_H
#define RADIX_IMPL_H

#include "radix.h"

/* Inline version of flint_mpn_divrem_2_1_preinv_unnorm for speed (used by
   radix_divmod_bn_1 and radix_mulmid_classical): divides {up, 2} by d of
   FLINT_BITS - norm bits, norm >= 1, with dinv = n_preinvert_limb(d << norm);
   writes the quotient to {qp, 2} (qp may be up) and returns the remainder. */
FLINT_FORCE_INLINE mp_limb_t
_radix_divrem_2_1_unnorm1(mp_ptr qp, mp_srcptr up, mp_limb_t d, mp_limb_t dinv, unsigned int norm)
{
    mp_limb_t u0, u1, r;

    FLINT_ASSERT(norm >= 1);

    u1 = up[1];
    u0 = up[0];
    if (u1 < d)
    {
        d <<= norm;
        qp[1] = 0;
        r = (u1 << norm) | (u0 >> (FLINT_BITS - norm));
    }
    else
    {
        /* Reached by long columns in radix_mulmid_classical. It is
           unreachable in radix_divmod_bn_1, but leaving it in gave faster
           code with GCC on Zen 3 (in the former radix_mul_1 loop). */
        d <<= norm;
        r = (u1 >> (FLINT_BITS - norm));
        udiv_qrnnd_preinv(qp[1], r, r, (u1 << norm) | (u0 >> (FLINT_BITS - norm)), d, dinv);
    }

    udiv_qrnnd_preinv(qp[0], r, r, u0 << norm, d, dinv);
    return r >> norm;
}

#endif
