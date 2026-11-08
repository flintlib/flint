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

/* Some carry loops have a variant that was measured faster on arm64
   (Apple M4) and slower on x86-64 (Zen 4, Ice Lake): the one-word path of
   radix_mul_1 and the three-limb path of radix_mulmid_classical. */
#if defined(__aarch64__) || defined(_M_ARM64)
# define RADIX_ARM64_CARRIES 1
#else
# define RADIX_ARM64_CARRIES 0
#endif

/* Inline version of flint_mpn_divrem_2_1_preinv_unnorm for speed (used by
   radix_mulmid_classical): divides {up, 2} by d of
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
        /* reached by long columns in radix_mulmid_classical */
        d <<= norm;
        r = (u1 >> (FLINT_BITS - norm));
        udiv_qrnnd_preinv(qp[1], r, r, (u1 << norm) | (u0 >> (FLINT_BITS - norm)), d, dinv);
    }

    udiv_qrnnd_preinv(qp[0], r, r, u0 << norm, d, dinv);
    return r >> norm;
}

/* Divides {up, 3} by B = mod.n of nu = FLINT_BITS - mod.norm bits, for
   3 <= nu <= FLINT_BITS - 1 and {up, 3} < 2^(nu + 2 FLINT_BITS - 2): the
   two top limbs, below 2^(nu + FLINT_BITS - 2), by NMOD_DIVREM2_BOUNDED,
   then the low limb by one Moller-Granlund step. Writes the quotient to
   {qp, 3} (qp[2] = 0; qp may be up) and returns the remainder. */
FLINT_FORCE_INLINE mp_limb_t
_radix_divrem_3_1_bounded(mp_ptr qp, mp_srcptr up, nmod_t mod)
{
    mp_limb_t q1, r1, r, u0;
    const unsigned int norm = mod.norm;

    FLINT_ASSERT(norm >= 1 && norm <= FLINT_BITS - 3);

    u0 = up[0];
    NMOD_DIVREM2_BOUNDED(q1, r1, up[2], up[1], mod);
    qp[2] = 0;
    qp[1] = q1;
    udiv_qrnnd_preinv(qp[0], r, (r1 << norm) | (u0 >> (FLINT_BITS - norm)), u0 << norm,
                      mod.n << norm, mod.ninv);
    return r >> norm;
}

#endif
