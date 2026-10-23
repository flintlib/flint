/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

#if FLINT_HAVE_NATIVE_mpn_mulmid_n

void
flint_mpn_mulmid_via_n_padded(mp_ptr z, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn,
                 mp_size_t zlo, mp_size_t zhi)
{
    mp_size_t zn = zhi - zlo;
    mp_size_t pa, qb, La, Lb, zlo2, n, half, alpha, beta;
    mp_srcptr ap, bp;
    mp_ptr A, B, rp;
    TMP_INIT;

    FLINT_ASSERT(an >= 1);
    FLINT_ASSERT(bn >= 1);
    FLINT_ASSERT(zlo >= 0);
    FLINT_ASSERT(zhi > zlo);
    FLINT_ASSERT(zhi <= an + bn);

    /* (1) limbs of either input at or above zhi cannot reach the window */
    an = FLINT_MIN(an, zhi);
    bn = FLINT_MIN(bn, zhi);

    /* (2) drop low limbs that cannot contribute to any in-window diagonal */
    pa = (zlo > bn - 1) ? zlo - (bn - 1) : 0;
    qb = (zlo > an - 1) ? zlo - (an - 1) : 0;
    La = an - pa;
    Lb = bn - qb;

    /* nothing survives: every contributing product was dropped or out of range */
    if (La <= 0 || Lb <= 0)
    {
        flint_mpn_zero(z, zn);
        return;
    }

    ap = a + pa;
    bp = b + qb;
    zlo2 = zlo - pa - qb;       /* window floor in sliced coordinates, >= 0 */

    /* longer sliced operand goes in the 2n-1 slot */
    if (La < Lb)
    {
        FLINT_SWAP(mp_srcptr, ap, bp);
        FLINT_SWAP(mp_size_t, La, Lb);
    }

    /* (3) smallest balanced size covering the window and holding both operands */
    n = zn;
    if (Lb > n)
        n = Lb;
    half = (La + Lb - zlo2 + 1) / 2;        /* ceil((La + Lb - zlo2) / 2) */
    if (half > n)
        n = half;

    alpha = (Lb - 1) - zlo2;                /* offset of the long operand in A */
    beta = n - Lb;                          /* short operand sits atop B       */

    FLINT_ASSERT(alpha >= 0 && beta >= 0);
    FLINT_ASSERT(zn <= n);
    FLINT_ASSERT(alpha + La <= 2 * n - 1);
    FLINT_ASSERT(beta + Lb <= n);

    TMP_START;
    A = TMP_ARRAY_ALLOC(2 * n - 1, mp_limb_t);
    B = TMP_ARRAY_ALLOC(n, mp_limb_t);
    rp = TMP_ARRAY_ALLOC(n + 2, mp_limb_t);

    flint_mpn_zero(A, 2 * n - 1);
    flint_mpn_copyi(A + alpha, ap, La);
    flint_mpn_zero(B, n);
    flint_mpn_copyi(B + beta, bp, Lb);

    flint_mpn_mulmid_n(rp, A, B, n);

    /* output limb k is product position zlo + k; copy the window verbatim */
    flint_mpn_copyi(z, rp, zn);

    TMP_END;
}

#endif
