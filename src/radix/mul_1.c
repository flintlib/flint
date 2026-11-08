/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "radix.h"
#include "impl.h"

ulong
radix_mul_1(nn_ptr z, nn_srcptr a, slong n, ulong c, const radix_t radix)
{
    slong i;
    ulong hi, lo;
    slong sbits;
    nmod_t mod = radix->B;

    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(c < mod.n);

    sbits = 2 * NMOD_BITS(mod);

    if (sbits <= FLINT_BITS)
    {
#if RADIX_ARM64_CARRIES
        /*
            Decoupled carries as below, with one-word products: a[i] c is
            below B^2 <= 2^FLINT_BITS and is split on its own as
            q_i B + r_i (q_i <= B - 2), and z[i] = (r_i + q_(i-1) + cb) mod B
            with a carry cb of 0 or 1; r_i + q_(i-1) + cb <= 2B - 2.
        */
        const ulong B = mod.n;
        ulong q, r, qp = 0, cb = 0;
        FLINT_ASSERT(mod.norm != 0);

        for (i = 0; i < n; i++)
        {
            r = n_divrem_preinv_unnorm(&q, a[i] * c, B, mod.ninv, mod.norm);
            r += qp + cb;
            cb = (r >= B);
            z[i] = r - (B & -cb);
            qp = q;
        }

        FLINT_ASSERT(qp + cb < B);
        return qp + cb;
#else
        ulong cy = 0;
        FLINT_ASSERT(mod.norm != 0);

        for (i = 0; i < n; i++)
        {
            cy += a[i] * c;
            z[i] = n_divrem_preinv_unnorm(&cy, cy, mod.n, mod.ninv, mod.norm);
        }

        FLINT_ASSERT(cy < mod.n);
        return cy;
#endif
    }
    else
    {
        /*
            Decoupled carries: each product is split on its own as
            a[i] c = q_i B + r_i, with q_i <= B - 2 since a[i] c <= (B - 1)^2,
            and z[i] = (r_i + q_(i-1) + cb) mod B, where the carry cb is 0
            or 1. The divisions are then independent, and the chain from
            one limb to the next is a few additions and comparisons.
        */
        const ulong B = mod.n;
        ulong q, r, qp = 0, cb = 0;

        if (mod.norm >= 2)
        {
            /* a[i] c < 2^(2 nu) <= 2^(nu + FLINT_BITS - 2), and
               r + qp + cb <= 2B - 2 fits a word */
            for (i = 0; i < n; i++)
            {
                umul_ppmm(hi, lo, a[i], c);
                NMOD_DIVREM2_BOUNDED(q, r, hi, lo, mod);
                r += qp + cb;
                cb = (r >= B);
                z[i] = r - (B & -cb);
                qp = q;
            }
        }
        else if (mod.norm == 1)
        {
            /* division of 2 a[i] c by the normalized 2B */
            const ulong d = B << 1;

            for (i = 0; i < n; i++)
            {
                umul_ppmm(hi, lo, a[i], c);
                udiv_qrnnd_preinv(q, r, (hi << 1) | (lo >> (FLINT_BITS - 1)), lo << 1, d, mod.ninv);
                r = (r >> 1) + qp + cb;
                cb = (r >= B);
                z[i] = r - (B & -cb);
                qp = q;
            }
        }
        else
        {
            /* r + qp + cb may not fit a word: compare r with B - (qp + cb),
               which is at least 1 */
            for (i = 0; i < n; i++)
            {
                umul_ppmm(hi, lo, a[i], c);
                udiv_qrnnd_preinv(q, r, hi, lo, B, mod.ninv);
                qp += cb;
                cb = (r >= B - qp);
                z[i] = r + qp - (B & -cb);
                qp = q;
            }
        }

        FLINT_ASSERT(qp + cb < B);
        return qp + cb;
    }
}
