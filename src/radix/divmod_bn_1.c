/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "radix.h"

/*
    Single-limb Hensel division: the bn == 1 case of radix_divmod_bn_classical,
    with b passed directly as one limb.

    Develops n limbs of the Hensel quotient q with q * b == a (mod B^n), least
    significant limb first. Returns 1 on success and 0 if b is not invertible
    modulo B (b not coprime to the radix), in which case nothing is written.

    This is the Hensel analogue of radix_divrem_1, and the dual of radix_mul_1:
    rather than propagating a carry while multiplying, we choose each quotient
    limb to cancel the current low limb, and propagate the resulting carry.

    With binv = b^(-1) mod B, maintain a single-limb carry c (initially 0) so
    that, after limb i, q_partial * b == (a mod B^{i+1}) + c * B^{i+1}. At limb i,

        w   = (a_i - c) mod B          (the limb of a still to be cancelled)
        q_i = w * binv mod B           (so that c + q_i*b == a_i (mod B))
        c   = (c + q_i*b) div B        (new carry; (c + q_i*b) mod B == a_i)

    The body is exactly radix_mul_1's "accumulate q_i*b into the carry and
    reduce" step, with two extra single-limb operations (the subtraction giving
    w and the multiplication giving q_i) per limb.

    q receives n limbs (room for n) and may alias a (a_i is read before q_i is
    written, and limbs of a at positions >= n are untouched).

    If rem != NULL it receives the resume remainder (a - q*b)/B^n mod B in
    rem[0]; after the loop c equals (q*b - (a mod B^n)) / B^n, so the remainder
    limb is (a_n - c) mod B -- the same "w" the next iteration would form. Then
    a == q*b + B^n*rem (mod B^{n+1}); exact division (with q of n limbs) holds iff
    rem[0] == 0 and a has no nonzero limbs at positions > n.
*/
int
radix_divmod_bn_1(nn_ptr q, nn_ptr rem, nn_srcptr a, slong an, ulong b,
    slong n, const radix_t radix)
{
    nmod_t mod = radix->B;
    const ulong B = mod.n;
    ulong binv, bb, k, m, r, hi, lo, cy = 0;
    slong i;

    FLINT_ASSERT(an >= 1);
    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(b < mod.n);

    /* binv = b^(-1) mod B (single limb); 0 iff b is a non-unit. */
    bb = b;
    if (!radix_invmod_bn(&binv, &bb, 1, 1, radix))
        return 0;

    /*
        The carry without a division by B: binv b = 1 + k B with k < B. For
        w < B, write w binv = m B + q_i with q_i = w binv mod B; then
        q_i b = w (1 + k B) - m b B = w + (w k - m b) B, so that
        floor(q_i b / B) = w k - m b, exactly, computed modulo 2^FLINT_BITS.
        With w = (a_i - c) mod B and the borrow t = [a_i < c], c + w = a_i + t B
        and c + q_i b = a_i + (t + w k - m b) B: the next carry is
        t + w k - m b <= B - 1. The division of w binv by B gives both q_i
        and m, and only that division remains on the chain.
    */
    umul_ppmm(hi, lo, binv, b);

    if (2 * NMOD_BITS(mod) <= FLINT_BITS)
    {
        /* w binv < B^2 fits a word */
        r = n_divrem_preinv_unnorm(&k, lo, B, mod.ninv, mod.norm);
        FLINT_ASSERT(hi == 0 && r == 1);

        for (i = 0; i < n; i++)
        {
            ulong ai = (i < an) ? a[i] : 0;
            ulong t = (ai < cy);
            ulong w = ai - cy + (B & -t);
            q[i] = n_divrem_preinv_unnorm(&m, w * binv, B, mod.ninv, mod.norm);
            cy = w * k - m * b + t;
        }
    }
    else if (mod.norm >= 2)
    {
        /* w binv < B^2 <= 2^(nu + FLINT_BITS - 2) */
        NMOD_DIVREM2_BOUNDED(k, r, hi, lo, mod);
        FLINT_ASSERT(r == 1);

        for (i = 0; i < n; i++)
        {
            ulong ai = (i < an) ? a[i] : 0;
            ulong t = (ai < cy);
            ulong w = ai - cy + (B & -t);
            umul_ppmm(hi, lo, w, binv);
            NMOD_DIVREM2_BOUNDED(m, q[i], hi, lo, mod);
            cy = w * k - m * b + t;
        }
    }
    else if (mod.norm == 1)
    {
        const ulong d = B << 1;

        udiv_qrnnd_preinv(k, r, (hi << 1) | (lo >> (FLINT_BITS - 1)), lo << 1, d, mod.ninv);
        FLINT_ASSERT((r >> 1) == 1);

        for (i = 0; i < n; i++)
        {
            ulong ai = (i < an) ? a[i] : 0;
            ulong t = (ai < cy);
            ulong w = ai - cy + (B & -t);
            umul_ppmm(hi, lo, w, binv);
            udiv_qrnnd_preinv(m, r, (hi << 1) | (lo >> (FLINT_BITS - 1)), lo << 1, d, mod.ninv);
            q[i] = r >> 1;
            cy = w * k - m * b + t;
        }
    }
    else
    {
        udiv_qrnnd_preinv(k, r, hi, lo, B, mod.ninv);
        FLINT_ASSERT(r == 1);

        for (i = 0; i < n; i++)
        {
            ulong ai = (i < an) ? a[i] : 0;
            ulong t = (ai < cy);
            ulong w = ai - cy + (B & -t);
            umul_ppmm(hi, lo, w, binv);
            udiv_qrnnd_preinv(m, q[i], hi, lo, B, mod.ninv);
            cy = w * k - m * b + t;
        }
    }

    (void) r;
    FLINT_ASSERT(cy < mod.n);

    if (rem != NULL)
        rem[0] = nmod_sub((n < an) ? a[n] : 0, cy, mod);

    return 1;
}
