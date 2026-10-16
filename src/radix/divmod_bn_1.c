/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "radix.h"

/* Inline version of flint_mpn_divrem_2_1_preinv_unnorm for speed
   (mirrors radix_mul_1). */
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
        d <<= norm;
        r = (u1 >> (FLINT_BITS - norm));
        udiv_qrnnd_preinv(qp[1], r, r, (u1 << norm) | (u0 >> (FLINT_BITS - norm)), d, dinv);
    }

    udiv_qrnnd_preinv(qp[0], r, r, u0 << norm, d, dinv);
    return r >> norm;
}

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
    ulong binv, bb;
    slong i;
    slong sbits;

    FLINT_ASSERT(an >= 1);
    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(b < mod.n);

    /* binv = b^(-1) mod B (single limb); 0 iff b is a non-unit. */
    bb = b;
    if (!radix_invmod_bn(&binv, &bb, 1, 1, radix))
        return 0;

    sbits = 2 * NMOD_BITS(mod);

    if (sbits <= FLINT_BITS)
    {
        ulong cy = 0;

        FLINT_ASSERT(mod.norm != 0);

        for (i = 0; i < n; i++)
        {
            ulong ai = (i < an) ? a[i] : 0;
            ulong w = nmod_sub(ai, cy, mod);
            ulong qi = nmod_mul(w, binv, mod);
            q[i] = qi;
            cy += qi * b;
            (void) n_divrem_preinv_unnorm(&cy, cy, mod.n, mod.ninv, mod.norm);
        }

        FLINT_ASSERT(cy < mod.n);

        if (rem != NULL)
            rem[0] = nmod_sub((n < an) ? a[n] : 0, cy, mod);
    }
    else if (mod.norm == 0)
    {
        ulong hi, lo, r;
        ulong cy0 = 0, cy1 = 0;

        for (i = 0; i < n; i++)
        {
            ulong ai = (i < an) ? a[i] : 0;
            ulong w = nmod_sub(ai, cy0, mod);
            ulong qi = nmod_mul(w, binv, mod);
            q[i] = qi;
            umul_ppmm(hi, lo, qi, b);
            add_ssaaaa(cy1, cy0, cy1, cy0, hi, lo);
            r = n_divrem_norm(&cy1, cy1, mod.n);
            udiv_qrnnd_preinv(cy0, r, r, cy0, mod.n, mod.ninv);
            (void) r;
        }

        FLINT_ASSERT(cy0 < mod.n);
        FLINT_ASSERT(cy1 == 0);

        if (rem != NULL)
            rem[0] = nmod_sub((n < an) ? a[n] : 0, cy0, mod);
    }
    else
    {
        ulong hi, lo;
        ulong cy[2] = { 0, 0 };

        for (i = 0; i < n; i++)
        {
            ulong ai = (i < an) ? a[i] : 0;
            ulong w = nmod_sub(ai, cy[0], mod);
            ulong qi = nmod_mul(w, binv, mod);
            q[i] = qi;
            umul_ppmm(hi, lo, qi, b);
            add_ssaaaa(cy[1], cy[0], cy[1], cy[0], hi, lo);
            (void) _radix_divrem_2_1_unnorm1(cy, cy, mod.n, mod.ninv, mod.norm);
        }

        FLINT_ASSERT(cy[0] < mod.n);
        FLINT_ASSERT(cy[1] == 0);

        if (rem != NULL)
            rem[0] = nmod_sub((n < an) ? a[n] : 0, cy[0], mod);
    }

    return 1;
}
