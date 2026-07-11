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

/*
    Set (z, zhi - zlo) to the limb window [zlo, zhi) of the integer product
    a * b, as a lower approximation: partial products a[p] * b[q] whose low
    limb lies below limb zlo (i.e. p + q < zlo) are dropped, so the carry they
    would propagate into limb zlo is not recovered.  This matches the contract
    of radix_mulmid_classical with radix B = 2^64 (and the schoolbook
    flint_mpn_mulmid_classical).  Equivalently

        z = ( sum_{p+q >= zlo} a[p]*b[q] * B^(p+q-zlo) ) mod B^(zhi-zlo).

    Requires an >= 1, bn >= 1 and 0 <= zlo < zhi <= an + bn.

    Implementation: reduce the window to a single *balanced* middle product and
    hand it to flint_mpn_mulmid_n (Karatsuba/Toom-42).  flint_mpn_mulmid_n
    computes, for inputs {A, 2n-1} and {B, n}, the band n-1 <= p'+q' < 2n-1 of
    A*B shifted down by n-1, dropping all products with p'+q' < n-1 -- which is
    exactly the same "drop everything below the floor" rule as our window's
    lower-approximation, so the deficits coincide and the result is the *same*
    value the schoolbook would produce (it is not merely another approximation).

    To map the window onto that shape we:

      (1) clamp an, bn to zhi (limbs at/above zhi cannot reach the window);

      (2) slice off the low limbs that can never contribute: a[p] with
          p + (bn-1) < zlo, and b[q] with q + (an-1) < zlo.  Every product
          removed here has p+q < zlo (already dropped) or low limb >= zhi
          (outside the window), so slicing changes nothing.  This also keeps
          the balanced size from being tied to the full an;

      (3) place the (longer) sliced operand of length La in the 2n-1 slot at
          offset alpha and the shorter one of length Lb in the n slot at the
          top (offset beta = n - Lb), with alpha + beta = n - 1 - zlo2 so that
          output limb k lands exactly at product position zlo + k.  Here zlo2
          is the window floor in sliced coordinates; one shows zlo2 <= Lb - 1,
          hence alpha = (Lb-1) - zlo2 >= 0.

    We choose the smallest n that fits both operands and covers the window:

        n = max( zn, Lb, ceil((La + Lb - zlo2) / 2) ).

    NOTE ON THE ALTERNATIVE.  This is the "generous n + zero-pad" route: one
    balanced call, no fix-up.  It is a good fit when the window is wide (zn
    comparable to Lb), i.e. exactly the regime where Karatsuba/Toom pays off.
    It is *not* a good fit for a narrow window sitting in a tall, skewed region:
    there n is forced up to ~Lb (B must hold all of b), so the balanced call
    costs ~O(Lb^2) while the window only depends on ~zn diagonals and plain
    schoolbook (flint_mpn_mulmid_classical) costs ~O(zn * Lb).  A
    production dispatcher should therefore fall back to schoolbook for narrow
    windows, or chunk the region the way the GMP-derived general middle product
    (flint_mpn_mulmid_unbalanced, mulmid_gmp.c) does with MULMID_CHUNK, and only enter this
    balanced reduction once a chunk is square enough to benefit.
*/
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
