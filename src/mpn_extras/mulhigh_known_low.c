/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

/* Below this product length the low product is formed with the
   assembly-optimised flint_mpn_mullow_n on zero-padded operands rather
   than with a windowed middle product. */
#ifndef FLINT_MPN_MULHIGH_KNOWN_LOW_MULLOW_CUTOFF
#define FLINT_MPN_MULHIGH_KNOWN_LOW_MULLOW_CUTOFF 12
#endif

/*
    out[0 .. khi-klo) = limbs [klo, khi) of the product x*y, given that the
    low klo limbs of x*y are already known and equal (kl, kl_len) (limbs at
    positions >= kl_len being zero).

    Port of _radix_mulhigh_known_low with B = 2^FLINT_BITS. When klo > 3 the
    high limbs are obtained from the windowed middle product over
    [klo-3, khi) with 3 guard limbs, corrected using the known low limbs:
    the middle product is a lower approximation whose deficit (a single
    carry from below the window) is bounded by min(xn, yn) B < B^2, so after
    subtracting kl[klo-3 .. klo) the guard limbs have true value zero, and the
    computed high part is either exact (guard limb 2 equal to zero) or one
    unit too small (guard limb 2 equal to B-1). Since the limb radix is
    always huge relative to the operand lengths, the radix module's
    LIMB_RADIX >= min(xn, yn) condition is unconditionally satisfied.

    If the product is shorter than khi limbs (xn + yn < khi), the window is
    computed up to xn + yn and zero-padded; if it does not reach klo at all,
    the output is zero.

    Otherwise (klo <= 3) a full low product to khi limbs is taken; for short
    products this uses flint_mpn_mullow_n.

    scratch must have room for khi limbs.
*/
void
_flint_mpn_mulhigh_known_low(mp_ptr out, mp_srcptr x, mp_size_t xn,
    mp_srcptr y, mp_size_t yn, mp_srcptr kl, mp_size_t kl_len,
    mp_size_t klo, mp_size_t khi, mp_ptr scratch)
{
    mp_size_t outn = khi - klo;
    mp_size_t pn = xn + yn;

    FLINT_ASSERT(xn >= 1);
    FLINT_ASSERT(yn >= 1);
    FLINT_ASSERT(0 <= klo);
    FLINT_ASSERT(klo < khi);

    if (klo > 3 && khi > FLINT_MPN_MULHIGH_KNOWN_LOW_MULLOW_CUTOFF)
    {
        mp_size_t hi = FLINT_MIN(khi, pn);
        mp_size_t wn, j;
        mp_limb_t g3[3];

        if (hi <= klo)
        {
            flint_mpn_zero(out, outn);
            return;
        }

        wn = hi - klo;

        for (j = 0; j < 3; j++)
        {
            mp_size_t idx = klo - 3 + j;
            g3[j] = (idx < kl_len) ? kl[idx] : 0;
        }

        flint_mpn_mulmid(scratch, x, xn, y, yn, klo - 3, hi);   /* wn + 3 limbs */
        mpn_sub(scratch, scratch, wn + 3, g3, 3);
        if (scratch[2] != 0)
            mpn_add_1(scratch + 3, scratch + 3, wn, 1);

        flint_mpn_copyi(out, scratch + 3, wn);
        flint_mpn_zero(out + wn, outn - wn);
    }
    else
    {
        mp_size_t pl = FLINT_MIN(khi, pn);

        if (pl <= FLINT_MPN_MULHIGH_KNOWN_LOW_MULLOW_CUTOFF)
        {
            /* short product: use mullow on zero-padded operands, writing
               pl (possibly redundant) low limbs */
            mp_limb_t xp[FLINT_MPN_MULHIGH_KNOWN_LOW_MULLOW_CUTOFF];
            mp_limb_t yp[FLINT_MPN_MULHIGH_KNOWN_LOW_MULLOW_CUTOFF];
            mp_srcptr xx = x, yy = y;

            if (xn < pl)
            {
                flint_mpn_copyi(xp, x, xn);
                flint_mpn_zero(xp + xn, pl - xn);
                xx = xp;
            }
            if (yn < pl)
            {
                if (y == x && yn == xn)
                {
                    yy = xx;
                }
                else
                {
                    flint_mpn_copyi(yp, y, yn);
                    flint_mpn_zero(yp + yn, pl - yn);
                    yy = yp;
                }
            }

            flint_mpn_mullow_n(scratch, xx, yy, pl);
        }
        else
        {
            flint_mpn_mulmid(scratch, x, xn, y, yn, 0, pl);
        }

        if (pl < khi)
            flint_mpn_zero(scratch + pl, khi - pl);
        flint_mpn_copyi(out, scratch + klo, outn);
    }
}
