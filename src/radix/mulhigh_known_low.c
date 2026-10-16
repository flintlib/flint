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
    out[0 .. khi-klo) = limbs [klo, khi) of the product x*y, given that the low
    klo limbs of x*y are already known and equal kl (with kl_len explicit limbs,
    the remaining low limbs being zero).

    When klo > 3, the product reaches khi (xn + yn >= khi) and the limb radix is
    large enough, the high limbs are obtained from an approximate (carry-less)
    middle product over [klo-3, khi) with 3 guard limbs, corrected using the
    known low limbs exactly as in radix_invmod_bn / _radix_divmod_bn_block:
    subtracting kl[klo-3..klo) gives guard limbs of true value zero, after which
    the high part is exact or 1 ulp too small, signalled by a nonzero top guard.
    Otherwise a full low product to khi is taken.

    scratch must have room for khi limbs (the fallback forms a low product up to
    khi limbs; the optimised branch uses only khi-klo+3 <= khi of them).
*/
void
_radix_mulhigh_known_low(nn_ptr out, nn_srcptr x, slong xn, nn_srcptr y, slong yn,
    nn_srcptr kl, slong kl_len, slong klo, slong khi, nn_ptr scratch,
    const radix_t radix)
{
    slong outn = khi - klo;

    if (klo > 3 && (xn + yn) >= khi
            && LIMB_RADIX(radix) >= (ulong) FLINT_MIN(xn, yn))
    {
        ulong g3[3], one = 1;
        slong j;

        for (j = 0; j < 3; j++)
        {
            slong idx = klo - 3 + j;
            g3[j] = (idx < kl_len) ? kl[idx] : 0;
        }

        radix_mulmid(scratch, x, xn, y, yn, klo - 3, khi, radix);  /* [klo-3, khi) */
        radix_sub(scratch, scratch, outn + 3, g3, 3, radix);       /* clean guards */
        if (scratch[2] != 0)
            radix_add(scratch + 3, scratch + 3, outn, &one, 1, radix);
        flint_mpn_copyi(out, scratch + 3, outn);
    }
    else
    {
        slong pl = FLINT_MIN(khi, xn + yn);
        radix_mulmid(scratch, x, xn, y, yn, 0, pl, radix);
        if (pl < khi)
            flint_mpn_zero(scratch + pl, khi - pl);
        flint_mpn_copyi(out, scratch + klo, outn);
    }
}

