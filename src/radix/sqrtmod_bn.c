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
    Square root modulo B^n by one Karp-Markstein refinement on top of a
    half-precision reciprocal square root, the Hensel analogue of the last step
    of _padic_sqrt_p and the square-root counterpart of
    radix_divmod_bn_karp_markstein.

    With m = ceil(n/2):

        y = x^(-1/2) mod B^m              (half-precision reciprocal square root)
        b = x * y    mod B^m             (approximate root: b^2 == x mod B^m)
        s = b + y (x - b^2) / 2 mod B^n  (one Karp-Markstein step to n limbs)

    Correctness: b = x y, so b^2 = x (x y^2) == x (mod B^m), hence x - b^2 == 0
    (mod B^m). Writing x y^2 = 1 + B^m r and x - b^2 = B^m c,

        s^2 = b^2 + b y (x - b^2) + y^2 (x - b^2)^2 / 4
            = b^2 + (x y)(x - b^2) + ...
            == b^2 + (x - b^2) == x         (mod B^{2m}),

    using x y^2 == 1 (mod B^m) to drop the higher terms, and 2m >= n.

    Since x - b^2 == 0 (mod B^m), the correction y(x - b^2)/2 vanishes modulo
    B^m, so s shares its low m limbs with b and only its high n-m limbs are new.

    In both cases b is written straight into the low limbs of the root, the
    high part (b^2)[m', n') comes from a guarded middle product (its low m'
    limbs are x mod B^m') via the shared _radix_mulhigh_known_low, and only
    the high difference limbs are formed, so no full-width square, x copy or
    root copy is needed.

    Division by 2:
      - p odd: modular half of y*(x-b^2)[m, n) (the same short division by two as
        in radix_rsqrtmod_bn); the high limbs are limb-aligned.
      - p = 2: 1/2 does not exist, but x - b^2 == 0 (mod 2^{e m'}) is even
        (with m' = m + 1, see below), so the correction y (x - b^2) / 2 is an
        exact halving. Since only x - b^2 == 0 (mod 2^{e m'}) is guaranteed
        -- not (mod 2^{e m' + 1}) -- the halved correction is aligned at bit
        e m' - 1 rather than at limb m': it equals

            B^{m'-1} * ((y * dh) * 2^{e-1}  mod B^{n-m'+1}),

        dh = (x - b^2)[m', n+1), and is added over limbs [m'-1, n) of the
        root, straddling the limb boundary by one bit. The multiplication by
        2^{e-1} is a left shift by e-1 p-adic digits; the bits it pushes out
        at the top, like the original guard-limb formulation, fall above
        limb n and are discarded.

        The reciprocal root is taken one limb beyond ceil(n/2) (m' = m + 1)
        so that the few bits lost to the 2-adic squaring still leave all n
        limbs correct.

    s receives n limbs (caller provides room). Returns 1 on success, 0 if x is
    not a square modulo B (nothing written in that case). s must not alias x.
*/

/* W[0,k) = V[0,k) / 2 mod B^k for ODD radix; overflow-safe short division by 2,
   adding the (odd) modulus once when V is odd. W may alias V. */
static void
_radix_halve_odd(nn_ptr W, nn_srcptr V, slong k, const radix_t radix)
{
    ulong B = LIMB_RADIX(radix);
    ulong Bhalf = B >> 1;
    ulong rem, par;
    slong i;

    par = 0;
    for (i = 0; i < k; i++)
        par ^= V[i];
    par &= 1;

    rem = par;
    for (i = k - 1; i >= 0; i--)
    {
        ulong v = V[i];
        if (rem)
        {
            W[i] = Bhalf + (v >> 1) + (v & 1);
            rem = 1 - (v & 1);
        }
        else
        {
            W[i] = v >> 1;
            rem = v & 1;
        }
    }
}

int
radix_sqrtmod_bn(nn_ptr res, nn_srcptr x, slong xn, slong n, const radix_t radix)
{
    slong m, nm;
    nn_ptr y;
    TMP_INIT;

    FLINT_ASSERT(xn >= 1);
    FLINT_ASSERT(n >= 1);

    m = (n + 1) / 2;                 /* ceil(n/2) */
    nm = n - m;                      /* floor(n/2) */

    TMP_START;

    if (DIGIT_RADIX(radix) == 2)
    {
        /* p = 2: limb-aligned products, bit-straddling correction (see above) */
        slong e = radix->exp;
        slong mm = m + 1;            /* reciprocal-root precision (limbs) */
        slong ah;
        ulong bo;
        nn_ptr b;

        y = TMP_ALLOC(mm * sizeof(ulong));
        if (!radix_rsqrtmod_bn(y, x, xn, mm, radix))
        {
            TMP_END;
            return 0;
        }

        b = TMP_ALLOC(mm * sizeof(ulong));

        /* b = x y mod B^mm */
        radix_mulmid(b, x, FLINT_MIN(xn, mm), y, mm, 0, mm, radix);
        flint_mpn_copyi(res, b, FLINT_MIN(mm, n));

        /* nm = n - m = (n + 1) - mm is both the number of limbs of
           (x - b^2)[mm, n+1) and the size n - (mm - 1) of the correction
           window [mm-1, n) in the root. */
        if (nm > 0)
        {
            slong w = n + 1;
            nn_ptr b2h, dh, scr;

            b2h = TMP_ALLOC(nm * sizeof(ulong));
            dh  = TMP_ALLOC(nm * sizeof(ulong));
            scr = TMP_ALLOC(w * sizeof(ulong));

            /* (b^2)[mm, n+1): the low mm limbs of b^2 are x mod B^mm (known),
               so this is a guarded middle product rather than a full square. */
            _radix_mulhigh_known_low(b2h, b, mm, b, mm, x, FLINT_MIN(xn, mm),
                mm, w, scr, radix);

            /* dh = (x - b^2)[mm, n+1), where x is zero above limb xn. Rather
               than materialise the zero-padded operand, subtract over the
               limbs where x is present and negate over the zero limbs, then
               fold in the borrow out of the low subtraction (radix_neg takes
               no borrow-in). No incoming borrow at limb mm:
               x - b^2 == 0 (mod B^mm). */
            ah = (xn > mm) ? FLINT_MIN(xn - mm, nm) : 0;
            bo = 0;
            if (ah > 0)
                bo = radix_sub(dh, x + mm, ah, b2h, ah, radix);
            if (ah < nm)
            {
                radix_neg(dh + ah, b2h + ah, nm - ah, radix);
                if (bo)
                    radix_sub(dh + ah, dh + ah, nm - ah, &bo, 1, radix);
            }

            /* s[mm-1, n) += (y dh 2^{e-1})[0, nm); only the low nm <= mm
               limbs of y enter, digits shifted past the window top fall
               above limb n and are discarded, and the carry out of the
               addition is the reduction mod B^n. */
            radix_mulmid(scr, y, FLINT_MIN(mm, nm), dh, nm, 0, nm, radix);
            radix_lshift_digits(scr, scr, nm, e - 1, radix);

            if (n > mm)
                flint_mpn_zero(res + mm, n - mm);
            radix_add(res + (mm - 1), res + (mm - 1), nm, scr, nm, radix);
        }

        TMP_END;
        return 1;
    }
    else
    {
        /* odd p: limb-aligned Karp-Markstein (mirrors radix_divmod_bn_km) */
        y = TMP_ALLOC(m * sizeof(ulong));
        if (!radix_rsqrtmod_bn(y, x, xn, m, radix))
        {
            TMP_END;
            return 0;
        }

        /* b = x y mod B^m, written straight into the low m limbs of the root. */
        radix_mulmid(res, x, FLINT_MIN(xn, m), y, m, 0, m, radix);

        if (nm > 0)
        {
            nn_ptr b2h, dh, scr;
            slong ah;
            ulong bo;

            b2h = TMP_ALLOC(nm * sizeof(ulong));
            dh  = TMP_ALLOC(nm * sizeof(ulong));
            scr = TMP_ALLOC(n * sizeof(ulong));

            /* (b^2)[m, n): the low m limbs of b^2 are x mod B^m (known), so this
               is a guarded middle product rather than a full square. */
            _radix_mulhigh_known_low(b2h, res, m, res, m, x, FLINT_MIN(xn, m),
                m, n, scr, radix);

            /* dh = (x - b^2)[m, n) = x[m, n) - b2h, where x is zero above limb
               xn. Rather than materialise the zero-padded operand, subtract over
               the limbs where x is present and negate over the zero limbs, then
               fold in the borrow out of the low subtraction (radix_neg takes no
               borrow-in). No incoming borrow at limb m: x - b^2 == 0 (mod B^m). */
            ah = (xn > m) ? FLINT_MIN(xn - m, nm) : 0;
            bo = 0;
            if (ah > 0)
                bo = radix_sub(dh, x + m, ah, b2h, ah, radix);
            if (ah < nm)
            {
                radix_neg(dh + ah, b2h + ah, nm - ah, radix);
                if (bo)
                    radix_sub(dh + ah, dh + ah, nm - ah, &bo, 1, radix);
            }

            /* s[m, n) = (y * dh / 2)[0, nm) -> straight into the high limbs;
               only the low nm <= m limbs of y enter. */
            radix_mulmid(res + m, y, nm, dh, nm, 0, nm, radix);
            _radix_halve_odd(res + m, res + m, nm, radix);
        }

        TMP_END;
        return 1;
    }
}
