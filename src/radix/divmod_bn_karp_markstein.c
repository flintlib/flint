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
static void
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

/*
    Karp-Markstein Hensel division.

    Develops n limbs of the p-adic (Hensel) quotient q with

        q * b == a   (mod B^n),

    using one Karp-Markstein refinement step on top of a half-precision Hensel
    inverse of b. This is the Hensel analogue of radix_div_approx; it is
    asymptotically fast (cost dominated by an inverse and a few multiplications,
    O(M(n)) when bn is comparable to n), and so is preferable to the schoolbook
    radix_divmod_bn_classical when the divisor b is long.

    With m = ceil(n/2):

        y  = b^(-1) mod B^m                 (half-precision inverse)
        q0 = a * y mod B^m                  (low m limbs of q)
        d  = a - b*q0                       (== 0 mod B^m, since b*q0 == a mod B^m)
        qh = y * (d / B^m) mod B^{n-m}      (high n-m limbs of q)
        q  = q0 + B^m * qh

    Correctness: writing b*y = 1 + B^m*s, with q = q0 + y*d and d = a - b*q0,

        b*q = b*q0 + b*y*d = b*q0 + (1 + B^m s) d = a + B^m s d == a (mod B^{2m}),

    and 2m >= n, so q*b == a (mod B^n). The standard "compute b^(-1) to full
    precision then multiply by a" would instead need the inverse to n limbs; the
    Karp-Markstein trick halves that.

    Two of the products only need their high limbs, whose low limbs are known in
    advance, so they use _radix_mulhigh_known_low (an approximate high product):
      - b*q0: only limbs [m, n) are used and (b*q0) == a (mod B^m);
      - q*b (for the remainder): only limbs [n, n+bn) are used and q*b == a
        (mod B^n).

    q receives n quotient limbs (caller provides room for n) and may alias a.
    Returns 1 on success and 0 if b[0] is not invertible modulo B (nothing is
    written in that case).

    If rem != NULL it receives the resume (Henselian) remainder (a - q*b)/B^n
    reduced to its low bn limbs (room for bn required), unnormalised; then
    a == q*b + B^n*rem (mod B^{n+bn}). Exactness can be tested by the caller as
    in radix_divmod_bn_classical.
*/
int
radix_divmod_bn_karp_markstein(nn_ptr q, nn_ptr rem, nn_srcptr a, slong an,
    nn_srcptr b, slong bn, slong n, const radix_t radix)
{
    slong m, nm;
    nn_ptr y, qf;
    TMP_INIT;

    FLINT_ASSERT(an >= 1);
    FLINT_ASSERT(bn >= 1);
    FLINT_ASSERT(n >= 1);

    TMP_START;

    m = (n + 1) / 2;            /* ceil(n/2) */
    nm = n - m;                 /* floor(n/2), 0 <= nm <= m */

    /* y = b^(-1) mod B^m; 0 iff b[0] is a non-unit. */
    y = TMP_ALLOC(m * sizeof(ulong));
    if (!radix_invmod_bn(y, b, bn, m, radix))
    {
        TMP_END;
        return 0;
    }

    /* Build the quotient in scratch so that q may alias a. */
    qf = TMP_ALLOC(n * sizeof(ulong));

    /* q0 = a * y mod B^m -> low m limbs (a low product) */
    radix_mulmid(qf, a, FLINT_MIN(an, m), y, m, 0, m, radix);

    if (nm > 0)
    {
        nn_ptr bq0h, dh, scratch;
        slong ah;

        scratch = TMP_ALLOC(n * sizeof(ulong));

        /* bq0h = (b*q0)[m, n); the low m limbs of b*q0 equal a mod B^m. */
        bq0h = TMP_ALLOC(nm * sizeof(ulong));
        _radix_mulhigh_known_low(bq0h, b, bn, qf, m, a, an, m, n, scratch, radix);

        /* dh = (a - b*q0)[m, n) = a[m, n) - bq0h (no incoming borrow, since the
           low m limbs of a - b*q0 vanish). */
        dh = TMP_ALLOC(nm * sizeof(ulong));
        ah = (an > m) ? FLINT_MIN(an - m, nm) : 0;
        if (ah > 0)
            flint_mpn_copyi(dh, a + m, ah);
        flint_mpn_zero(dh + ah, nm - ah);
        radix_sub(dh, dh, nm, bq0h, nm, radix);

        /* qh = y * dh mod B^{n-m} -> high n-m limbs (a low product) */
        radix_mulmid(qf + m, y, m, dh, nm, 0, nm, radix);
    }

    if (rem != NULL)
    {
        nn_ptr prodh, scratch;
        slong ac;

        scratch = TMP_ALLOC((n + bn) * sizeof(ulong));

        /* prodh = (q*b)[n, n+bn); the low n limbs of q*b equal a mod B^n. */
        prodh = TMP_ALLOC(bn * sizeof(ulong));
        _radix_mulhigh_known_low(prodh, qf, n, b, bn, a, an, n, n + bn, scratch, radix);

        /* rem = (a - q*b)/B^n mod B^bn = a[n, n+bn) - prodh (no incoming borrow). */
        ac = (an > n) ? FLINT_MIN(an - n, bn) : 0;
        if (ac > 0)
            flint_mpn_copyi(rem, a + n, ac);
        flint_mpn_zero(rem + ac, bn - ac);
        radix_sub(rem, rem, bn, prodh, bn, radix);
    }

    flint_mpn_copyi(q, qf, n);

    TMP_END;
    return 1;
}
