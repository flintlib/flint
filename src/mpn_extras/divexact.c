/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

/* Hensel division beats GMP's mpn_divexact from about 1024 limbs for
   balanced operands and earlier for long quotients (profile/p-divexact.c) */
#ifndef FLINT_MPN_DIVEXACT_UNBALANCED_CUTOFF   /* n >= 4 bn */
#define FLINT_MPN_DIVEXACT_UNBALANCED_CUTOFF 64
#endif

/*
    Exact division by Hensel (2-adic) division: with b = 2^v B^k b' for odd
    b', a = 2^v B^k a' and q = a' / b' = a' b'^(-1) mod B^n where n is the
    number of quotient limbs. Requires that b divides a.
*/
void
_flint_mpn_divexact_hensel(mp_ptr q, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
{
    mp_size_t k = 0, n;
    unsigned int v;
    mp_ptr as, bs;
    int truncated = 0;
    TMP_INIT;

    FLINT_ASSERT(an >= bn);
    FLINT_ASSERT(bn >= 1);
    FLINT_ASSERT(b[bn - 1] != 0);

    n = an - bn + 1;

    /* the quotient has fewer limbs when the top limb of a is below that of b */
    if (n > 1 && a[an - 1] < b[bn - 1])
    {
        q[n - 1] = 0;
        n--;
    }

    /* strip zero low limbs */
    while (b[k] == 0)
        k++;
    FLINT_ASSERT(flint_mpn_zero_p(a, k));
    a += k;
    an -= k;
    b += k;
    bn -= k;

    v = flint_ctz(b[0]);

    /* only the low n + bn limbs of a enter a Hensel division */
    if (an > n + bn + 1)
    {
        an = n + bn + 1;
        truncated = 1;
    }

    if (v == 0)
    {
        flint_mpn_bdiv_q(q, a, an, b, bn, n);
        return;
    }

    TMP_START;
    as = TMP_ALLOC((an + bn) * sizeof(mp_limb_t));
    bs = as + an;
    mpn_rshift(as, a, an, v);
    mpn_rshift(bs, b, bn, v);
    /* the top shifted limb of a is incomplete when a was truncated above;
       it is never needed */
    flint_mpn_bdiv_q(q, as, an - truncated, bs, bn, n);
    TMP_END;
}

/*
    A bidirectional variant (low half of the quotient by Hensel division,
    high half by Euclidean division, so that two half-size problems replace
    one full-size one) was implemented and measured: it was about 1.2x
    slower than the pure Hensel division for balanced operands and 2.5x
    slower when the quotient is shorter than the divisor, since the
    Euclidean half (an inverse plus three half-size products) costs more
    than the Hensel half (an inverse plus two). Exact division is therefore
    purely 2-adic.
*/

void
_flint_mpn_divexact(mp_ptr q, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
{
    FLINT_ASSERT(an >= bn);
    FLINT_ASSERT(bn >= 1);
    FLINT_ASSERT(b[bn - 1] != 0);

    if (bn == 1)
    {
        mpn_divexact_1(q, a, an, b[0]);
        return;
    }

#if FLINT_HAVE_NATIVE_mpn_divexact
    {
        mp_size_t n = an - bn + 1;

        if (!(FLINT_MIN(bn, n) >= FLINT_MPN_DIVEXACT_NEWTON_CUTOFF
            || (bn >= FLINT_MPN_DIVEXACT_UNBALANCED_CUTOFF && n >= 4 * bn)))
        {
            mpn_divexact(q, a, an, b, bn);
            return;
        }
    }
#endif

    _flint_mpn_divexact_hensel(q, a, an, b, bn);
}
