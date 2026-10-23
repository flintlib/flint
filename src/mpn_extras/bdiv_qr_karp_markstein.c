/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

/*
    Karp-Markstein Hensel division, port of radix_divmod_bn_karp_markstein
    with B = 2^64. With m = ceil(n/2):

        y  = b^(-1) mod B^m
        q0 = a y mod B^m                   (low m limbs of q)
        d  = a - b q0                      (== 0 mod B^m)
        qh = y (d / B^m) mod B^(n-m)       (high n-m limbs of q)

    Writing b y = 1 + B^m s, b (q0 + y d) = a + B^m s d == a (mod B^(2m)).
    The high limbs of b q0 and of q b are obtained from products whose low
    limbs are known (a mod B^m resp. a mod B^n).

    q may alias a. If r != NULL it receives (a - q b) / B^n mod B^bn.
*/
void
flint_mpn_bdiv_qr_karp_markstein(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t an,
    mp_srcptr b, mp_size_t bn, mp_size_t n)
{
    mp_size_t m, nm;
    mp_ptr y, qf, scratch;
    TMP_INIT;

    FLINT_ASSERT(an >= 1);
    FLINT_ASSERT(bn >= 1);
    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(b[0] & 1);

    m = (n + 1) / 2;
    nm = n - m;

    TMP_START;

    y = TMP_ALLOC((m + n + (n + bn) + nm) * sizeof(mp_limb_t));
    qf = y + m;
    scratch = qf + n;

    flint_mpn_binv(y, b, bn, m);

    /* q0 = a y mod B^m */
    flint_mpn_mulmid(qf, a, FLINT_MIN(an, m), y, m, 0, m);

    if (nm > 0)
    {
        mp_ptr dh = scratch + (n + bn);
        mp_size_t ah;

        /* dh = (b q0)[m, n); the low m limbs of b q0 are a mod B^m */
        _flint_mpn_mulhigh_known_low(dh, b, bn, qf, m, a, an, m, n, scratch);

        /* dh = (a - b q0)[m, n), with a zero above limb an */
        mpn_neg(dh, dh, nm);
        ah = (an > m) ? FLINT_MIN(an - m, nm) : 0;
        if (ah > 0)
            mpn_add(dh, dh, nm, a + m, ah);

        /* qh = y dh mod B^(n-m); only the low nm <= m limbs of y enter */
        flint_mpn_mulmid(qf + m, y, nm, dh, nm, 0, nm);
    }

    if (r != NULL)
    {
        mp_size_t ac;

        /* r = (q b)[n, n+bn); the low n limbs of q b are a mod B^n */
        _flint_mpn_mulhigh_known_low(r, qf, n, b, bn, a, an, n, n + bn, scratch);

        /* r = (a - q b) / B^n mod B^bn */
        mpn_neg(r, r, bn);
        ac = (an > n) ? FLINT_MIN(an - n, bn) : 0;
        if (ac > 0)
            mpn_add(r, r, bn, a + n, ac);
    }

    flint_mpn_copyi(q, qf, n);

    TMP_END;
}
