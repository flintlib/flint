/*
    Copyright (C) 2015 Kushagra Singh
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Opus 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
   Pollard rho with Brent's cycle detection, specialised to two-word
   moduli.  The inner loop uses Montgomery arithmetic with a fully
   unrolled two-limb REDC, which avoids the generic mpn call overhead
   and the normalisation/shift bookkeeping of
   flint_mpn_factor_pollard_brent_single.

   The iteration is y <- y^2 + c (mod n).  Montgomery representation is
   compatible with this: if yh = y*R then REDC(yh*yh) = y^2*R, and adding
   ch = c*R gives (y^2+c)*R.  The accumulated product of differences is
   scaled by a power of R, which is a unit mod n, so the gcd is unaffected.
*/

#include "longlong.h"
#include "mpn_extras.h"
#include "ulong_extras.h"

#include "ll_mont2.h"

#if FLINT_BITS == 64

/*
   Find a nontrivial factor of the two-word odd composite n = {nhi,nlo}.
   Returns 1 and writes the (one or two word) factor to factor[0..1] on
   success, 0 on failure.  max_iters bounds Brent's outer doubling
   parameter, so roughly 2*max_iters iterations of the map are performed
   per try.
*/
int
n_ll_factor_rho(nn_ptr factor, ulong nhi, ulong nlo, ulong max_tries,
                ulong max_iters)
{
    ulong n0 = nlo, n1 = nhi, np;
    ulong r2[2], nvec[2], gvec[2], qvec[2];
    ulong one1, one0;
    ulong try_i;
    mp_size_t gsize;

    if (n1 == 0 || (n0 & 1) == 0)
        return 0;

    np = LL_NEG_NINV(n0);
    nvec[0] = n0;
    nvec[1] = n1;

    /* R^2 mod n, where R = 2^128 */
    {
        ulong tmp[5], q[4];
        tmp[0] = tmp[1] = tmp[2] = tmp[3] = 0;
        tmp[4] = 1;                       /* 2^256 */
        mpn_tdiv_qr(q, r2, 0, tmp, 5, nvec, 2);
    }

    /* 1 in Montgomery form = R mod n */
    LL_MULMOD(one1, one0, 0, 1, r2[1], r2[0], n1, n0, np);

    for (try_i = 0; try_i < max_tries; try_i++)
    {
        ulong x1, x0, y1, y0, q1, q0, c1, c0;
        ulong r, k, i;

        /* c = try_i + 1, in Montgomery form */
        LL_MULMOD(c1, c0, 0, try_i + 1, r2[1], r2[0], n1, n0, np);

        y1 = one1; y0 = one0;
        q1 = one1; q0 = one0;
        x1 = y1;   x0 = y0;

        for (r = 1; r <= max_iters; r *= 2)
        {
            x1 = y1; x0 = y0;

            for (i = 0; i < r; i++)
            {
                LL_SQRMOD(y1, y0, y1, y0, n1, n0, np);
                LL_ADDMOD(y1, y0, y1, y0, c1, c0, n1, n0);
            }

            for (k = 0; k < r; k += 128)
            {
                ulong lim = FLINT_MIN(128, r - k);

                for (i = 0; i < lim; i++)
                {
                    ulong d1, d0;

                    LL_SQRMOD(y1, y0, y1, y0, n1, n0, np);
                    LL_ADDMOD(y1, y0, y1, y0, c1, c0, n1, n0);

                    if (x1 > y1 || (x1 == y1 && x0 >= y0))
                        sub_ddmmss(d1, d0, x1, x0, y1, y0);
                    else
                        sub_ddmmss(d1, d0, y1, y0, x1, x0);

                    /* d == 0 means the cycle closed with no factor split */
                    if ((d1 | d0) == 0)
                        goto next_try;

                    LL_MULMOD(q1, q0, q1, q0, d1, d0, n1, n0, np);
                }

                if ((q1 | q0) == 0)
                    goto next_try;

                qvec[0] = q0;
                qvec[1] = q1;
                gsize = flint_mpn_gcd_full(gvec, nvec, 2, qvec, q1 ? 2 : 1);

                if (gsize == 1)
                {
                    if (gvec[0] != 1)
                    {
                        factor[0] = gvec[0];
                        factor[1] = 0;
                        return 1;
                    }
                }
                else
                {
                    if (!(gvec[0] == n0 && gvec[1] == n1))
                    {
                        factor[0] = gvec[0];
                        factor[1] = gvec[1];
                        return 1;
                    }
                    /* gcd == n: the batch swallowed all factors, retry */
                    goto next_try;
                }
            }
        }
next_try: ;
    }

    return 0;
}

#endif
