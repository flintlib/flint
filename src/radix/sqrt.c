/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "radix.h"

void radix_rsqrt_1_approx_basecase(nn_ptr Y, ulong a, slong n, const radix_t radix)
{
    nn_ptr U, V;
    slong nU, nV;
    TMP_INIT;

    TMP_START;
    U = TMP_ALLOC((4 * n + 2) * sizeof(ulong));
    V = U + 2 * n + 1;

    flint_mpn_zero(U, 2 * n);
    U[2 * n] = a;

    nV = radix_get_mpn(V, U, 2 * n + 1, radix);
    mpn_sqrtrem(U, NULL, V, nV);
    nU = (nV + 1) / 2;
    MPN_NORM(U, nU);
    mpn_divrem_1(U, 0, U, nU, a);
    MPN_NORM(U, nU);

    nV = radix_set_mpn(V, U, nU, radix);

    flint_mpn_copyi(Y, V, nV);
    flint_mpn_zero(Y + nV, n - nV);

    TMP_END;
}

// todo: allocs
// todo: fast mul_1
// todo: fast divrem_two

/*
The error is smaller than 2 ulp, and for large radix, 1+eps ulp. Proof sketch:
suppose on input y = 1/sqrt(a) + eps where we assume that |eps| <= C*B^(-m)
for some C < 2. The mathematical error after the Newton iteration is

    1/sqrt(a) - 3*sqrt(a)/2 * eps^2 - a/2 * eps^3.

The calculation of 1 - a*y^2 is done exactly and generates a fixed-point number
with 2m-limb precision. Subsequently dividing by 2 generates up to B^(-2m)
error. Rounding the result of the final multiplication by y generates up to
B^(-n) error, plus B^(-2m) error if mulhigh is used.

Note that m >= n/2 + 1. Setting e.g. C = 1.7 one can check explicitly that
the error is bounded by C*B^(-n) in the worst case B = 3.
For larger B one can choose C closer to 1.
*/

void radix_rsqrt_1_approx(nn_ptr Y, ulong a, slong n, const radix_t radix)
{
    if (n <= 4)
    {
        radix_rsqrt_1_approx_basecase(Y, a, n, radix);
    }
    else
    {
        slong m, Un;
        ulong cy;
        nn_ptr T, U;
        slong Talloc, Ualloc;
        nn_ptr Yhi, Thi;
        TMP_INIT;

        m = (n + 1) / 2 + 1;

        Yhi = Y + n - m;
        radix_rsqrt_1_approx(Yhi, a, m, radix);
        flint_mpn_zero(Y, n - m);

        /* TODO: reduce temporary space */
        if (LIMB_RADIX(radix) > m)
            Talloc = m + 3;   /* In case of high product */
        else
            Talloc = 2 * m + 2;   /* In case of full product */
        Ualloc = m + 2;

        TMP_START;
        T = TMP_ALLOC((Talloc + Ualloc) * sizeof(ulong));
        U = T + Talloc;

        /* T = Yhi^2 with 2m fraction limbs */
        radix_mulmid(T, Yhi, m, Yhi, m, 0, m + 2, radix);
        radix_mulmid(U, T, m + 2, &a, 1, 0, m + 2, radix);
        cy = (U[m + 1] == 0);
        /* a*y^2-1 or 1-a*y^2 */
        if (!cy)
            radix_neg(U, U, m + 2, radix);

        radix_divrem_two(U, U, m + 2, radix);
        Un = m + 2;
        MPN_NORM(U, Un);
        if (Un == 0)
            goto cleanup;

        if (LIMB_RADIX(radix) > m)
        {
            radix_mulmid(T, Yhi, m, U, Un, m - 2, Un + m, radix);
            Thi = T + 2;
        }
        else
        {
            radix_mulmid(T, Yhi, m, U, Un, 0, Un + m, radix);
            Thi = T + m;
        }

        if (cy)
            radix_sub(Y, Y, n, Thi + 2 * m - n, Un + n - 2 * m, radix);
        else
            radix_add(Y, Y, n, Thi + 2 * m - n, Un + n - 2 * m, radix);

cleanup:
        TMP_END;
    }
}

