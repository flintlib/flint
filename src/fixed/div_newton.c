/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "mpn_extras.h"
#include "fixed.h"

/* Approximate fixed-point inversion and division, ported from
   radix_inv_approx / radix_div_approx onto mpn arithmetic (the
   radix-B generality drops out: with B = 2^64 only the
   middle-product branches survive, division by two is a shift, and
   the binary basecases need no radix conversion).  The _newton suffix
   flags the contracts: not ulp-accurate, two-sided errors of a few
   ulps as documented in fixed.h.

   The Newton and Karp-Markstein iterations are carried verbatim;
   see radix/div.c for the error analysis (with the LIMB_RADIX(radix)
   large / mulhigh case applying throughout). */

void
fixed_inv_newton_basecase(nn_ptr Q, nn_srcptr A, slong An, slong n)
{
    nn_ptr U, R;
    nn_srcptr V;
    slong Un, Vn;
    TMP_INIT;

    Vn = FLINT_MIN(An, n + 1);
    Un = Vn + n + 1;

    TMP_START;
    U = TMP_ALLOC((Un + Vn) * sizeof(ulong));
    R = U + Un;
    V = A + An - Vn;

    flint_mpn_zero(U, Un - 1);
    U[Un - 1] = 1;

    if (Vn == 1)
        mpn_divrem_1(Q, 0, U, Un, V[0]);
    else
        mpn_tdiv_qr(Q, R, 0, U, Un, V, Vn);

    TMP_END;
}

/* Tuned on this machine (see dev/notes); must be at least 4 for the
   indexing/padding operations in the code to be valid. */
#ifndef FIXED_INV_NEWTON_CUTOFF
#define FIXED_INV_NEWTON_CUTOFF 14
#endif

void
fixed_inv_newton(nn_ptr Q, nn_srcptr A, slong An, slong n)
{
    if (n <= FIXED_INV_NEWTON_CUTOFF)
    {
        fixed_inv_newton_basecase(Q, A, An, n);
    }
    else
    {
        nn_ptr U, V, T, Vhigh, Uhigh;
        nn_srcptr Ahigh;
        slong m, Uhighn, Tn, Ahighn, Vhighn;
        int Unegative;
        TMP_INIT;

        /* T ~= 1/A with m fraction limbs, 2 integral limbs, stored in
           the high part of Q */
        m = (n + 1) / 2 + 1;
        fixed_inv_newton(Q + n - m, A, An, m);

        TMP_START;

        T = Q + n - m;
        Tn = m + 1 + (T[m + 1] != 0);

        flint_mpn_zero(Q, n - m);

        /* U ~= 1 - T*A' for the high n+1 limbs A' of A, as a
           fixed-point number with n fraction limbs: a middle product
           plus one control limb, the higher fraction limbs being 0 or
           B-1 a priori */
        slong control_limb = n / 2 + 2;
        slong A_low_zeroes = 0;

        if (An >= n + 1)
        {
            Ahighn = n + 1;
            Ahigh = A + An - Ahighn;
        }
        else
        {
            Ahighn = An;
            Ahigh = A;
            A_low_zeroes = n + 1 - An;
        }

        /* the full product has n+m+1 fraction limbs; remove m-1 to
           keep 2 guard limbs below the n retained ones */
        slong low_trunc = m - 1;

        if (A_low_zeroes <= low_trunc)
        {
            U = TMP_ALLOC((2 + (control_limb + 1)) * sizeof(ulong));
            flint_mpn_mulmid(U, Ahigh, Ahighn, T, Tn,
                low_trunc - A_low_zeroes,
                low_trunc - A_low_zeroes + 2 + (control_limb + 1));
            Uhigh = U + 2;
        }
        else
        {
            /* A is short: the window starts below the product of the
               zero-extended operand; fill the gap with zeros */
            slong gap = A_low_zeroes - low_trunc;
            U = TMP_ALLOC((2 + (control_limb + 1)) * sizeof(ulong));
            flint_mpn_mulmid(U + gap, Ahigh, Ahighn, T, Tn, 0,
                2 + (control_limb + 1) - gap);
            flint_mpn_zero(U, gap);
            Uhigh = U + 2;
        }

        Unegative = (Uhigh[control_limb] != 0);
        if (Unegative)
            mpn_neg(Uhigh, Uhigh, control_limb);

        Uhighn = control_limb;
        while (Uhighn > 0 && Uhigh[Uhighn - 1] == 0)
            Uhighn--;

        if (Uhighn != 0)
        {
            /* V = T * |1 - A'T| with n+2 fraction limbs */
            V = TMP_ALLOC((Tn + Uhighn - (m - 2)) * sizeof(ulong));
            flint_mpn_mulmid(V, T, Tn, Uhigh, Uhighn,
                m - 2, Tn + Uhighn);
            Vhigh = V + 2;
            Vhighn = Tn + Uhighn - (m - 2) - 2;

            /* T + T(1 - A'T) resp. T - T(A'T - 1) */
            if (Unegative)
                mpn_add(Q, Q, n + 2, Vhigh, Vhighn);
            else
                mpn_sub(Q, Q, n + 2, Vhigh, Vhighn);
        }

        TMP_END;
    }
}

void
fixed_div_newton_invmul(nn_ptr Q, nn_srcptr B, slong Bn, nn_srcptr A,
    slong An, slong n)
{
    nn_ptr T;
    TMP_INIT;
    TMP_START;

    if (Bn > n)
    {
        B = B + Bn - n;
        Bn = n;
    }

    T = TMP_ALLOC((Bn + n + 2) * sizeof(ulong));

    fixed_inv_newton(Q, A, An, n);
    /* should be a high multiplication */
    if (n + 2 >= Bn)
        flint_mpn_mul(T, Q, n + 2, B, Bn);
    else
        flint_mpn_mul(T, B, Bn, Q, n + 2);
    flint_mpn_copyi(Q, T + Bn, n + 2);

    TMP_END;
}

#ifndef FIXED_DIV_NEWTON_CUTOFF
#define FIXED_DIV_NEWTON_CUTOFF 12
#endif

/* Karp-Markstein division: evaluates TB + T * (B - A * TB) with
   T ~= 1/A carried at half precision */
void
fixed_div_newton(nn_ptr Q, nn_srcptr B, slong Bn, nn_srcptr A,
    slong An, slong n)
{
    nn_ptr U, V, T, Vhigh, Uhigh;
    nn_srcptr Ahigh;
    slong m, Uhighn, Tn, Ahighn, Vhighn;
    int Unegative;
    TMP_INIT;

    if (n <= FIXED_DIV_NEWTON_CUTOFF)
    {
        fixed_div_newton_invmul(Q, B, Bn, A, An, n);
        return;
    }

    m = (n + 1) / 2 + 1;

    fixed_inv_newton(Q + n - m, A, An, m);

    TMP_START;
    T = Q + n - m;
    Tn = m + 1 + (T[m + 1] != 0);

    flint_mpn_zero(Q, n - m);

    /* TB ~= B/A with m fraction limbs */
    nn_ptr TB, TBhigh;
    slong Bn2 = FLINT_MIN(m, Bn);
    nn_srcptr B2 = B + Bn - Bn2;

    if (Bn2 > 2)
    {
        TB = TMP_ALLOC((Tn + 2) * sizeof(ulong));
        flint_mpn_mulmid(TB, T, Tn, B2, Bn2, Bn2 - 2, Bn2 + Tn);
        TBhigh = TB + 2;
    }
    else
    {
        TB = TMP_ALLOC((Bn2 + Tn) * sizeof(ulong));
        if (Tn >= Bn2)
            flint_mpn_mul(TB, T, Tn, B2, Bn2);
        else
            flint_mpn_mul(TB, B2, Bn2, T, Tn);
        TBhigh = TB + Bn2;
    }

    slong control_limb = n / 2 + 2;
    slong A_low_zeroes = 0;

    if (An >= n + 1)
    {
        Ahighn = n + 1;
        Ahigh = A + An - Ahighn;
    }
    else
    {
        Ahighn = An;
        Ahigh = A;
        A_low_zeroes = n + 1 - An;
    }

    slong low_trunc = m - 1;

    if (A_low_zeroes <= low_trunc)
    {
        U = TMP_ALLOC((2 + (control_limb + 1)) * sizeof(ulong));
        flint_mpn_mulmid(U, Ahigh, Ahighn, TBhigh, Tn,
            low_trunc - A_low_zeroes,
            low_trunc - A_low_zeroes + 2 + (control_limb + 1));
        Uhigh = U + 2;
    }
    else
    {
        slong gap = A_low_zeroes - low_trunc;
        U = TMP_ALLOC((2 + (control_limb + 1)) * sizeof(ulong));
        flint_mpn_mulmid(U + gap, Ahigh, Ahighn, TBhigh, Tn, 0,
            2 + (control_limb + 1) - gap);
        flint_mpn_zero(U, gap);
        Uhigh = U + 2;
    }

    slong Un = control_limb + 1;

    /* (B - A*TB) with n fraction limbs, with sign */
    if (Bn > n - Un)
    {
        if (Bn >= n)
        {
            nn_srcptr Bhigh = B + Bn - n;
            mpn_sub_n(Uhigh, Bhigh, Uhigh, Un);
        }
        else
        {
            slong nz = n - Bn;
            ulong cy = mpn_neg(Uhigh, Uhigh, nz);
            mpn_sub_n(Uhigh + nz, B, Uhigh + nz, Un - nz);
            mpn_sub_1(Uhigh + nz, Uhigh + nz, Un - nz, cy);
        }

        Unegative = 1;
        if (Uhigh[Un - 1] >> (FLINT_BITS - 1))
        {
            mpn_neg(Uhigh, Uhigh, Un);
            Unegative = 0;
        }
    }
    else
    {
        Unegative = 0;
        if (Uhigh[Un - 1] >> (FLINT_BITS - 1))
        {
            mpn_neg(Uhigh, Uhigh, Un);
            Unegative = 1;
        }
    }

    Uhighn = Un;
    while (Uhighn > 0 && Uhigh[Uhighn - 1] == 0)
        Uhighn--;

    if (Uhighn != 0)
    {
        V = TMP_ALLOC((Tn + Uhighn - (m - 2)) * sizeof(ulong));
        flint_mpn_mulmid(V, T, Tn, Uhigh, Uhighn, m - 2, Tn + Uhighn);
        Vhigh = V + 2;
        Vhighn = Tn + Uhighn - (m - 2) - 2;

        /* overwrite T with TB */
        Q[n + 1] = 0;
        flint_mpn_copyi(Q + n - m, TBhigh, Tn);

        if (Unegative)
            mpn_add(Q, Q, n + 2, Vhigh, Vhighn);
        else
            mpn_sub(Q, Q, n + 2, Vhigh, Vhighn);
    }
    else
    {
        /* overwrite T with TB */
        Q[n + 1] = 0;
        flint_mpn_copyi(Q + n - m, TBhigh, Tn);
    }

    TMP_END;
}
