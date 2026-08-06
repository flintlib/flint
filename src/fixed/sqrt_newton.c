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

/* Approximate fixed-point reciprocal square roots and square roots,
   ported from radix_rsqrt_1_approx / radix_rsqrt_approx /
   radix_sqrt_approx onto mpn arithmetic; see div_newton.c for the
   porting conventions and radix/sqrt.c for the error analyses. */

void
fixed_rsqrt_ui_newton_basecase(nn_ptr Y, ulong a, slong n)
{
    nn_ptr U;
    slong nU;
    TMP_INIT;

    TMP_START;
    U = TMP_ALLOC((2 * n + 1) * sizeof(ulong));

    /* floor(sqrt(a B^(2n))) = floor(B^n sqrt(a)), then divide by a */
    flint_mpn_zero(U, 2 * n);
    U[2 * n] = a;

    mpn_sqrtrem(U, NULL, U, 2 * n + 1);
    nU = n + 1;
    while (nU > 0 && U[nU - 1] == 0)
        nU--;
    mpn_divrem_1(U, 0, U, nU, a);
    while (nU > 0 && U[nU - 1] == 0)
        nU--;

    FLINT_ASSERT(nU <= n);
    flint_mpn_copyi(Y, U, nU);
    flint_mpn_zero(Y + nU, n - nU);

    TMP_END;
}

#ifndef FIXED_RSQRT_UI_NEWTON_CUTOFF
#define FIXED_RSQRT_UI_NEWTON_CUTOFF 24
#endif

void
fixed_rsqrt_ui_newton(nn_ptr Y, ulong a, slong n)
{
    if (n <= FIXED_RSQRT_UI_NEWTON_CUTOFF)
    {
        fixed_rsqrt_ui_newton_basecase(Y, a, n);
    }
    else
    {
        slong m, Un;
        ulong cy;
        nn_ptr T, U;
        nn_ptr Yhi, Thi;
        TMP_INIT;

        m = (n + 1) / 2 + 1;

        Yhi = Y + n - m;
        fixed_rsqrt_ui_newton(Yhi, a, m);
        flint_mpn_zero(Y, n - m);

        TMP_START;
        T = TMP_ALLOC((m + 3 + m + 2) * sizeof(ulong));
        U = T + m + 3;

        /* T = low m+2 window of Yhi^2 (2m fraction limbs) */
        flint_mpn_mulmid(T, Yhi, m, Yhi, m, 0, m + 2);
        cy = mpn_mul_1(U, T, m + 2, a);   /* mod B^(m+2) */
        (void) cy;
        cy = (U[m + 1] == 0);
        /* a*y^2 - 1 or 1 - a*y^2 */
        if (!cy)
            mpn_neg(U, U, m + 2);

        mpn_rshift(U, U, m + 2, 1);
        Un = m + 2;
        while (Un > 0 && U[Un - 1] == 0)
            Un--;
        if (Un + n - 2 * m <= 0)
            goto cleanup;

        flint_mpn_mulmid(T, Yhi, m, U, Un, m - 2, Un + m);
        Thi = T + 2;

        if (cy)
            mpn_sub(Y, Y, n, Thi + 2 * m - n, Un + n - 2 * m);
        else
            mpn_add(Y, Y, n, Thi + 2 * m - n, Un + n - 2 * m);

cleanup:
        TMP_END;
    }
}

void
fixed_rsqrt_newton_basecase(nn_ptr Q, nn_srcptr A, slong An, slong n)
{
    nn_ptr W, U, bq, br;
    nn_srcptr V;
    slong Vn, g, Wn, Un, bsn, bqn;
    slong guard_limbs = 2;
    TMP_INIT;

    Vn = FLINT_MIN(An, n + guard_limbs);
    g = FLINT_MAX(0, n + guard_limbs - Vn);
    V = A + An - Vn;

    Wn = 2 * (Vn + g);
    Un = n + Vn + g + 1;

    TMP_START;
    W = TMP_ALLOC((Wn + Un + (Wn + 3) / 2 + Un) * sizeof(ulong));
    U = W + Wn;
    bq = U + Un;
    br = bq + (Wn + 3) / 2;

    /* W = V B^(Vn + 2g); its top limb may be zero when A < 1/B */
    flint_mpn_zero(W, Vn + 2 * g);
    flint_mpn_copyi(W + Vn + 2 * g, V, Vn);

    while (Wn > 0 && W[Wn - 1] == 0)
        Wn--;
    FLINT_ASSERT(Wn >= 2 * (Vn + g) - 1);

    flint_mpn_zero(U, Un - 1);
    U[Un - 1] = 1;

    /* bq = floor(B^(Un-1) / floor(sqrt(W))); the quotient has at most
       n+2 limbs since sqrt(W) >= B^(Vn+g-1) */
    bsn = (Wn + 1) / 2;
    mpn_sqrtrem(br, NULL, W, Wn);       /* sqrt into br as scratch */
    bqn = Un - bsn + 1;
    FLINT_ASSERT(bqn <= n + 2);
    mpn_tdiv_qr(bq, W, 0, U, Un, br, bsn);   /* remainder into W */

    while (bqn > 0 && bq[bqn - 1] == 0)
        bqn--;
    flint_mpn_copyi(Q, bq, bqn);
    flint_mpn_zero(Q + bqn, (n + 2) - bqn);

    TMP_END;
}

/* Tuned on this machine (see dev/notes); must be at least 5 so that
   m - 2 >= 2 and 2m - n - 2 >= 0 hold in the recursive cases of
   fixed_rsqrt_newton and fixed_sqrt_newton. */
#ifndef FIXED_RSQRT_NEWTON_CUTOFF
#define FIXED_RSQRT_NEWTON_CUTOFF 20
#endif

void
fixed_rsqrt_newton(nn_ptr Q, nn_srcptr A, slong An, slong n)
{
    if (n <= FIXED_RSQRT_NEWTON_CUTOFF)
    {
        fixed_rsqrt_newton_basecase(Q, A, An, n);
    }
    else
    {
        nn_ptr U, V, W, T, Vhigh, Uhigh;
        nn_srcptr Ahigh;
        slong m, Uhighn, Tn, Ahighn, Vhighn;
        int Unegative;
        TMP_INIT;

        /* T ~= 1/sqrt(A) with m fraction limbs, 2 integral limbs, in
           the high part of Q */
        m = (n + 1) / 2 + 1;
        fixed_rsqrt_newton(Q + n - m, A, An, m);

        TMP_START;

        T = Q + n - m;
        Tn = m + 1 + (T[m + 1] != 0);

        flint_mpn_zero(Q, n - m);

        /* W = T^2 with 2m fraction limbs (full product; nearly all
           its limbs are consumed below) */
        W = TMP_ALLOC(2 * Tn * sizeof(ulong));
        flint_mpn_sqr(W, T, Tn);

        /* U ~= 1 - A'*T^2 with n fraction limbs, via the middle
           product + control limb construction of fixed_inv_newton */
        slong control_limb = n / 2 + 2;
        slong A_low_zeroes = 0;

        if (An >= n + 2)
        {
            Ahighn = n + 2;
            Ahigh = A + An - Ahighn;
        }
        else
        {
            Ahighn = An;
            Ahigh = A;
            A_low_zeroes = n + 2 - An;
        }

        /* full product: n+2+2m fraction limbs; keep 2 guard limbs
           under the n retained ones */
        slong low_trunc = 2 * m;

        /* A_low_zeroes <= n + 1 <= 2m - 1, so the band always starts
           at limb >= 1 */
        FLINT_ASSERT(A_low_zeroes <= low_trunc - 1);

        U = TMP_ALLOC((2 + (control_limb + 1)) * sizeof(ulong));
        flint_mpn_mulmid(U, Ahigh, Ahighn, W, 2 * Tn,
            low_trunc - A_low_zeroes,
            low_trunc - A_low_zeroes + 2 + (control_limb + 1));
        Uhigh = U + 2;

        Unegative = (Uhigh[control_limb] != 0);
        if (Unegative)
            mpn_neg(Uhigh, Uhigh, control_limb);

        Uhighn = control_limb;
        while (Uhighn > 0 && Uhigh[Uhighn - 1] == 0)
            Uhighn--;

        if (Uhighn != 0)
        {
            /* correction T * (1 - A T^2) / 2 */
            mpn_rshift(Uhigh, Uhigh, Uhighn, 1);
            while (Uhighn > 0 && Uhigh[Uhighn - 1] == 0)
                Uhighn--;
        }

        if (Uhighn != 0)
        {
            V = TMP_ALLOC((Tn + Uhighn - (m - 2)) * sizeof(ulong));
            flint_mpn_mulmid(V, T, Tn, Uhigh, Uhighn,
                m - 2, Tn + Uhighn);
            Vhigh = V + 2;
            Vhighn = Tn + Uhighn - (m - 2) - 2;

            if (Unegative)
                mpn_add(Q, Q, n + 2, Vhigh, Vhighn);
            else
                mpn_sub(Q, Q, n + 2, Vhigh, Vhighn);
        }

        TMP_END;
    }
}

void
fixed_sqrt_newton_rsqrtmul(nn_ptr Q, nn_srcptr A, slong An, slong n)
{
    nn_ptr T;
    nn_srcptr A2;
    slong A2n;
    TMP_INIT;
    TMP_START;

    A2n = FLINT_MIN(An, n + 2);
    A2 = A + An - A2n;

    T = TMP_ALLOC((A2n + n + 2) * sizeof(ulong));

    fixed_rsqrt_newton(Q, A, An, n);
    /* should be a high multiplication */
    if (n + 2 >= A2n)
        flint_mpn_mul(T, Q, n + 2, A2, A2n);
    else
        flint_mpn_mul(T, A2, A2n, Q, n + 2);
    flint_mpn_copyi(Q, T + A2n, n + 2);

    TMP_END;
}

#ifndef FIXED_SQRT_NEWTON_CUTOFF
#define FIXED_SQRT_NEWTON_CUTOFF 12
#endif

/* Karp-Markstein square root: S + T * (A - S^2) / 2 with S = T * A,
   T ~= 1/sqrt(A) carried at half precision */
void
fixed_sqrt_newton(nn_ptr Q, nn_srcptr A, slong An, slong n)
{
    nn_ptr U, V, T, TS, Vhigh, Uhigh, TShigh;
    slong m, Un, Uhighn, Tn, Vhighn;
    int Unegative;
    TMP_INIT;

    if (n <= FIXED_SQRT_NEWTON_CUTOFF)
    {
        fixed_sqrt_newton_rsqrtmul(Q, A, An, n);
        return;
    }

    m = (n + 1) / 2 + 1;

    fixed_rsqrt_newton(Q + n - m, A, An, m);

    TMP_START;
    T = Q + n - m;
    Tn = m + 1 + (T[m + 1] != 0);

    flint_mpn_zero(Q, n - m);

    /* S = T * A2 ~= sqrt(A) with m fraction limbs; two extra limbs of
       the truncated operand are needed, the truncation error being
       multiplied by T <= B */
    slong A2n = FLINT_MIN(m + 2, An);
    nn_srcptr A2 = A + An - A2n;

    if (A2n > 2)
    {
        TS = TMP_ALLOC((Tn + 2) * sizeof(ulong));
        flint_mpn_mulmid(TS, T, Tn, A2, A2n, A2n - 2, A2n + Tn);
        TShigh = TS + 2;
    }
    else
    {
        TS = TMP_ALLOC((A2n + Tn) * sizeof(ulong));
        if (Tn >= A2n)
            flint_mpn_mul(TS, T, Tn, A2, A2n);
        else
            flint_mpn_mul(TS, A2, A2n, T, Tn);
        TShigh = TS + A2n;
    }

    slong control_limb = n / 2 + 2;

    Un = control_limb + 1;

    /* U = band of S^2 at weights B^-n .. B^(control_limb - n), plus 2
       guard limbs; the full square has 2m fraction limbs, so the band
       starts 2m - n limbs up */
    FLINT_ASSERT(2 * m - n - 2 >= 0);

    U = TMP_ALLOC((2 + Un) * sizeof(ulong));
    flint_mpn_mulmid(U, TShigh, Tn, TShigh, Tn,
        2 * m - n - 2, 2 * m - n + Un);
    Uhigh = U + 2;

    /* (A - S^2) with n fraction limbs, with sign */
    if (An > n - Un)
    {
        if (An >= n)
        {
            nn_srcptr Ahigh = A + An - n;
            mpn_sub_n(Uhigh, Ahigh, Uhigh, Un);
        }
        else
        {
            slong nz = n - An;
            ulong cy = mpn_neg(Uhigh, Uhigh, nz);
            mpn_sub_n(Uhigh + nz, A, Uhigh + nz, Un - nz);
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
        /* correction T * (A - S^2) / 2 */
        mpn_rshift(Uhigh, Uhigh, Uhighn, 1);
        while (Uhighn > 0 && Uhigh[Uhighn - 1] == 0)
            Uhighn--;
    }

    if (Uhighn != 0)
    {
        V = TMP_ALLOC((Tn + Uhighn - (m - 2)) * sizeof(ulong));
        flint_mpn_mulmid(V, T, Tn, Uhigh, Uhighn, m - 2, Tn + Uhighn);
        Vhigh = V + 2;
        Vhighn = Tn + Uhighn - (m - 2) - 2;

        /* overwrite T with S */
        Q[n + 1] = 0;
        flint_mpn_copyi(Q + n - m, TShigh, Tn);

        if (Unegative)
            mpn_add(Q, Q, n + 2, Vhigh, Vhighn);
        else
            mpn_sub(Q, Q, n + 2, Vhigh, Vhighn);
    }
    else
    {
        Q[n + 1] = 0;
        flint_mpn_copyi(Q + n - m, TShigh, Tn);
    }

    TMP_END;
}
