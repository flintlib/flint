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

    flint_mpn_sqrtrem(U, NULL, U, 2 * n + 1);
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
    flint_mpn_sqrtrem(br, NULL, W, Wn);       /* sqrt into br as scratch */
    bqn = Un - bsn + 1;
    FLINT_ASSERT(bqn <= n + 2);
    flint_mpn_tdiv_qr(bq, W, U, Un, br, bsn);   /* remainder into W */

    while (bqn > 0 && bq[bqn - 1] == 0)
        bqn--;
    flint_mpn_copyi(Q, bq, bqn);
    flint_mpn_zero(Q + bqn, (n + 2) - bqn);

    TMP_END;
}

/* Reciprocal square roots by third-order Newton steps: each step starts
   from m ~ n/3 limbs and applies

       1/sqrt(a) = y (1 + v)^(-1/2) = y (1 - v/2 + 3v^2/8 - ...),

   the tail (5/16)|v|^3 y being below the rounding for 3m >= n + 3.
   The residual v = a y^2 - 1 is formed as a band of a y^2 with the unit
   wrapped into a control limb (the construction of fixed_inv_newton);
   v^2 needs only the top ~n - 2m limbs, and the whole correction
   w = (v -+ 3v^2/4)/2 is formed in units two limbs below the output
   before those guard limbs are dropped, so that it is rounded to one
   unit of the output, plus a fraction B^-2 for v^2.  A step costs the
   square T^2 of the start (m limbs), the residual band A T^2, a square
   of about n/3 limbs for v^2 and the product T w.  Against a
   second-order step (the same square, band and product at m = n/2,
   covering half the precision instead of a third) this is 8-34%
   faster from 32 limbs on. */

#ifndef FIXED_RSQRT_UI_NEWTON_CUTOFF
#define FIXED_RSQRT_UI_NEWTON_CUTOFF 24
#endif

void
fixed_rsqrt_ui_newton(nn_ptr Y, ulong a, slong n)
{
    slong m, V0n, wn, sqn, sh, i;
    int pos;
    nn_ptr T, V0, W, SQ, Yhi, Vp;
    TMP_INIT;

    if (n <= FIXED_RSQRT_UI_NEWTON_CUTOFF)
    {
        fixed_rsqrt_ui_newton_basecase(Y, a, n);
        return;
    }

    /* 3m >= n + 5: the tail 20 B^(3/2 - 3m) is far below B^-n */
    m = (n + 5 + 2) / 3;

    Yhi = Y + n - m;
    fixed_rsqrt_ui_newton(Yhi, a, m);
    flint_mpn_zero(Y, n - m);

    TMP_START;
    T = TMP_ALLOC((m + 2) * sizeof(ulong));
    V0 = TMP_ALLOC((m + 2) * sizeof(ulong));

    /* v B^(2m) = a Yhi^2 - B^(2m) exactly, from the low m + 2 limbs of
       Yhi^2 (|v| < 4 B^(1/2 - m), so |v B^(2m)| < B^(m+1)) */
    flint_mpn_mulmid(T, Yhi, m, Yhi, m, 0, m + 2);
    mpn_mul_1(V0, T, m + 2, a);
    pos = (V0[m + 1] == 0);     /* v >= 0: y too large */
    if (!pos)
        mpn_neg(V0, V0, m + 2);
    V0n = m + 2;
    while (V0n > 0 && V0[V0n - 1] == 0)
        V0n--;
    if (V0n == 0)
        goto cleanup;

    /* w (units B^-(n+2)) = (|v| -+ 3 v^2/4) / 2 with
       |v| = V0 B^(n + 2 - 2m) and v^2 = V0^2 B^(n + 2 - 4m) */
    sh = n + 2 - 2 * m;
    wn = sh + V0n + 1;
    W = TMP_ALLOC(wn * sizeof(ulong));
    flint_mpn_zero(W, sh);
    flint_mpn_copyi(W + sh, V0, V0n);
    W[wn - 1] = 0;

    SQ = TMP_ALLOC((2 * V0n + 1) * sizeof(ulong));
    flint_mpn_sqr(SQ, V0, V0n);
    sqn = 2 * V0n;
    i = 4 * m - n - 2;          /* limbs of V0^2 below the unit */
    if (sqn > i)
    {
        nn_ptr sq = SQ + i;
        slong l = sqn - i;
        ulong cy;
        /* 3 v^2 / 4 */
        cy = mpn_mul_1(sq, sq, l, 3);
        sq[l] = cy;
        l++;
        mpn_rshift(sq, sq, l, 2);
        while (l > 0 && sq[l - 1] == 0)
            l--;
        if (l > 0)
        {
            if (pos)
                mpn_sub(W, W, wn, sq, l);
            else
                mpn_add(W, W, wn, sq, l);
        }
    }
    mpn_rshift(W, W, wn, 1);
    while (wn > 0 && W[wn - 1] == 0)
        wn--;
    if (wn == 0)
        goto cleanup;

    /* Y -+= Yhi w: limb i of the product has weight B^(i - m - n - 2);
       from index m (two guard limbs), Vp + 2 is at B^-n */
    Vp = TMP_ALLOC((wn + 2) * sizeof(ulong));
    flint_mpn_mulmid(Vp, Yhi, m, W, wn, m, m + wn);
    {
        slong vn = wn - 2;
        if (vn > 0)
        {
            if (pos)
                mpn_sub(Y, Y, n, Vp + 2, FLINT_MIN(vn, n));
            else
                mpn_add(Y, Y, n, Vp + 2, FLINT_MIN(vn, n));
        }
    }

cleanup:
    TMP_END;
}

#ifndef FIXED_RSQRT_NEWTON_CUTOFF
#define FIXED_RSQRT_NEWTON_CUTOFF 20
#endif

void
fixed_rsqrt_newton(nn_ptr Q, nn_srcptr A, slong An, slong n)
{
    nn_ptr U, V, W, T, Ug, SQ, Vhigh;
    nn_srcptr Ahigh;
    slong m, Tn, Ahighn, Ugn, Vhighn, control_limb, low_trunc, sqn, s, hn;
    slong A_low_zeroes = 0;
    int Unegative;
    TMP_INIT;

    if (n <= FIXED_RSQRT_NEWTON_CUTOFF)
    {
        fixed_rsqrt_newton_basecase(Q, A, An, n);
        return;
    }

    /* 3m >= n + 4: |u| < 8 B^-m, tail 160 B^-3m T far below B^-n T */
    m = (n + 4 + 2) / 3;
    fixed_rsqrt_newton(Q + n - m, A, An, m);

    TMP_START;

    T = Q + n - m;
    Tn = m + 1 + (T[m + 1] != 0);
    flint_mpn_zero(Q, n - m);

    /* W = T^2 with 2m fraction limbs */
    W = TMP_ALLOC(2 * Tn * sizeof(ulong));
    flint_mpn_sqr(W, T, Tn);

    /* the band of A T^2 at weights B^-(n+2) .. B^-(m-2), with the unit
       wrapped into a control limb */
    control_limb = n - m + 2;
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
    low_trunc = 2 * m;

    /* a short A can start the band below the product's bottom (2m ~
       2n/3 can be below the n + 2 - An missing limbs): those band limbs
       are zero */
    U = TMP_ALLOC((2 + (control_limb + 1)) * sizeof(ulong));
    {
        slong lo = low_trunc - A_low_zeroes;
        slong hi = lo + 2 + (control_limb + 1);
        if (lo >= 0)
            flint_mpn_mulmid(U, Ahigh, Ahighn, W, 2 * Tn, lo, hi);
        else
        {
            flint_mpn_zero(U, FLINT_MIN(-lo, hi - lo));
            if (hi > 0)
                flint_mpn_mulmid(U - lo, Ahigh, Ahighn, W, 2 * Tn, 0, hi);
        }
    }

    /* |u| with its two guard limbs, as an integer in units B^-(n+2);
       Unegative: A T^2 < 1, u = 1 - A T^2 > 0, T too small */
    Unegative = (U[2 + control_limb] != 0);
    Ug = U;
    Ugn = control_limb + 2;
    if (Unegative)
        mpn_neg(Ug, Ug, Ugn);
    while (Ugn > 0 && Ug[Ugn - 1] == 0)
        Ugn--;
    if (Ugn == 0)
        goto cleanup;

    /* w = (|u| +- 3 u^2/4) / 2 (+ when T is too small), in units
       B^-(n+2): u^2 = Ug^2 B^-(n+2), from the top limbs of Ug (the s
       dropped ones, s = m - 4, move it by under 2 B^-2 units) */
    s = FLINT_MAX(0, FLINT_MIN(m - 4, Ugn - 1));
    hn = Ugn - s;
    SQ = TMP_ALLOC((2 * hn + 1) * sizeof(ulong));
    flint_mpn_sqr(SQ, Ug + s, hn);
    sqn = 2 * hn;
    {
        slong i = n + 2 - 2 * s;
        if (sqn > i)
        {
            nn_ptr sq = SQ + i;
            slong l = sqn - i;
            ulong cy;
            cy = mpn_mul_1(sq, sq, l, 3);
            sq[l] = cy;
            l++;
            mpn_rshift(sq, sq, l, 2);
            while (l > 0 && sq[l - 1] == 0)
                l--;
            if (l > 0)
            {
                Ug[Ugn] = 0;
                if (Unegative)
                    mpn_add(Ug, Ug, Ugn + 1, sq, l);
                else
                    mpn_sub(Ug, Ug, Ugn + 1, sq, l);
                Ugn++;
            }
        }
    }
    mpn_rshift(Ug, Ug, Ugn, 1);
    /* drop the guard limbs: w at weight B^-n */
    Ug += 2;
    Ugn -= 2;
    while (Ugn > 0 && Ug[Ugn - 1] == 0)
        Ugn--;

    if (Ugn > 0)
    {
        /* Q +-= T w */
        V = TMP_ALLOC((Tn + Ugn - (m - 2)) * sizeof(ulong));
        flint_mpn_mulmid(V, T, Tn, Ug, Ugn, m - 2, Tn + Ugn);
        Vhigh = V + 2;
        Vhighn = Tn + Ugn - (m - 2) - 2;

        if (Unegative)
            mpn_add(Q, Q, n + 2, Vhigh, Vhighn);
        else
            mpn_sub(Q, Q, n + 2, Vhigh, Vhighn);
    }

cleanup:
    TMP_END;
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
   T ~= 1/sqrt(A) carried at half precision (fixed_rsqrt_newton) */
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
