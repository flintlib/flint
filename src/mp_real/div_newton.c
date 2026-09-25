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
#include "mp_real.h"
#include "impl.h"

/* Approximate fixed-point inversion and division, ported from
   radix_inv_approx / radix_div_approx onto mpn arithmetic (the
   radix-B generality drops out: with B = 2^64 only the
   middle-product branches survive, division by two is a shift, and
   the binary basecases need no radix conversion).  The _newton suffix
   flags the contracts: not ulp-accurate, two-sided errors of a few
   ulps as documented in mp_real.h.

   The Newton and Karp-Markstein iterations are carried verbatim;
   see radix/div.c for the error analysis (with the LIMB_RADIX(radix)
   large / mulhigh case applying throughout). */

void
_mp_real_inv_newton_basecase(nn_ptr Q, nn_srcptr A, slong An, slong n)
{
    slong Vn = FLINT_MIN(An, n + 1);

    /* floor(B^(Vn + n) / V) for the top Vn limbs V of A, with n + 2
       limbs, by division: unlike flint_mpn_inv, _flint_mpn_inv_basecase
       never calls _mp_real_inv_newton at the same or a higher precision */
    _flint_mpn_inv_basecase(Q, A + An - Vn, Vn, Vn + n);
}

/* Tuned on this machine (see dev/notes); must be at least 4 for the
   indexing/padding operations in the code to be valid. */
#ifndef MP_REAL_INV_NEWTON_CUTOFF
#define MP_REAL_INV_NEWTON_CUTOFF 14
#endif

/* Precisions taking a third-order step (below); the others take the
   second-order one.  Both have the same contract, so a recursion can
   mix them freely. */
#ifndef MP_REAL_INV_NEWTON_3_LOW_MIN
#define MP_REAL_INV_NEWTON_3_LOW_MIN 24
#endif
#ifndef MP_REAL_INV_NEWTON_3_LOW_MAX
#define MP_REAL_INV_NEWTON_3_LOW_MAX 80
#endif
#ifndef MP_REAL_INV_NEWTON_3_MIN
#define MP_REAL_INV_NEWTON_3_MIN 1000
#endif

/* Second-order step: T ~= 1/A at m ~ n/2 limbs, Q = T + T (1 - A T). */
static void
_mp_real_inv_newton_step2(nn_ptr Q, nn_srcptr A, slong An, slong n)
{
    nn_ptr U, V, T, Vhigh, Uhigh;
    nn_srcptr Ahigh;
    slong m, Uhighn, Tn, Ahighn, Vhighn;
    int Unegative;
    TMP_INIT;

    /* T ~= 1/A with m fraction limbs, 2 integral limbs, stored in
       the high part of Q */
    m = (n + 1) / 2 + 1;
    _mp_real_inv_newton(Q + n - m, A, An, m);

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

/* Third-order step: T ~= 1/A at m ~ n/3 limbs, e = A T - 1,

       1/A = T / (1 + e) = T (1 - e + e^2 - ...),

   with |e| <= 4 B^-m and the tail |T e^3| <= 64 B^(1 - 3m) below the
   rounding for 3m >= n + 3.  e is formed as in the second-order step
   (the band of A T with the unit wrapped into a control limb, merely
   wider: weights B^-(n+2) to B^-(m-2)), e^2 from its top ~n - 2m
   limbs, and w = |e| -+ e^2 in units two limbs below the output before
   the guard limbs are dropped; Q = T -+ T w by a middle product.  The
   step costs the band (A times an m-limb T, 2n/3 limbs kept), a square
   of about n/3 limbs and the product T w (n/3 by 2n/3 limbs): 1.1 to
   1.5 times a second-order step (whose residual is a single middle
   product), against recursion covering a third of the precision
   instead of half, so it pays off only at some precisions. */
static void
_mp_real_inv_newton_step3(nn_ptr Q, nn_srcptr A, slong An, slong n)
{
    nn_ptr U, V, T, Ug, SQ;
    nn_srcptr Ahigh;
    slong m, Tn, Ahighn, Ugn, control_limb, low_trunc, sqn, s, hn;
    slong A_low_zeroes = 0;
    int Unegative;
    TMP_INIT;

    /* 3m >= n + 3: |e| <= 4 B^-m, tail |T e^3| <= 64 B^(1-3m) */
    m = (n + 3 + 2) / 3;
    _mp_real_inv_newton(Q + n - m, A, An, m);

    TMP_START;
    T = Q + n - m;
    Tn = m + 1 + (T[m + 1] != 0);
    flint_mpn_zero(Q, n - m);

    /* band of A' T at weights B^-(n+2) .. B^-(m-2) plus control limb */
    control_limb = n - m + 2;
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
    low_trunc = m - 1;

    U = TMP_ALLOC((2 + (control_limb + 1) + 1) * sizeof(ulong));
    {
        slong lo = low_trunc - A_low_zeroes;
        slong hi = lo + 2 + (control_limb + 1);
        if (lo >= 0)
            flint_mpn_mulmid(U, Ahigh, Ahighn, T, Tn, lo, hi);
        else
        {
            flint_mpn_zero(U, FLINT_MIN(-lo, hi - lo));
            if (hi > 0)
                flint_mpn_mulmid(U - lo, Ahigh, Ahighn, T, Tn, 0, hi);
        }
    }

    /* Unegative: A T < 1, e < 0 */
    Unegative = (U[2 + control_limb] != 0);
    Ug = U;
    Ugn = control_limb + 2;
    if (Unegative)
        mpn_neg(Ug, Ug, Ugn);
    while (Ugn > 0 && Ug[Ugn - 1] == 0)
        Ugn--;
    if (Ugn == 0)
        goto cleanup;

    /* w = |e| + e^2 (e < 0) or |e| - e^2, units B^-(n+2) */
    s = FLINT_MAX(0, FLINT_MIN(m - 4, Ugn - 1));
    hn = Ugn - s;
    SQ = TMP_ALLOC(2 * hn * sizeof(ulong));
    flint_mpn_sqr(SQ, Ug + s, hn);
    sqn = 2 * hn;
    {
        slong i = n + 2 - 2 * s;
        if (sqn > i)
        {
            nn_ptr sq = SQ + i;
            slong l = sqn - i;
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
    Ug += 2;
    Ugn -= 2;
    while (Ugn > 0 && Ug[Ugn - 1] == 0)
        Ugn--;

    if (Ugn > 0)
    {
        V = TMP_ALLOC((Tn + Ugn - (m - 2)) * sizeof(ulong));
        flint_mpn_mulmid(V, T, Tn, Ug, Ugn, m - 2, Tn + Ugn);
        if (Unegative)
            mpn_add(Q, Q, n + 2, V + 2, Tn + Ugn - (m - 2) - 2);
        else
            mpn_sub(Q, Q, n + 2, V + 2, Tn + Ugn - (m - 2) - 2);
    }

cleanup:
    TMP_END;
}

void
_mp_real_inv_newton(nn_ptr Q, nn_srcptr A, slong An, slong n)
{
    if (n <= MP_REAL_INV_NEWTON_CUTOFF)
        _mp_real_inv_newton_basecase(Q, A, An, n);
    else if ((n >= MP_REAL_INV_NEWTON_3_LOW_MIN && n <= MP_REAL_INV_NEWTON_3_LOW_MAX)
            || n >= MP_REAL_INV_NEWTON_3_MIN)
        _mp_real_inv_newton_step3(Q, A, An, n);
    else
        _mp_real_inv_newton_step2(Q, A, An, n);
}

void
_mp_real_div_newton_invmul(nn_ptr Q, nn_srcptr B, slong Bn, nn_srcptr A,
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

    _mp_real_inv_newton(Q, A, An, n);
    /* should be a high multiplication */
    if (n + 2 >= Bn)
        flint_mpn_mul(T, Q, n + 2, B, Bn);
    else
        flint_mpn_mul(T, B, Bn, Q, n + 2);
    flint_mpn_copyi(Q, T + Bn, n + 2);

    TMP_END;
}

#ifndef MP_REAL_DIV_NEWTON_CUTOFF
#define MP_REAL_DIV_NEWTON_CUTOFF 12
#endif

/* Karp-Markstein division: evaluates TB + T * (B - A * TB) with
   T ~= 1/A carried at half precision */
void
_mp_real_div_newton(nn_ptr Q, nn_srcptr B, slong Bn, nn_srcptr A,
    slong An, slong n)
{
    nn_ptr U, V, T, Vhigh, Uhigh;
    nn_srcptr Ahigh;
    slong m, Uhighn, Tn, Ahighn, Vhighn;
    int Unegative;
    TMP_INIT;

    if (n <= MP_REAL_DIV_NEWTON_CUTOFF)
    {
        _mp_real_div_newton_invmul(Q, B, Bn, A, An, n);
        return;
    }

    m = (n + 1) / 2 + 1;

    _mp_real_inv_newton(Q + n - m, A, An, m);

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
