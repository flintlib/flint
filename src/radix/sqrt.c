/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
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
        radix_mul_1(U, T, m + 2, a, radix);   /* mod B^(m+2): carry discarded */
        cy = (U[m + 1] == 0);
        /* a*y^2-1 or 1-a*y^2 */
        if (!cy)
            radix_neg(U, U, m + 2, radix);

        radix_divrem_two(U, U, m + 2, radix);
        Un = m + 2;
        MPN_NORM(U, Un);
        if (Un + n - 2 * m <= 0)
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


/* Integer square root with remainder via binary (machine) radix
   conversion; mirrors radix_divrem_via_mpn. */
void
radix_sqrtrem_via_mpn(nn_ptr s, nn_ptr r, nn_srcptr a, slong an, const radix_t radix)
{
    nn_ptr ba, br, bs, tmp;
    slong ban, brn, bsn, sn;
    TMP_INIT;

    FLINT_ASSERT(an >= 1);
    FLINT_ASSERT(a[an - 1] != 0);

    sn = (an + 1) / 2;

    TMP_START;
    tmp = TMP_ALLOC((2 * an + (an + 3) / 2) * sizeof(ulong));
    ba = tmp;
    br = ba + an;
    bs = br + an;

    ban = radix_get_mpn(ba, a, an, radix);
    bsn = (ban + 1) / 2;

    brn = mpn_sqrtrem(bs, (r == NULL) ? NULL : br, ba, ban);

    /* Need to do radix conversion in temporary space as radix conversion may
       need an extra output scratch limb. */
    slong qn, rn, qn2, rn2;
    nn_ptr tmp2, s2, r2;

    qn2 = radix_set_mpn_need_alloc(bsn, radix);
    tmp2 = TMP_ALLOC(qn2 * sizeof(ulong));
    s2 = tmp2;
    qn = radix_set_mpn(s2, bs, bsn, radix);
    FLINT_ASSERT(qn <= sn);
    flint_mpn_copyi(s, s2, qn);
    flint_mpn_zero(s + qn, sn - qn);

    if (r != NULL)
    {
        if (brn == 0)
        {
            flint_mpn_zero(r, sn + 1);
        }
        else
        {
            rn2 = radix_set_mpn_need_alloc(brn, radix);
            tmp2 = TMP_ALLOC(rn2 * sizeof(ulong));
            r2 = tmp2;
            rn = radix_set_mpn(r2, br, brn, radix);
            FLINT_ASSERT(rn <= sn + 1);
            flint_mpn_copyi(r, r2, rn);
            flint_mpn_zero(r + rn, (sn + 1) - rn);
        }
    }

    TMP_END;
}

/*
Input: A is a fixed-point number with An fraction limbs representing
a value alpha in [B^-2, 1), i.e. at least one of the top two limbs of A
is nonzero. (The wider normalization, compared to [B^-1, 1) for
radix_inv_approx, is what radix_sqrtrem needs to handle inputs with an
odd number of limbs, since the square-root scaling B^(1/2) is not a limb
shift.)

Output: n + 2 limbs approximating 1/sqrt(alpha) in (1, B] with n fraction
limbs and 2 integral limbs. The second integral limb is nonzero only when
alpha is close to B^-2.

Q = floor(B^(n+Vn+g) / floor(sqrt(V * B^(Vn+2g)))) where V is the top
Vn = min(An, n+2) limbs of A and g = max(0, n + 2 - Vn) supplies working
precision when A is short: the relative error of floor(sqrt(W)) is up to
1/sqrt(W), so W must have ~2(n+2) limbs regardless of An. Two guard limbs
(rather than inv_approx's one) are needed because the relative effect of
truncating A is e/(2*alpha) <= B^2*e/2 under the wider normalization.
*/
void
radix_rsqrt_approx_basecase(nn_ptr Q, nn_srcptr A, slong An, slong n, const radix_t radix)
{
    nn_ptr W, U, bW, bU, bs, bq, br, q2;
    nn_srcptr V;
    slong Vn, g, Wn, Un, bWn, bUn, bsn, bqn, qn, qn2;
    slong guard_limbs = 2;
    TMP_INIT;

    Vn = FLINT_MIN(An, n + guard_limbs);
    g = FLINT_MAX(0, n + guard_limbs - Vn);
    V = A + An - Vn;

    Wn = 2 * (Vn + g);
    Un = n + Vn + g + 1;

    TMP_START;
    /* Radix scratch for W = V * B^(Vn + 2g) and U = B^(n + Vn + g).
       The square root and the reciprocal are both computed in the binary
       (machine) radix, avoiding radix conversion roundtrips for the
       intermediate square root. */
    W = TMP_ALLOC((Wn + Un) * sizeof(ulong));
    U = W + Wn;

    flint_mpn_zero(W, Vn + 2 * g);
    flint_mpn_copyi(W + Vn + 2 * g, V, Vn);

    /* The top limb of A (hence of W) may be zero when alpha < B^-1. */
    while (Wn > 0 && W[Wn - 1] == 0)
        Wn--;
    FLINT_ASSERT(Wn >= 2 * (Vn + g) - 1);

    flint_mpn_zero(U, Un - 1);
    U[Un - 1] = 1;

    /* bs = floor(sqrt(W)) */
    bW = TMP_ALLOC((Wn + Un + (Wn + 3) / 2 + Un + (Wn + 3) / 2) * sizeof(ulong));
    bU = bW + Wn;
    bs = bU + Un;
    bq = bs + (Wn + 3) / 2;
    br = bq + Un;

    bWn = radix_get_mpn(bW, W, Wn, radix);
    bsn = (bWn + 1) / 2;
    mpn_sqrtrem(bs, NULL, bW, bWn);

    /* bq = floor(U / bs); the quotient is < B^(n+2) since
       bs >= B^(Vn+g-1) */
    bUn = radix_get_mpn(bU, U, Un, radix);
    bqn = bUn - bsn + 1;
    FLINT_ASSERT(bUn >= bsn);
    mpn_tdiv_qr(bq, br, 0, bU, bUn, bs, bsn);

    /* Need to do radix conversion in temporary space as radix conversion may
       need an extra output scratch limb. */
    qn2 = radix_set_mpn_need_alloc(bqn, radix);
    q2 = TMP_ALLOC(qn2 * sizeof(ulong));
    qn = radix_set_mpn(q2, bq, bqn, radix);
    FLINT_ASSERT(qn <= n + 2);
    flint_mpn_copyi(Q, q2, qn);
    flint_mpn_zero(Q + qn, (n + 2) - qn);

    TMP_END;
}

/* Must be large enough that m - 2 >= 2 and 2m - n - 2 >= 0 hold in the
   recursive cases of radix_rsqrt_approx and radix_sqrt_approx.
   Empirically 8 is close to optimal for word-size limb radices; for
   small-digit radices (where the basecase's binary operands are much
   smaller than the radix limb count) a larger cutoff would be slightly
   better at small n. */
#define RADIX_RSQRT_NEWTON_CUTOFF 8

/*
Claim: the absolute error is bounded by 4*B^(-n)/sqrt(A). This is not tight.

Proof sketch: let m = (n+1)/2 + 1 + [B <= 3], so that B^(-2m) <= B^(-n-2)
(<= B^(-n-4) for B = 3).

The initial approximation is assumed to satisfy

    Y = (1/sqrt(A)) * (1 + e1)  where  |e1| <= C*B^(-m)  with C = 4.

We compute Y' = Y + Y * (1 - (A + e2) * Y^2 + e3) / 2 + e4 where

    e2 = B^(-n-2) is the error from truncating A (two guard limbs are
    needed since the fixed point of the iteration moves by
    e2/(2*A) <= B^2*e2/2 relative),

    e3 = B^(-n) + B^(-n)/2 + {(2m+4)*B^(-n-1)}
    e4 = B^(-n) + {(m+2)*B^(-n-1)}

are the errors from the truncated products, the halving of the residual,
and the final multiplication by Y; the terms in curly brackets appear only
when mulhigh is used, which requires B > 2*(m+2). All contributions to the
relative error of Y' are

    (3/2)*e1^2 + (1/2)*|e1|^3 + (1 + e1)*(Y^2*e2 + e3)/2 + e4*sqrt(A)

and since Y^2*e2 <= B^(-n)/2 * (1+e1)^2 and B^(-2m) <= B^(-n-2), verifying
that this is below 4*B^(-n) for all B >= 3 is straightforward numerics.
(The basecase satisfies the bound with C < 2.)
*/
void
radix_rsqrt_approx(nn_ptr Q, nn_srcptr A, slong An, slong n, const radix_t radix)
{
    if (n <= RADIX_RSQRT_NEWTON_CUTOFF)
    {
        radix_rsqrt_approx_basecase(Q, A, An, n, radix);
    }
    else
    {
        nn_ptr U, V, W, T, Vhigh, Uhigh;
        nn_srcptr Ahigh;
        slong m, Uhighn, Tn, Wn, Ahighn, Vhighn;
        int Unegative;
        TMP_INIT;

        /*
        Compute T ~= 1 / sqrt(A) as a fixed-point number with m fraction
        limbs, 2 integral limbs and store in the high part of Q.
        */
        m = (n + 1) / 2 + 1 + (LIMB_RADIX(radix) <= 3);
        radix_rsqrt_approx(Q + n - m, A, An, m, radix);

        TMP_START;

        T = Q + n - m;
        Tn = m + 1 + (T[m + 1] != 0);

        flint_mpn_zero(Q, n - m);

        /* W = T^2 with 2m fraction limbs (full product; all but O(1) of its
           low limbs are needed below, so truncation is not worthwhile). */
        Wn = 2 * Tn;
        W = TMP_ALLOC(Wn * sizeof(ulong));
        radix_mul(W, T, Tn, T, Tn, radix);

        /*
        Let A' be the high n+2 limbs of A. We compute U ~= 1 - A'*T^2 as a
        fixed-point number with n fraction limbs, using the same middle
        product + control limb construction as radix_inv_approx: the low
        limbs of the product can be truncated because we only need n
        fraction limbs, and the high limbs are known a priori to be all
        0 or B-1 since A'*T^2 approximates 1.
        */
        slong control_limb = n / 2 + 2;

        /* If An < n + 2, zero-extend A. */
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

        // The full product has n+2+2m fraction limbs; remove 2m+2 fraction
        // limbs to get n fraction limbs.
        slong low_trunc_without_mulhigh = 2 * m + 2;
        // When using mulhigh, include 2 guard limbs to get an accurate high
        // product.
        slong low_trunc_with_mulhigh = 2 * m;

        /* A_low_zeroes <= n <= 2m - 2, so unlike in radix_inv_approx no
           third case for very short A is needed. */
        FLINT_ASSERT(A_low_zeroes <= low_trunc_with_mulhigh - 2);

        if (LIMB_RADIX(radix) > 2 * (m + 2))
        {
            U = TMP_ALLOC((2 + (control_limb + 1)) * sizeof(ulong));
            radix_mulmid(U, Ahigh, Ahighn, W, Wn,
                low_trunc_with_mulhigh - A_low_zeroes,
                low_trunc_with_mulhigh - A_low_zeroes + 2 + (control_limb + 1),
                radix);
            Uhigh = U + 2;
        }
        else
        {
            // Emulate high product by computing the full product and
            // throwing away excess fraction limbs
            U = TMP_ALLOC((low_trunc_without_mulhigh - A_low_zeroes + (control_limb) + 1) * sizeof(ulong));
            radix_mulmid(U, Ahigh, Ahighn, W, Wn, 0,
                low_trunc_without_mulhigh - A_low_zeroes + (control_limb) + 1,
                radix);
            Uhigh = U + low_trunc_without_mulhigh - A_low_zeroes;
        }

        Unegative = (Uhigh[control_limb] != 0);
        if (Unegative)
            radix_neg(Uhigh, Uhigh, control_limb, radix);

        /* There may be additional cancellation. */
        Uhighn = control_limb;
        while (Uhighn > 0 && Uhigh[Uhighn - 1] == 0)
            Uhighn--;

        if (Uhighn != 0)
        {
            /* The correction is T * (1 - A*T^2) / 2. */
            radix_divrem_two(Uhigh, Uhigh, Uhighn, radix);
            while (Uhighn > 0 && Uhigh[Uhighn - 1] == 0)
                Uhighn--;
        }

        if (Uhighn != 0)
        {
            /* Compute V = T * U with n fraction limbs. */
            if (LIMB_RADIX(radix) > 2 * (m + 2))
            {
                // V = T * |1 - A' * T^2| / 2 with n+2 fraction limbs
                V = TMP_ALLOC((Tn + Uhighn - (m - 2)) * sizeof(ulong));
                radix_mulmid(V, T, Tn, Uhigh, Uhighn, m - 2, Tn + Uhighn, radix);
                Vhigh = V + 2;
                Vhighn = Tn + Uhighn - (m - 2) - 2;
            }
            else
            {
                // V = T * |1 - A' * T^2| / 2 with n+m fraction limbs
                V = TMP_ALLOC((Tn + Uhighn) * sizeof(ulong));
                radix_mul(V, T, Tn, Uhigh, Uhighn, radix);
                Vhigh = V + m;
                Vhighn = Tn + Uhighn - m;
            }

            /* Finally, compute T + T * (1 - A'T^2)/2 or T - T * (A'T^2 - 1)/2
            depending on the sign of U. */
            if (Unegative)
                radix_add(Q, Q, n + 2, Vhigh, Vhighn, radix);
            else
                radix_sub(Q, Q, n + 2, Vhigh, Vhighn, radix);
        }

        TMP_END;
    }
}

void
radix_sqrt_approx_rsqrtmul(nn_ptr Q, nn_srcptr A, slong An, slong n, const radix_t radix)
{
    nn_ptr T;
    nn_srcptr A2;
    slong A2n;
    TMP_INIT;
    TMP_START;

    A2n = FLINT_MIN(An, n + 2);
    A2 = A + An - A2n;

    T = TMP_ALLOC((A2n + n + 2) * sizeof(ulong));

    radix_rsqrt_approx(Q, A, An, n, radix);
    /* Should be a high multiplication. */
    radix_mul(T, Q, n + 2, A2, A2n, radix);
    flint_mpn_copyi(Q, T + A2n, n + 2);

    TMP_END;
}

/* Karp-Markstein square root: analogous to radix_div_approx
   but evaluates S + T * (A - S^2) / 2 with S = T * A instead of
   TB + T * (B - A * TB).

   Input as for radix_rsqrt_approx: A has An fraction limbs with value
   alpha in [B^-2, 1). Output: n + 2 limbs approximating sqrt(alpha) in
   [B^-1, 1) with n fraction limbs; the integral limbs are normally zero
   but the value can round to 1.000... when alpha is close to 1.

   Claim: the absolute error is bounded by 4*B^(-n)/sqrt(A). Note that,
   as for radix_rsqrt_approx, the error is proportional to 1/sqrt(A)
   rather than to the output sqrt(A): the corrections of the last
   iteration are multiplied by T ~= 1/sqrt(A) <= B, which amplifies
   their B^(-n)-size rounding errors accordingly. */
void
radix_sqrt_approx(nn_ptr Q, nn_srcptr A, slong An, slong n, const radix_t radix)
{
    nn_ptr U, V, T, TS, Vhigh, Uhigh, TShigh;
    slong m, Un, Uhighn, Tn, Vhighn;
    int Unegative;
    TMP_INIT;

    if (n <= 4)
    {
        radix_sqrt_approx_rsqrtmul(Q, A, An, n, radix);
        return;
    }

    m = (n + 1) / 2 + 1 + (LIMB_RADIX(radix) <= 3);

    radix_rsqrt_approx(Q + n - m, A, An, m, radix);

    TMP_START;
    T = Q + n - m;
    Tn = m + 1 + (T[m + 1] != 0);

    flint_mpn_zero(Q, n - m);

    // Compute S = T * A2 ~= sqrt(A) as a fixed-point number with m fraction
    // limbs. Unlike in radix_div_approx, two extra limbs of the truncated
    // operand are needed: the truncation error is multiplied by T <= B.
    slong A2n = FLINT_MIN(m + 2, An);
    nn_srcptr A2 = A + An - A2n;

    if (LIMB_RADIX(radix) > 2 * (m + 2) && A2n > 2)
    {
        TS = TMP_ALLOC((Tn + 2) * sizeof(ulong));
        radix_mulmid(TS, T, Tn, A2, A2n, A2n - 2, A2n + Tn, radix);
        TShigh = TS + 2;
    }
    else
    {
        TS = TMP_ALLOC((A2n + Tn) * sizeof(ulong));
        radix_mul(TS, T, Tn, A2, A2n, radix);
        TShigh = TS + A2n;
    }

    slong control_limb = n / 2 + 2;

    Un = control_limb + 1;

    /* U = band of S^2 at weights B^(-n) .. B^(control_limb - n), plus
       2 guard limbs in the mulhigh case. The full square has 2m fraction
       limbs, so the band starts 2m - n limbs up. */
    if (LIMB_RADIX(radix) > 2 * (m + 2) && 2 * m - n - 2 >= 0)
    {
        U = TMP_ALLOC((2 + Un) * sizeof(ulong));
        radix_mulmid(U, TShigh, Tn, TShigh, Tn, 2 * m - n - 2, 2 * m - n + Un, radix);
        Uhigh = U + 2;
    }
    else
    {
        U = TMP_ALLOC((2 * m - n + Un) * sizeof(ulong));
        radix_mulmid(U, TShigh, Tn, TShigh, Tn, 0, 2 * m - n + Un, radix);
        Uhigh = U + 2 * m - n;
    }

    // Compute (A - S^2) with n fraction limbs, with sign
    if (An > n - Un)
    {
        if (An >= n)
        {
            nn_srcptr Ahigh = A + An - n;
            radix_sub(Uhigh, Ahigh, Un, Uhigh, Un, radix);
        }
        else
        {
            slong nz = n - An;
            ulong cy = radix_neg(Uhigh, Uhigh, nz, radix);
            radix_sub(Uhigh + nz, A, Un - nz, Uhigh + nz, Un - nz, radix);
            radix_sub(Uhigh + nz, Uhigh + nz, Un - nz, &cy, 1, radix);
        }

        Unegative = 1;
        if (radix_cmp_bn_half(Uhigh, Un, radix) > 0)
        {
            radix_neg(Uhigh, Uhigh, Un, radix);
            Unegative = 0;
        }
    }
    else
    {
        Unegative = 0;
        if (radix_cmp_bn_half(Uhigh, Un, radix) > 0)
        {
            radix_neg(Uhigh, Uhigh, Un, radix);
            Unegative = 1;
        }
    }

    Uhighn = Un;
    while (Uhighn > 0 && Uhigh[Uhighn - 1] == 0)
        Uhighn--;

    if (Uhighn != 0)
    {
        /* The correction is T * (A - S^2) / 2. */
        radix_divrem_two(Uhigh, Uhigh, Uhighn, radix);
        while (Uhighn > 0 && Uhigh[Uhighn - 1] == 0)
            Uhighn--;
    }

    if (Uhighn != 0)
    {
        if (LIMB_RADIX(radix) > 2 * (m + 2))
        {
            V = TMP_ALLOC((Tn + Uhighn - (m - 2)) * sizeof(ulong));
            radix_mulmid(V, T, Tn, Uhigh, Uhighn, m - 2, Tn + Uhighn, radix);
            Vhigh = V + 2;
            Vhighn = Tn + Uhighn - (m - 2) - 2;
        }
        else
        {
            V = TMP_ALLOC((Tn + Uhighn) * sizeof(ulong));
            radix_mul(V, T, Tn, Uhigh, Uhighn, radix);
            Vhigh = V + m;
            Vhighn = Tn + Uhighn - m;
        }

        // Overwrite T with S
        Q[n + 1] = 0;
        flint_mpn_copyi(Q + n - m, TShigh, Tn);

        if (Unegative)
            radix_add(Q, Q, n + 2, Vhigh, Vhighn, radix);
        else
            radix_sub(Q, Q, n + 2, Vhigh, Vhighn, radix);
    }
    else
    {
        // Overwrite T with S
        Q[n + 1] = 0;
        flint_mpn_copyi(Q + n - m, TShigh, Tn);
    }

    TMP_END;
}

/* Instead of computing sqrt(A) with ~1 ulp error and doing a full
   multiplication to check the remainder, compute a few fraction limbs.
   If the first fraction limb is not too close to the limb boundary, we
   have the correct root and the remainder can be determined by a low
   multiplication (or omitted if we only want the root).

   Three extra fraction limbs are used (one more than division's
   PREINV2_EXTRA_LIMBS): the sqrt_approx error is C*B^(-n2)/sqrt(alpha)
   with alpha as small as B^-2 when an is odd, i.e. up to C*B^(-n2+1),
   so n2 = sn + 3 leaves error <= C*B^(-2) at the integer scale.

   The fast acceptance is gated on LIMB_RADIX > 8: with the error bound
   C*B^(-2), C = 4, the check q[-1] in [2, B-2] certifies the integer
   part only when C < B. For smaller radices we always run the (always
   correct) verification loops. */
void
radix_sqrtrem_newton_karp_markstein(nn_ptr s, nn_ptr r, nn_srcptr a, slong an, const radix_t radix)
{
    nn_ptr S, P, t, q, Apad;
    nn_srcptr Aview;
    slong sn, n2, viewn;
    ulong one = 1;
    TMP_INIT;

    FLINT_ASSERT(an >= 1);
    FLINT_ASSERT(a[an - 1] != 0);

    sn = (an + 1) / 2;

    TMP_START;

    /* View a as a fixed-point number with 2*sn fraction limbs (zero-padded
       on top when an is odd), representing alpha = a * B^(-2sn) in
       [B^-2, 1). Then sqrt(a) = sqrt(alpha) * B^sn. */
    viewn = 2 * sn;
    if (an == viewn)
    {
        Aview = a;
    }
    else
    {
        Apad = TMP_ALLOC(viewn * sizeof(ulong));
        flint_mpn_copyi(Apad, a, an);
        Apad[viewn - 1] = 0;
        Aview = Apad;
    }

    n2 = sn + 3;

    S = TMP_ALLOC((n2 + 2) * sizeof(ulong));
    radix_sqrt_approx(S, Aview, viewn, n2, radix);

    FLINT_ASSERT(S[n2 + 1] == 0);
    /* q[0], ..., q[sn-1] are the candidate digits of floor(sqrt(a));
       q[-1] is the first fraction digit; q[sn] the integral overflow. */
    q = S + n2 - sn;

    if (LIMB_RADIX(radix) > 8 && q[sn] == 0
        && q[-1] > 1 && q[-1] < LIMB_RADIX(radix) - 1)
    {
        /* Certified: q = floor(sqrt(a)). */
        if (r != NULL)
        {
            P = TMP_ALLOC(FLINT_MAX(an, sn + 1) * sizeof(ulong));
            radix_mulmid(P, q, sn, q, sn, 0, an, radix);
            radix_sub(P, a, an, P, an, radix);
            FLINT_ASSERT(an <= sn + 1 || flint_mpn_zero_p(P + sn + 1, an - (sn + 1)));
            if (an < sn + 1)
                flint_mpn_zero(P + an, (sn + 1) - an);
            flint_mpn_copyi(r, P, sn + 1);
        }
        flint_mpn_copyi(s, q, sn);
        TMP_END;
        return;
    }

    if (q[sn] != 0)
        radix_sub(q, q, sn + 1, &one, 1, radix);

    /* Verification: adjust q by O(1) steps. P holds q^2, later a - q^2,
       over an an+1 limb window; t holds 2q+1 over sn+1 limbs. */
    P = TMP_ALLOC((an + 1) * sizeof(ulong));
    t = TMP_ALLOC((sn + 1) * sizeof(ulong));

    flint_mpn_zero(P, an + 1);
    radix_mulmid(P, q, sn, q, sn, 0, FLINT_MIN(2 * sn, an + 1), radix);

    while (P[an] != 0 || mpn_cmp(P, a, an) > 0)
    {
        radix_sub(q, q, sn, &one, 1, radix);
        /* t = 2q + 1 for the new q; q_old^2 - (2q+1) = q^2 */
        t[sn] = radix_add(t, q, sn, q, sn, radix);
        radix_add(t, t, sn + 1, &one, 1, radix);
        radix_sub(P, P, an + 1, t, sn + 1, radix);
    }

    /* P = a - q^2 (nonnegative; P[an] is zero at this point) */
    radix_sub(P, a, an, P, an, radix);

    for (;;)
    {
        /* t = 2q + 1 */
        t[sn] = radix_add(t, q, sn, q, sn, radix);
        radix_add(t, t, sn + 1, &one, 1, radix);

        /* stop when a - q^2 < 2q + 1, i.e. (q+1)^2 > a */
        if (flint_mpn_zero_p(P + sn + 1, (an + 1) - (sn + 1))
            && mpn_cmp(P, t, sn + 1) < 0)
            break;

        radix_add(q, q, sn, &one, 1, radix);
        radix_sub(P, P, an + 1, t, sn + 1, radix);
    }

    if (r != NULL)
        flint_mpn_copyi(r, P, sn + 1);
    flint_mpn_copyi(s, q, sn);

    TMP_END;
}

/* Below this size, converting to binary and using mpn_sqrtrem wins;
   above it, Newton-Karp-Markstein is faster on average. The crossover
   grows as the limbs get narrower, since the binary operands shrink
   proportionally: with b bits per limb, an an-limb input is only about
   an*b/64 binary limbs. The formula below approximates measured
   crossovers of ~40 limbs at b = 63, ~100 at b = 30, ~450 at b = 10.
   (As for radix_divrem, the crossover is not monotone: with B = 2^63
   the mpn path is competitive again in a window around an ~ 256..1024
   before the radix multiplication switches to fft_small, so a
   dual-range dispatch in the style of radix_divrem could be used.) */
RADIX_INLINE slong
_radix_sqrtrem_km_cutoff(const radix_t radix)
{
    slong b = FLINT_BIT_COUNT(LIMB_RADIX(radix));
    return 2500 / b + 20000 / (b * b);
}

/* s = floor(sqrt(a)), sn = (an+1)/2 limbs. If r is not NULL, it is set
   to a - s^2 (at most 2s), zero-padded to sn + 1 limbs. */
void
radix_sqrtrem(nn_ptr s, nn_ptr r, nn_srcptr a, slong an, const radix_t radix)
{
    FLINT_ASSERT(an >= 1);
    FLINT_ASSERT(a[an - 1] != 0);

    if (an == 1)
    {
        ulong s0 = n_sqrt(a[0]);
        s[0] = s0;
        if (r != NULL)
        {
            r[0] = a[0] - s0 * s0;
            r[1] = 0;
        }
        return;
    }

    if (an < _radix_sqrtrem_km_cutoff(radix))
        radix_sqrtrem_via_mpn(s, r, a, an, radix);
    else
        radix_sqrtrem_newton_karp_markstein(s, r, a, an, radix);
}

/* s = floor(sqrt(a)), sn = (an+1)/2 limbs. Returns 1 iff a is a
   perfect square. */
int
radix_sqrt(nn_ptr s, nn_srcptr a, slong an, const radix_t radix)
{
    slong sn;
    int exact;
    nn_ptr r;
    TMP_INIT;

    FLINT_ASSERT(an >= 1);
    FLINT_ASSERT(a[an - 1] != 0);

    if (an == 1)
    {
        ulong s0 = n_sqrt(a[0]);
        s[0] = s0;
        return s0 * s0 == a[0];
    }

    sn = (an + 1) / 2;

    TMP_START;
    r = TMP_ALLOC((sn + 1) * sizeof(ulong));
    radix_sqrtrem(s, r, a, an, radix);
    exact = flint_mpn_zero_p(r, sn + 1);
    TMP_END;

    return exact;
}
