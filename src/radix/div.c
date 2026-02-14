/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "longlong.h"
#include "radix.h"

/* todo: optimise */
ulong
radix_divrem_1(nn_ptr res, nn_srcptr a, slong an, ulong d, const radix_t radix)
{
    slong i;
    ulong q, r, hi, lo, B = LIMB_RADIX(radix);

    q = a[an - 1] / d;
    r = a[an - 1] % d;
    res[an - 1] = q;

    for (i = an - 2; i >= 0; i--)
    {
        umul_ppmm(hi, lo, r, B);
        add_ssaaaa(hi, lo, hi, lo, 0, a[i]);
        udiv_qrnnd(q, r, hi, lo, d);
        res[i] = q;
    }

    return r;
}

/* todo: optimise */
void
radix_divexact_1(nn_ptr res, nn_srcptr a, slong an, ulong d, const radix_t radix)
{
    ulong r;
    r = radix_divrem_1(res, a, an, d, radix);
    FLINT_ASSERT(r == 0);
    (void) r;
}

void
radix_divrem_via_mpn(nn_ptr q, nn_ptr r, nn_srcptr a, slong an, nn_srcptr b, slong bn, const radix_t radix)
{
    /* Binary (machine) radix */
    nn_ptr bq, br, ba, bb, tmp;
    slong bqn, brn, ban, bbn, alloc;
    TMP_INIT;

    FLINT_ASSERT(an >= bn);
    FLINT_ASSERT(bn >= 1);
    FLINT_ASSERT(b[bn - 1] != 0);

    /* Ensure that ban >= bbn */
    if (flint_mpn_zero_p(a + bn, an - bn) && mpn_cmp(a, b, bn) < 0)
    {
        if (r != NULL)
            flint_mpn_copyi(r, a, bn);
        if (q != NULL)
            flint_mpn_zero(q, an - bn + 1);
        return;
    }

    /* todo: tighten allocations when radix is much smaller than the machine radix */
    alloc = ((an - bn + 2) + (bn + 1) + (an + 1) + (bn + 1));
    TMP_START;
    tmp = TMP_ALLOC(alloc * sizeof(ulong));
    bq = tmp;
    br = bq + (an - bn + 2);
    ba = br + (bn + 1);
    bb = ba + (an + 1);

    ban = radix_get_mpn(ba, a, an, radix);
    bbn = radix_get_mpn(bb, b, bn, radix);
    bqn = ban - bbn + 1;
    brn = bbn;

    FLINT_ASSERT(ban >= bbn);

    mpn_tdiv_qr(bq, br, 0, ba, ban, bb, bbn);

    /* Need to do radix conversion in temporary space as radix conversion may
       need an extra output scratch limb. */

    slong qn, rn, qn2, rn2;
    nn_ptr tmp2, q2, r2;

    if (q != NULL)
    {
        qn2 = radix_set_mpn_need_alloc(bqn, radix);
        tmp2 = TMP_ALLOC(qn2 * sizeof(ulong));
        q2 = tmp2;

        qn = radix_set_mpn(q2, bq, bqn, radix);

        if (r != NULL)
        {
            /* r = a - bq */
            /* Use br as scratch space to allow aliasing a with r */
            radix_mulmid(br, b, bn, q2, qn, 0, bn, radix);
            radix_sub(r, a, bn, br, bn, radix);
        }

        flint_mpn_copyi(q, q2, qn);
        flint_mpn_zero(q + qn, (an - bn + 1) - qn);
    }
    else
    {
        rn2 = radix_set_mpn_need_alloc(brn, radix);
        tmp2 = TMP_ALLOC(rn2 * sizeof(ulong));
        r2 = tmp2;

        rn = radix_set_mpn(r2, br, brn, radix);
        flint_mpn_copyi(r, r2, rn);
        flint_mpn_zero(r + rn, bn - rn);
    }

    TMP_END;
}

void
radix_inv_approx_basecase(nn_ptr Q, nn_srcptr A, slong An, slong n, const radix_t radix)
{
    nn_ptr U;
    nn_srcptr V;
    slong Un, Vn;
    slong guard_limbs;
    TMP_INIT;

    guard_limbs = 1;
    Vn = FLINT_MIN(An, n + guard_limbs);
    Un = Vn + n + 1;

    TMP_START;
    U = TMP_ALLOC(Un * sizeof(ulong));
    V = A + An - Vn;

    flint_mpn_zero(U, Un - 1);
    U[Un - 1] = 1;

    if (Vn == 1)
        radix_divrem_1(Q, U, Un, V[0], radix);
    else
        radix_divrem_via_mpn(Q, NULL, U, Un, V, Vn, radix);

    TMP_END;
}

/* Must be at least 4 */
#define RADIX_INV_NEWTON_CUTOFF 8

void
radix_inv_approx(nn_ptr Q, nn_srcptr A, slong An, slong n, const radix_t radix)
{
    if (n <= RADIX_INV_NEWTON_CUTOFF)
    {
        radix_inv_approx_basecase(Q, A, An, n, radix);
    }
    else
    {
        nn_ptr U, V, T, Vhigh, Uhigh;
        nn_srcptr Ahigh;
        slong m, Uhighn, Tn, Ahighn, Vhighn;
        int Unegative;
        TMP_INIT;

        /*
        Compute T ~= 1 / A as a fixed-point number with m fraction limbs,
        2 integral limbs and store in the high part of Q.
        */
        m = n / 2 + 1 + (LIMB_RADIX(radix) == 2);
        radix_inv_approx(Q + n - m, A, An, m, radix);

        TMP_START;

        /* Tn = number of limbs in T, including m fraction limbs and
           1 or 2 integral limbs. */
        T = Q + n - m;
        Tn = m + 1 + (T[m + 1] != 0);

        /* Set the low limbs of Q to zero. TODO: optimize out. */
        flint_mpn_zero(Q, n - m);

        /*
        Let A' be the high n+1 limbs of A. We compute U ~= 1 - T*A' as a
        fixed-point number with n fraction limbs.

        The full product T*A' has n+m+1 fraction limbs, 2 integral limbs.
        We can replace this ~1.5n-limb full product with a middle product
        of ~n/2 limbs, truncating both ~n/2 high and low limbs.

        The reason why we can truncate the low ~n/2 limbs is that we only
        need ~n fraction limbs in the output.

        The reason why we can truncate the high ~n/2 limbs is that we know
        a priori that T*A' approximates 1 to within this accuracy. The higher
        fraction limbs are thus all 0 or B-1 and need not be computed
        explicitly. Computing just one such limb allows determining the
        sign of U.

        TODO: the algorithm could be simplified by choosing the initial
        approximation T so that TA' <= 1 always holds, so that we can
        do away with the sign check.
        */

        slong control_limb = n / 2 + 2;

        /* If An < n + 1, zero-extend A. */
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

        // The full product has n+m+1 fraction limbs; remove m+1 fraction limbs
        // to get n fraction limbs.
        slong low_trunc_without_mulhigh = m + 1;
        // When using mulhigh, include 2 guard limbs to get an accurate high
        // product.
        slong low_trunc_with_mulhigh = m - 1;

        if (LIMB_RADIX(radix) > m + 2 && A_low_zeroes <= low_trunc_with_mulhigh)
        {
            U = TMP_ALLOC((2 + (control_limb + 1)) * sizeof(ulong));
            radix_mulmid(U, Ahigh, Ahighn, T, Tn,
                low_trunc_with_mulhigh - A_low_zeroes,
                low_trunc_with_mulhigh - A_low_zeroes + 2 + (control_limb + 1),
                radix);
            Uhigh = U + 2;
        }
        else if (A_low_zeroes <= low_trunc_without_mulhigh)
        {
            // Emulate high product by computing the full product and
            // throwing away excess fraction limbs
            U = TMP_ALLOC((low_trunc_without_mulhigh - A_low_zeroes + (control_limb) + 1) * sizeof(ulong));
            radix_mulmid(U, Ahigh, Ahighn, T, Tn, 0,
                low_trunc_without_mulhigh - A_low_zeroes + (control_limb) + 1,
                radix);
            Uhigh = U + low_trunc_without_mulhigh - A_low_zeroes;
        }
        else
        {
            // A is short; we need to zero-extend to get n fraction limbs
            U = TMP_ALLOC((control_limb + 1) * sizeof(ulong));
            radix_mulmid(U + A_low_zeroes - low_trunc_without_mulhigh,
                Ahigh, Ahighn, T, Tn, 0,
                control_limb + 1 - (A_low_zeroes - low_trunc_without_mulhigh),
                radix);
            flint_mpn_zero(U, A_low_zeroes - low_trunc_without_mulhigh);
            Uhigh = U;
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
            /* Compute V = T * U with n fraction limbs. */
            if (LIMB_RADIX(radix) > m + 2)
            {
                // V = T * |1 - A' * T| with n+2 fraction limbs
                V = TMP_ALLOC((Tn + Uhighn - (m - 2)) * sizeof(ulong));
                radix_mulmid(V, T, Tn, Uhigh, Uhighn, m - 2, Tn + Uhighn, radix);
                Vhigh = V + 2;
                Vhighn = Tn + Uhighn - (m - 2) - 2;
            }
            else
            {
                // V = T * |1 - A' * T| with n+m fraction limbs
                V = TMP_ALLOC((Tn + Uhighn) * sizeof(ulong));
                radix_mul(V, T, Tn, Uhigh, Uhighn, radix);
                Vhigh = V + m;
                Vhighn = Tn + Uhighn - m;
            }

            /* Finally, compute T + T * (1 - A'T) or T - T * (A'T - 1) depending
            on the sign of U. */
            if (Unegative)
                radix_add(Q, Q, n + 2, Vhigh, Vhighn, radix);
            else
                radix_sub(Q, Q, n + 2, Vhigh, Vhighn, radix);
        }

        TMP_END;
    }
}

void
radix_divrem_preinv(nn_ptr Q, nn_ptr R, nn_srcptr A, slong An, nn_srcptr B, slong Bn, nn_srcptr Binv, slong Binvn, const radix_t radix)
{
    nn_ptr T, U, q, r;
    slong n;
    ulong one = 1;
    TMP_INIT;
    TMP_START;

    n = An - Bn + 1;

    /* Todo: use the unbalanced division algorithm from radix_divrem */
    if (Binvn < n)
        flint_throw(FLINT_ERROR, "radix_divrem_preinv: inverse has too few limbs");

    if (LIMB_RADIX(radix) > n + 2)
    {
        U = TMP_ALLOC((n + 3) * sizeof(ulong));
        radix_mulmid(U, A + Bn - 1, n, Binv + Binvn - n, n + 2, An - Bn, 2 * n + 2, radix);
        q = U + 2;
    }
    else
    {
        U = TMP_ALLOC((2 * n + 2) * sizeof(ulong));
        radix_mul(U, A + Bn - 1, n, Binv + Binvn - n, n + 2, radix);
        q = U + An - Bn + 2;
    }

    T = TMP_ALLOC((Bn + n) * sizeof(ulong));
    r = T;

    if (q[n] != 0)
        radix_sub(q, q, n + 1, &one, 1, radix);

    radix_mul(r, q, n, B, Bn, radix);

    while (r[An] != 0 || mpn_cmp(r, A, An) > 0)
    {
        radix_sub(r, r, An + 1, B, Bn, radix);
        radix_sub(q, q, n, &one, 1, radix);
    }

    radix_sub(r, A, An, r, An, radix);

    while (r[Bn] != 0 || mpn_cmp(r, B, Bn) >= 0)
    {
        radix_add(q, q, n, &one, 1, radix);
        radix_sub(r, r, Bn + 1, B, Bn, radix);
    }

    flint_mpn_copyi(Q, q, An - Bn + 1);
    flint_mpn_copyi(R, r, Bn);

    TMP_END;
}

void
radix_divrem_newton(nn_ptr Q, nn_ptr R, nn_srcptr A, slong An, nn_srcptr B, slong Bn, const radix_t radix)
{
    nn_ptr T;
    slong n;
    TMP_INIT;
    TMP_START;

    if (flint_mpn_zero_p(A + Bn, An - Bn) && mpn_cmp(A, B, Bn) < 0)
    {
        flint_mpn_copyi(R, A, Bn);
        flint_mpn_zero(Q, An - Bn + 1);
        return;
    }

    n = An - Bn + 1;
    T = TMP_ALLOC((n + 2) * sizeof(ulong));
    radix_inv_approx(T, B, Bn, n, radix);
    radix_divrem_preinv(Q, R, A, An, B, Bn, T, n, radix);
    TMP_END;
    return;
}

void
radix_divrem(nn_ptr q, nn_ptr r, nn_srcptr a, slong an, nn_srcptr b, slong bn, const radix_t radix)
{
    FLINT_ASSERT(an >= bn);
    FLINT_ASSERT(bn >= 1);
    FLINT_ASSERT(b[bn - 1] != 0);

    if (bn == 1)
    {
        r[0] = radix_divrem_1(q, a, an, b[0], radix);
        return;
    }

    if (bn >= 48)
    {
        radix_divrem_newton(q, r, a, an, b, bn, radix);
        return;
    }

    /* Unbalanced division: view the an x bn division as an N x 1
       division with length-bn limbs. */
    if (an > 2 * bn)
    {
        nn_ptr T, R;
        slong i;
        slong antop;
        ulong cy;
        TMP_INIT;

        FLINT_ASSERT(an > 2 * bn);

        TMP_START;
        R = TMP_ALLOC(3 * bn * sizeof(ulong));
        T = R + bn;

        i = (an + bn - 1) / bn - 2;
        antop = an - i * bn;
        FLINT_ASSERT(antop >= bn);

        radix_divrem(q + i * bn, R, a + i * bn, antop, b, bn, radix);
        cy = q[i * bn];

        for (i--; i >= 0; i--)
        {
            flint_mpn_copyi(T, a + i * bn, bn);
            flint_mpn_copyi(T + bn, R, bn);
            radix_divrem(q + i * bn, R, T, 2 * bn, b, bn, radix);

            /* The recursive call writes 2 bn + 1 limbs, but we know that the partial
               quotient fits in bn limbs. Restore the limb that was overwritten with
               a zero. */
            q[(i + 1) * bn] = cy;
            cy = q[i * bn];
        }

        flint_mpn_copyi(r, R, bn);
        TMP_END;
        return;
    }

    radix_divrem_via_mpn(q, r, a, an, b, bn, radix);
}

