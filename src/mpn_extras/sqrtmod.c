/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

/* A length-n mpn read as an mpz; GMP wants the top limb to be nonzero. */
static void
_mpn_roinit(mpz_t z, nn_srcptr x, mp_size_t n)
{
    while (n > 0 && x[n - 1] == 0)
        n--;

    mpz_roinit_n(z, (mp_srcptr) x, n);
}

/*
    GMP exposes the Jacobi symbol only on mpz, so this is the one place
    where a read-only view is taken; nothing is copied or allocated.
*/
int
flint_mpn_is_square_mod(nn_srcptr a, nn_srcptr d, mp_size_t n)
{
    mpz_t az, dz;

    if (n == 1)
        return n_jacobi_unsigned(a[0], d[0]) != -1;

    _mpn_roinit(az, a, n);
    _mpn_roinit(dz, d, n);

    return mpz_jacobi(az, dz) != -1;
}

/* To and from the representation flint_mpn_mulmod_preinvn works in. */
static void
_scale(nn_ptr r, nn_srcptr x, mp_size_t n, flint_bitcnt_t norm)
{
    if (norm)
        mpn_lshift(r, x, n, norm);
    else
        flint_mpn_copyi(r, x, n);
}

static void
_unscale(nn_ptr r, nn_srcptr x, mp_size_t n, flint_bitcnt_t norm)
{
    if (norm)
        mpn_rshift(r, x, n, norm);
    else
        flint_mpn_copyi(r, x, n);
}

/*
    Writing d - 1 = q 2^s with q odd, s = 1 and s = 2 are closed forms and
    otherwise we walk down the 2-part of the group (Tonelli and Shanks)
    using the least quadratic nonresidue.

    Returns 1 on success, 0 when the Jacobi symbol rules out a root, and
    -1 when the algorithm itself failed, which means d is not prime.
*/
static int
_flint_mpn_sqrtmod(nn_ptr res, nn_srcptr a, nn_srcptr d, mp_size_t n,
        nn_srcptr dinv, flint_bitcnt_t norm)
{
    nn_ptr dnormed, dinv_tmp, A, ONE, R, T, C, B, U, q, e;
    slong i, j, s, m, iter;
    mp_size_t qn;
    ulong k;
    int success = 1;
    TMP_INIT;

    /* a Jacobi symbol of -1 rules out a root without knowing d to be prime */
    if (!flint_mpn_is_square_mod(a, d, n))
    {
        flint_mpn_zero(res, n);
        return 0;
    }

    if (a[0] <= 1 && flint_mpn_zero_p(a + 1, n - 1))
    {
        flint_mpn_copyi(res, a, n);
        return 1;
    }

    TMP_START;

    dnormed = TMP_ALLOC((11 * n + 2) * sizeof(ulong));
    dinv_tmp = dnormed + n;
    A = dinv_tmp + n; ONE = A + n; R = ONE + n; T = R + n;
    C = T + n; B = C + n; U = B + n; q = U + n; e = q + n;

    _scale(dnormed, d, n, norm);

    if (dinv == NULL)
    {
        flint_mpn_preinvn(dinv_tmp, dnormed, n);
        dinv = dinv_tmp;
    }

    flint_mpn_zero(ONE, n);
    ONE[0] = UWORD(1) << norm;

    /* d - 1 = q 2^s, q odd */
    mpn_sub_1(q, d, n, 1);

    for (i = 0; q[i] == 0; i++)
        ;

    s = i * FLINT_BITS + flint_ctz(q[i]);
    qn = n - i;

    if (s % FLINT_BITS)
        mpn_rshift(q, q + i, qn, s % FLINT_BITS);
    else
        flint_mpn_copyi(q, q + i, qn);

    if (s == 1)
    {
        /* d = 3 mod 4, so the root is a^((d+1)/4) */
        e[n] = mpn_add_1(e, d, n, 1);
        mpn_rshift(e, e, n + 1, 2);

        _scale(A, a, n, norm);
        flint_mpn_powmod_preinvn(R, A, e, n + 1, n, dnormed, dinv, norm);
    }
    else if (s == 2)
    {
        /*
            d = 5 mod 8. Atkin: with b = (2a)^((d-5)/8) and y = 2 a b^2,
            a root is a b (y - 1). This is one exponentiation, where the
            variant in the old fmpz_sqrtmod needed a second one half of
            the time.
        */
        mpn_sub_1(e, d, n, 5);
        mpn_rshift(e, e, n, 3);

        _scale(A, a, n, norm);
        flint_mpn_addmod_n(T, A, A, dnormed, n);                    /* T = 2a */
        flint_mpn_powmod_preinvn(B, T, e, n, n, dnormed, dinv, norm);
        flint_mpn_mulmod_preinvn(U, B, B, n, dnormed, dinv, norm);
        flint_mpn_mulmod_preinvn(U, U, T, n, dnormed, dinv, norm);   /* U = y */
        flint_mpn_submod_n(U, U, ONE, dnormed, n);
        flint_mpn_mulmod_preinvn(B, B, A, n, dnormed, dinv, norm);
        flint_mpn_mulmod_preinvn(R, B, U, n, dnormed, dinv, norm);
    }
    else
    {
        /*
            a^((q-1)/2) yields both a^((q+1)/2) and a^q for a multiplication
            each, so the descent is set up with two exponentiations, not
            three.
        */
        mpn_sub_1(e, q, qn, 1);
        mpn_rshift(e, e, qn, 1);

        _scale(A, a, n, norm);
        flint_mpn_powmod_preinvn(U, A, e, qn, n, dnormed, dinv, norm);
        flint_mpn_mulmod_preinvn(R, U, A, n, dnormed, dinv, norm);   /* a^((q+1)/2) */
        flint_mpn_mulmod_preinvn(T, R, U, n, dnormed, dinv, norm);   /* a^q */

        /*
            The least quadratic nonresidue, which is well within reach.
            Only odd candidates need testing: 2 is a residue for
            d = 1 mod 8, and then so is any even number below the first
            odd nonresidue.
        */
        for (k = 3; ; k += 2)
        {
            flint_mpn_zero(B, n);
            B[0] = k;

            if (!flint_mpn_is_square_mod(B, d, n))
                break;

            if (k >= (UWORD(1) << 20))
            {
                /* only reachable if d is not actually prime */
                success = -1;
                goto cleanup;
            }
        }

        _scale(B, B, n, norm);
        flint_mpn_powmod_preinvn(C, B, q, qn, n, dnormed, dinv, norm);

        m = s;

        for (iter = 0; mpn_cmp(T, ONE, n) != 0; iter++)
        {
            /* at most s - 1 descents if d is prime */
            if (iter >= s)
            {
                success = -1;
                goto cleanup;
            }

            /* the least i with 0 < i < m and T^(2^i) = 1 */
            flint_mpn_copyi(U, T, n);

            for (i = 1; i < m; i++)
            {
                flint_mpn_mulmod_preinvn(U, U, U, n, dnormed, dinv, norm);

                if (mpn_cmp(U, ONE, n) == 0)
                    break;
            }

            if (i >= m)
            {
                /* again, only reachable for a modulus that is not prime */
                success = -1;
                goto cleanup;
            }

            /* B = C^(2^(m - i - 1)) */
            flint_mpn_copyi(B, C, n);

            for (j = 0; j < m - i - 1; j++)
                flint_mpn_mulmod_preinvn(B, B, B, n, dnormed, dinv, norm);

            m = i;
            flint_mpn_mulmod_preinvn(C, B, B, n, dnormed, dinv, norm);
            flint_mpn_mulmod_preinvn(T, T, C, n, dnormed, dinv, norm);
            flint_mpn_mulmod_preinvn(R, R, B, n, dnormed, dinv, norm);
        }
    }

    _unscale(res, R, n, norm);

cleanup:
    if (success != 1)
        flint_mpn_zero(res, n);

    TMP_END;

    return success;
}

int
flint_mpn_sqrtmod_preinv(nn_ptr res, nn_srcptr a, nn_srcptr d, mp_size_t n,
        nn_srcptr dinv, flint_bitcnt_t norm)
{
    if (n == 1)
    {
        ulong r;

        if (n_jacobi_unsigned(a[0], d[0]) == -1)
        {
            res[0] = 0;
            return 0;
        }

        r = n_sqrtmod(a[0], d[0]);
        res[0] = r;

        /* the symbol did not rule a root out, so only a lie about d can */
        return (r != 0 || a[0] == 0) ? 1 : -1;
    }

    return _flint_mpn_sqrtmod(res, a, d, n, dinv, norm);
}

int
flint_mpn_sqrtmod(nn_ptr res, nn_srcptr a, nn_srcptr d, mp_size_t n)
{
    /* only the descent wants the inverse, and it computes its own */
    return flint_mpn_sqrtmod_preinv(res, a, d, n, NULL, flint_clz(d[n - 1]));
}
