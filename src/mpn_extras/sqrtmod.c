/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmpcompat.h"
#include "mpn_extras.h"

/* A length-n mpn read as an mpz; GMP wants the top limb to be nonzero. */
static void
_mpn_roinit(mpz_t z, nn_srcptr x, mp_size_t n)
{
    while (n > 0 && x[n - 1] == 0)
        n--;

    mpz_roinit_n(z, (mp_srcptr) x, n);
}

/* The reverse, zero-padding to exactly n limbs; z must be in [0, 2^(n B)). */
static void
_mpn_set_mpz(nn_ptr r, mp_size_t n, mpz_srcptr z)
{
    mp_size_t zn = z->_mp_size;

    flint_mpn_copyi(r, z->_mp_d, zn);
    flint_mpn_zero(r + zn, n - zn);
}

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

/*
    Writing d - 1 = q 2^s with q odd, s = 1 and s = 2 are closed forms and
    otherwise we walk down the 2-part of the group (Tonelli and Shanks)
    using the least quadratic nonresidue.

    Returns 1 on success, 0 when the Jacobi symbol rules out a root, and
    -1 when the algorithm itself failed, which means d is not prime.

    The exponentiations go to mpz_powm, which is hard to beat, but the
    descent is worth doing in mpn arithmetic with a precomputed inverse.
    It runs on values scaled by 2^norm, the representation in which
    flint_mpn_mulmod_preinvn works; scaling is exact in both directions
    because a scaled value is by construction below the scaled modulus.
*/
static int
_flint_mpn_sqrtmod(nn_ptr res, nn_srcptr a, nn_srcptr d, mp_size_t n,
        nn_srcptr dinv, flint_bitcnt_t norm)
{
    mpz_t az, dz, q, e, b, c, r, t, u;
    slong s, i, j, m, iter;
    ulong k;
    int success = 1;

    _mpn_roinit(az, a, n);
    _mpn_roinit(dz, d, n);

    /* a Jacobi symbol of -1 rules out a root without knowing d to be prime */
    if (mpz_jacobi(az, dz) == -1)
    {
        flint_mpn_zero(res, n);
        return 0;
    }

    if (flint_mpz_cmp_ui(az, 1) <= 0)
    {
        flint_mpn_copyi(res, a, n);
        return 1;
    }

    mpz_init(q);
    mpz_init(e);
    mpz_init(b);
    mpz_init(c);
    mpz_init(r);
    mpz_init(t);
    mpz_init(u);

    /* d - 1 = q 2^s, q odd */
    flint_mpz_sub_ui(q, dz, 1);
    s = mpz_scan1(q, 0);
    mpz_tdiv_q_2exp(q, q, s);

    if (s == 1)
    {
        /* d = 3 mod 4, so the root is a^((d+1)/4) */
        flint_mpz_add_ui(e, dz, 1);
        mpz_tdiv_q_2exp(e, e, 2);
        mpz_powm(r, az, e, dz);
    }
    else if (s == 2)
    {
        /*
            d = 5 mod 8. Atkin: with b = (2a)^((d-5)/8) and y = 2 a b^2,
            a root is a b (y - 1). This is one exponentiation, where the
            variant in the old fmpz_sqrtmod needed a second one half of
            the time.
        */
        flint_mpz_sub_ui(e, dz, 5);
        mpz_tdiv_q_2exp(e, e, 3);

        mpz_mul_2exp(t, az, 1); mpz_mod(t, t, dz);          /* t = 2a */
        mpz_powm(b, t, e, dz);
        mpz_mul(u, b, b); mpz_mod(u, u, dz);
        mpz_mul(u, u, t); mpz_mod(u, u, dz);                /* u = y */
        flint_mpz_sub_ui(u, u, 1);
        mpz_mul(b, b, az); mpz_mod(b, b, dz);
        mpz_mul(r, b, u); mpz_mod(r, r, dz);
    }
    else
    {
        nn_ptr dnormed, dinv_tmp, R, T, C, B, U, ONE;
        TMP_INIT;

        /*
            a^((q-1)/2) yields both a^((q+1)/2) and a^q for a multiplication
            each, so the descent is set up with two exponentiations, not
            three.
        */
        flint_mpz_sub_ui(e, q, 1);
        mpz_tdiv_q_2exp(e, e, 1);

        mpz_powm(u, az, e, dz);                             /* u = a^((q-1)/2) */
        mpz_mul(r, u, az); mpz_mod(r, r, dz);               /* r = a^((q+1)/2) */
        mpz_mul(t, r, u); mpz_mod(t, t, dz);                /* t = a^q */

        /*
            The least quadratic nonresidue, which is well within reach.
            Only odd candidates need testing: 2 is a residue for
            d = 1 mod 8, and then so is any even number below the first
            odd nonresidue.
        */
        for (k = 3; ; k += 2)
        {
            flint_mpz_set_ui(b, k);

            if (mpz_jacobi(b, dz) == -1)
                break;

            if (k >= (UWORD(1) << 20))
            {
                /* only reachable if d is not actually prime */
                success = -1;
                goto cleanup;
            }
        }

        mpz_powm(c, b, q, dz);                              /* c = k^q */

        TMP_START;

        dnormed = TMP_ALLOC((8 * n) * sizeof(ulong));
        R = dnormed + n; T = R + n; C = T + n;
        B = C + n; U = B + n; ONE = U + n; dinv_tmp = ONE + n;

        if (norm)
            mpn_lshift(dnormed, d, n, norm);
        else
            flint_mpn_copyi(dnormed, d, n);

        if (dinv == NULL)
        {
            flint_mpn_preinvn(dinv_tmp, dnormed, n);
            dinv = dinv_tmp;
        }

        _mpn_set_mpz(R, n, r);
        _mpn_set_mpz(T, n, t);
        _mpn_set_mpz(C, n, c);

        if (norm)
        {
            mpn_lshift(R, R, n, norm);
            mpn_lshift(T, T, n, norm);
            mpn_lshift(C, C, n, norm);
        }

        flint_mpn_zero(ONE, n);
        ONE[0] = UWORD(1) << norm;

        m = s;

        for (iter = 0; mpn_cmp(T, ONE, n) != 0; iter++)
        {
            /* at most s - 1 descents if d is prime */
            if (iter >= s)
            {
                success = -1;
                break;
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
                break;
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

        if (success == 1)
        {
            if (norm)
                mpn_rshift(R, R, n, norm);

            flint_mpn_copyi(res, R, n);
        }

        TMP_END;
        goto cleanup;
    }

    _mpn_set_mpz(res, n, r);

cleanup:
    if (success != 1)
        flint_mpn_zero(res, n);

    mpz_clear(q);
    mpz_clear(e);
    mpz_clear(b);
    mpz_clear(c);
    mpz_clear(r);
    mpz_clear(t);
    mpz_clear(u);

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
