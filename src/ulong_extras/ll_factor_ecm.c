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
   The elliptic curve method specialised to two-word moduli.

   Curve selection (Suyama parametrisation) is done once per curve with
   mpn arithmetic, which is negligible against the cost of the stages and
   keeps the modular inversion simple.  Stage 1 and stage 2 run entirely
   in two-word Montgomery arithmetic; this is where essentially all the
   time goes.

   The algorithm mirrors fmpz_factor_ecm: Montgomery curve XZ-only
   arithmetic, a binary ladder for stage 1, and the standard continuation
   with a primorial step P for stage 2.
*/

#include "longlong.h"
#include "mpn_extras.h"
#include "ulong_extras.h"

#include "ll_mont2.h"

#if FLINT_BITS == 64

typedef struct
{
    ulong n0, n1, np;
    ulong one[2];    /* R mod n */
    ulong R2[2];     /* R^2 mod n */
}
ll_mont_s;

typedef ll_mont_s ll_mont_t[1];

#define MMUL(r, a, b) \
    LL_MULMOD((r)[1], (r)[0], (a)[1], (a)[0], (b)[1], (b)[0], M->n1, M->n0, M->np)
#define MADD(r, a, b) \
    LL_ADDMOD((r)[1], (r)[0], (a)[1], (a)[0], (b)[1], (b)[0], M->n1, M->n0)
#define MSUB(r, a, b) \
    LL_SUBMOD((r)[1], (r)[0], (a)[1], (a)[0], (b)[1], (b)[0], M->n1, M->n0)
#define MSQR(r, a) \
    LL_SQRMOD((r)[1], (r)[0], (a)[1], (a)[0], M->n1, M->n0, M->np)

static void
ll_mont_init(ll_mont_t M, ulong n1, ulong n0)
{
    ulong num[5], quot[4], nvec[2];

    M->n0 = n0;
    M->n1 = n1;
    M->np = LL_NEG_NINV(n0);

    nvec[0] = n0;
    nvec[1] = n1;

    /* R^2 mod n, where R = 2^(2*FLINT_BITS) */
    num[0] = num[1] = num[2] = num[3] = 0;
    num[4] = 1;
    mpn_tdiv_qr(quot, M->R2, 0, num, 5, nvec, 2);

    /* R mod n = R^2 / R */
    LL_MULMOD(M->one[1], M->one[0], M->R2[1], M->R2[0], 0, 1, n1, n0, M->np);
}

/* convert the two-word value {a1,a0} in [0,n) into Montgomery form */
static void
ll_to_mont(nn_ptr r, ulong a1, ulong a0, const ll_mont_t M)
{
    LL_MULMOD(r[1], r[0], a1, a0, M->R2[1], M->R2[0], M->n1, M->n0, M->np);
}

/*
   Montgomery curve arithmetic on (X : Z) coordinates.  In both routines
   every input is consumed before any output is written, so the
   destination may alias the inputs other than the difference point.
*/

static void
ll_xz_double(nn_ptr X2, nn_ptr Z2, nn_srcptr X, nn_srcptr Z,
             nn_srcptr a24, const ll_mont_t M)
{
    ulong u[2], v[2], w[2];

    MADD(u, X, Z);  MSQR(u, u);      /* u = (X+Z)^2 */
    MSUB(v, X, Z);  MSQR(v, v);      /* v = (X-Z)^2 */
    MMUL(X2, u, v);                     /* X2 = u*v */
    MSUB(w, u, v);                      /* w = 4*X*Z */
    MMUL(u, w, a24);
    MADD(u, u, v);
    MMUL(Z2, w, u);
}

/* (X3:Z3) = (X1:Z1) + (X2:Z2) where (X0:Z0) is their difference */
static void
ll_xz_add(nn_ptr X3, nn_ptr Z3, nn_srcptr X1, nn_srcptr Z1,
          nn_srcptr X2, nn_srcptr Z2, nn_srcptr X0, nn_srcptr Z0,
          const ll_mont_t M)
{
    ulong u[2], v[2], w[2];

    MSUB(u, X2, Z2);
    MADD(v, X1, Z1);
    MMUL(u, u, v);
    MSUB(v, X1, Z1);
    MADD(w, X2, Z2);
    MMUL(v, v, w);
    MADD(w, u, v);
    MSUB(v, v, u);
    MSQR(w, w);
    MSQR(v, v);
    MMUL(X3, Z0, w);
    MMUL(Z3, X0, v);
}

/* (RX:RZ) = k*(X:Z) by the Montgomery ladder */
static void
ll_xz_mul(nn_ptr RX, nn_ptr RZ, nn_srcptr X, nn_srcptr Z, ulong k,
          nn_srcptr a24, const ll_mont_t M)
{
    ulong x1[2], z1[2], x2[2], z2[2];
    slong i;

    if (k == 0)
    {
        RX[0] = RX[1] = 0;
        RZ[0] = RZ[1] = 0;
        return;
    }

    x1[0] = X[0]; x1[1] = X[1];
    z1[0] = Z[0]; z1[1] = Z[1];

    if (k == 1)
    {
        RX[0] = x1[0]; RX[1] = x1[1];
        RZ[0] = z1[0]; RZ[1] = z1[1];
        return;
    }

    ll_xz_double(x2, z2, X, Z, a24, M);

    /* invariant: (x2:z2) - (x1:z1) = (X:Z) */
    for (i = (slong) FLINT_BIT_COUNT(k) - 2; i >= 0; i--)
    {
        if ((k >> i) & 1)
        {
            ll_xz_add(x1, z1, x1, z1, x2, z2, X, Z, M);
            ll_xz_double(x2, z2, x2, z2, a24, M);
        }
        else
        {
            ll_xz_add(x2, z2, x1, z1, x2, z2, X, Z, M);
            ll_xz_double(x1, z1, x1, z1, a24, M);
        }
    }

    RX[0] = x1[0]; RX[1] = x1[1];
    RZ[0] = z1[0]; RZ[1] = z1[1];
}

/* gcd(v, n); returns 1 and writes a proper divisor to factor */
static int
ll_gcd_factor(nn_ptr factor, nn_srcptr v, ulong n1, ulong n0)
{
    ulong nvec[2], vvec[2], g[2];
    mp_size_t gs;

    if ((v[0] | v[1]) == 0)
        return 0;

    nvec[0] = n0; nvec[1] = n1;
    vvec[0] = v[0]; vvec[1] = v[1];

    gs = flint_mpn_gcd_full(g, nvec, 2, vvec, v[1] ? 2 : 1);

    if (gs == 1)
    {
        if (g[0] == 1)
            return 0;
        factor[0] = g[0];
        factor[1] = 0;
        return 1;
    }

    if (g[0] == n0 && g[1] == n1)
        return 0;

    factor[0] = g[0];
    factor[1] = g[1];
    return 1;
}

/* stage 2 step; P is a primorial, tab[(i-mmin)*(maxj+1) + j] flags a prime hit */
#define LL_ECM_P 210

static int
ll_ecm_curve(nn_ptr factor, ulong n1, ulong n0, ulong sigma, ulong B1,
             const ll_mont_t M, const unsigned char * tab, ulong mmin, ulong mmax)
{
    ulong X[2], Z[2], a24[2];
    ulong sig[2], u[2], v[2], t[2], num[2], den[2], c5[2], c3[2], c16[2];
    int ret = 0;

    /*
       Suyama's parametrisation, in Montgomery form throughout:
           u = sigma^2 - 5,  v = 4*sigma,
           x0 = u^3,  z0 = v^3,
           a24 = (A+2)/4 = (v-u)^3 * (3u+v) / (16 * u^3 * v).
       Only the final division needs an inversion.
    */
    ll_to_mont(sig, 0, sigma, M);
    ll_to_mont(c5, 0, 5, M);
    ll_to_mont(c3, 0, 3, M);
    ll_to_mont(c16, 0, 16, M);

    MSQR(u, sig);
    MSUB(u, u, c5);                       /* u = sigma^2 - 5 */

    MADD(v, sig, sig);
    MADD(v, v, v);                        /* v = 4*sigma */

    MSQR(X, u);  MMUL(X, X, u);           /* x0 = u^3 */
    MSQR(Z, v);  MMUL(Z, Z, v);           /* z0 = v^3 */

    MSUB(t, v, u);
    MSQR(num, t);  MMUL(num, num, t);     /* (v-u)^3 */
    MMUL(t, u, c3);
    MADD(t, t, v);                        /* 3u + v */
    MMUL(num, num, t);

    MMUL(den, c16, X);
    MMUL(den, den, v);                    /* 16 * u^3 * v */

    /*
       Invert den.  mpn_gcdext wants both operands given as exactly
       n limbs with the modulus normalised, as in mpn_mod_inv; the gcd it
       returns is the factor of n when the inverse does not exist.
    */
    {
        ulong dord[2], nvec[2], tmp[2], one_ord[2], g[4], sv[4];
        mp_size_t gsize, ssize;

        /* leave Montgomery form: den * 1 / R gives the ordinary residue */
        one_ord[0] = 1; one_ord[1] = 0;
        MMUL(dord, den, one_ord);

        if ((dord[0] | dord[1]) == 0)
            goto cleanup;                 /* degenerate curve, skip it */

        nvec[0] = M->n0; nvec[1] = M->n1;
        tmp[0] = dord[0]; tmp[1] = dord[1];

        gsize = mpn_gcdext(g, sv, &ssize, tmp, 2, nvec, 2);

        if (!(gsize == 1 && g[0] == 1))
        {
            /* den shares a factor with n */
            if (!(gsize == 2 && g[0] == M->n0 && g[1] == M->n1))
            {
                factor[0] = g[0];
                factor[1] = (gsize > 1) ? g[1] : 0;
                ret = 1;
            }
            goto cleanup;
        }

        flint_mpn_zero(sv + FLINT_ABS(ssize), 2 - FLINT_ABS(ssize));
        if (ssize < 0)
            flint_mpn_negmod_n(sv, sv, nvec, 2);

        /* back into Montgomery form and combine */
        ll_to_mont(t, sv[1], sv[0], M);
        MMUL(a24, num, t);
    }

    /* ---------------- stage 1 ---------------- */
    {
        n_primes_t iter;
        ulong p;

        n_primes_init(iter);
        while ((p = n_primes_next(iter)) <= B1)
        {
            ulong q = p;
            while (q <= B1 / p)
                q *= p;
            ll_xz_mul(X, Z, X, Z, q, a24, M);
        }
        n_primes_clear(iter);
    }

    if (ll_gcd_factor(factor, Z, n1, n0))
    {
        ret = 1;
        goto cleanup;
    }

    if ((Z[0] | Z[1]) == 0)     /* point at infinity, nothing more to learn */
        goto cleanup;

    /* ---------------- stage 2 ---------------- */
    {
        const ulong P = LL_ECM_P;
        ulong maxj = (P + 1) / 2;
        ulong nj = (maxj >> 1) + 1;
        ulong row = maxj + 1;
        nn_ptr arrx, arrz;
        ulong Q0x2[2], Q0z2[2];
        ulong Qx[2], Qz[2], Rx[2], Rz[2], Qdx[2], Qdz[2];
        ulong g[4][2], a[2], b[2], sx[2], sz[2];
        ulong i, j;
        int gi = 0;

        arrx = flint_malloc(2 * nj * sizeof(ulong));
        arrz = flint_malloc(2 * nj * sizeof(ulong));

        /* arr[j>>1] = j*Q0 for odd j */
        arrx[0] = X[0]; arrx[1] = X[1];
        arrz[0] = Z[0]; arrz[1] = Z[1];
        ll_xz_double(Q0x2, Q0z2, X, Z, a24, M);
        ll_xz_add(arrx + 2, arrz + 2, Q0x2, Q0z2, arrx, arrz, arrx, arrz, M);
        for (j = 2; j < nj; j++)
            ll_xz_add(arrx + 2*j, arrz + 2*j, arrx + 2*(j-1), arrz + 2*(j-1),
                      Q0x2, Q0z2, arrx + 2*(j-2), arrz + 2*(j-2), M);

        ll_xz_mul(Qx, Qz, X, Z, P, a24, M);
        ll_xz_mul(Rx, Rz, Qx, Qz, mmin, a24, M);
        ll_xz_mul(Qdx, Qdz, Qx, Qz, mmin - 1, a24, M);

        for (j = 0; j < 4; j++)
        {
            g[j][0] = M->one[0];
            g[j][1] = M->one[1];
        }

        for (i = mmin; i <= mmax; i++)
        {
            const unsigned char * trow = tab + (i - mmin) * row;

            for (j = 1; j <= maxj; j += 2)
            {
                if (!trow[j])
                    continue;

                MMUL(a, Rx, arrz + 2*(j >> 1));
                MMUL(b, Rz, arrx + 2*(j >> 1));
                MSUB(a, a, b);
                MMUL(g[gi], g[gi], a);
                gi = (gi + 1) & 3;
            }

            sx[0] = Rx[0]; sx[1] = Rx[1];
            sz[0] = Rz[0]; sz[1] = Rz[1];
            ll_xz_add(Rx, Rz, Rx, Rz, Qx, Qz, Qdx, Qdz, M);
            Qdx[0] = sx[0]; Qdx[1] = sx[1];
            Qdz[0] = sz[0]; Qdz[1] = sz[1];
        }

        MMUL(g[0], g[0], g[1]);
        MMUL(g[2], g[2], g[3]);
        MMUL(g[0], g[0], g[2]);

        if (ll_gcd_factor(factor, g[0], n1, n0))
            ret = 1;

        flint_free(arrx);
        flint_free(arrz);
    }

cleanup:
    return ret;
}

int
n_ll_factor_ecm(nn_ptr factor, ulong nhi, ulong nlo, ulong curves,
                ulong B1, ulong B2, flint_rand_t state)
{
    ll_mont_t M;
    unsigned char * tab;
    ulong mmin, mmax, maxj, row, c;
    int ret = 0;

    if (nhi == 0 || (nlo & 1) == 0 || B1 < 2)
        return 0;

    if (B2 < B1)
        B2 = 100 * B1;

    mmin = (B1 + LL_ECM_P / 2) / LL_ECM_P;
    if (mmin == 0)
        mmin = 1;
    mmax = ((B2 - LL_ECM_P / 2) + LL_ECM_P - 1) / LL_ECM_P;
    if (mmax < mmin)
        mmax = mmin;
    maxj = (LL_ECM_P + 1) / 2;
    row = maxj + 1;

    /*
       Mark the (i,j) pairs for which i*P-j or i*P+j is a prime in (B1,B2].
       Walking the primes is cheaper than sieving to B2, and the table is
       built once for all curves.
    */
    tab = flint_calloc((mmax - mmin + 1) * row, sizeof(unsigned char));
    {
        n_primes_t iter;
        ulong p;

        n_primes_init(iter);
        while ((p = n_primes_next(iter)) <= B2)
        {
            ulong i, j;

            if (p <= B1)
                continue;

            i = (p + LL_ECM_P / 2) / LL_ECM_P;   /* nearest multiple of P */
            if (i < mmin || i > mmax)
                continue;
            j = (p > i * LL_ECM_P) ? (p - i * LL_ECM_P) : (i * LL_ECM_P - p);
            if (j <= maxj && (j & 1))
                tab[(i - mmin) * row + j] = 1;
        }
        n_primes_clear(iter);
    }

    ll_mont_init(M, nhi, nlo);

    for (c = 0; c < curves; c++)
    {
        ulong sigma = 6 + n_randint(state, UWORD(1) << 30);

        if (ll_ecm_curve(factor, nhi, nlo, sigma, B1, M, tab, mmin, mmax))
        {
            ret = 1;
            break;
        }
    }

    flint_free(tab);
    return ret;
}

#endif
