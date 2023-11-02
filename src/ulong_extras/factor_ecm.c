/*
    Copyright (C) 2015 Kushagra Singh

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod_vec.h"

static
ulong n_ecm_primorial[] =
{
#ifdef FLINT64

    UWORD(2), UWORD(6), UWORD(30), UWORD(210), UWORD(2310), UWORD(30030),
    UWORD(510510), UWORD(9699690), UWORD(223092870), UWORD(6469693230),
    UWORD(200560490130), UWORD(7420738134810), UWORD(304250263527210),
    UWORD(13082761331670030), UWORD(614889782588491410)

#else

    UWORD(2), UWORD(6), UWORD(30), UWORD(210), UWORD(2310), UWORD(30030),
    UWORD(510510), UWORD(9699690)

#endif
};

#ifdef FLINT64
#define num_n_ecm_primorials 15
#else
#define num_n_ecm_primorials 9
#endif

int
n_factor_ecm(mp_limb_t *f, mp_limb_t curves, mp_limb_t B1, mp_limb_t B2,
             flint_rand_t state, mp_limb_t n)
{
    mp_limb_t P, num, maxD, mmin, mmax, mdiff, prod, maxj, sig;
    int i, j, ret;
    n_ecm_t n_ecm_inf;

    const mp_limb_t *prime_array;

    n_ecm_inf->normbits = flint_clz(n);
    n <<= n_ecm_inf->normbits;
    n_ecm_inf->ninv = n_preinvert_limb(n);
    n_ecm_inf->one = UWORD(1) << n_ecm_inf->normbits;

    ret = 0;

    /************************ STAGE I PRECOMPUTATIONS ************************/

    num = n_prime_pi(B1);   /* number of primes under B1 */

    /* compute list of primes under B1 for stage I */
    prime_array = n_primes_arr_readonly(num);

    /************************ STAGE II PRECOMPUTATIONS ***********************/

    maxD = n_sqrt(B2);

    /* Selecting primorial */

    j = 1;
    while ((j < num_n_ecm_primorials) && (n_ecm_primorial[j] < maxD))
        j += 1;

    P = n_ecm_primorial[j - 1];

    mmin = (B1 + (P/2)) / P;
    mmax = ((B2 - P/2) + P - 1)/P;      /* ceil */
    maxj = (P + 1)/2;
    mdiff = mmax - mmin + 1;

    /* compute GCD_table */

    n_ecm_inf->GCD_table = flint_malloc(maxj + 1);

    for (j = 1; j <= maxj; j += 2)
    {
        if ((j%2) && n_gcd(j, P) == 1)
            n_ecm_inf->GCD_table[j] = 1;
        else
            n_ecm_inf->GCD_table[j] = 0;
    }

    /* compute prime table */

    n_ecm_inf->prime_table = flint_malloc(mdiff * sizeof(unsigned char*));

    for (i = 0; i < mdiff; i++)
        n_ecm_inf->prime_table[i] = flint_malloc((maxj + 1) * sizeof(unsigned char));

    for (i = 0; i < mdiff; i++)
    {
        for (j = 1; j <= maxj; j += 2)
        {
            n_ecm_inf->prime_table[i][j] = 0;

            /* if (i + mmin)*D + j
               is prime, mark 1. Can be possibly prime
               only if gcd(j, D) = 1 */

            if (n_ecm_inf->GCD_table[j] == 1)
            {
                prod = (i + mmin)*P + j;
                if (n_is_prime(prod))
                    n_ecm_inf->prime_table[i][j] = 1;

                prod = (i + mmin)*P - j;
                if (n_is_prime(prod))
                    n_ecm_inf->prime_table[i][j] = 1;
            }
        }
    }

    /****************************** TRY "CURVES" *****************************/

    for (j = 0; j < curves; j++)
    {
        sig = n_randint(state, n >> n_ecm_inf->normbits);
        sig = n_addmod(sig, 7, n >> n_ecm_inf->normbits);
        sig <<= n_ecm_inf->normbits;

        /************************ SELECT CURVE ************************/

        if (n_factor_ecm_select_curve(f, sig, n, n_ecm_inf))
        {
            /* Found factor while selecting curve,
               very very lucky :) */
            (*f) >>= n_ecm_inf->normbits;
            ret = -1;
            goto cleanup;
        }

        /************************** STAGE I ***************************/

        ret = n_factor_ecm_stage_I(f, prime_array, num, B1, n, n_ecm_inf);

        if (ret)
        {
            /* Found factor after stage I */
            (*f) >>= n_ecm_inf->normbits;
            ret = 1;
            goto cleanup;
        }

        /************************** STAGE II ***************************/

        ret = n_factor_ecm_stage_II(f, B1, B2, P, n, n_ecm_inf);

        if (ret)
        {
            /* Found factor after stage II */
            (*f) >>= n_ecm_inf->normbits;
            ret = 2;
            goto cleanup;
        }
    }

    cleanup:

    flint_free(n_ecm_inf->GCD_table);

    for (i = 0; i < mdiff; i++)
        flint_free(n_ecm_inf->prime_table[i]);

    flint_free(n_ecm_inf->prime_table);

    return ret;
}

/* P (x : z) = P1 (x1 : z1) + P2 (x2 : z2) where P0 (x0 : zo) is P - Q */

/*    Coordinates of P :

        x = 4 * z0 * (x1 * x2 - z1 * z2)^2 mod n
        z = 4 * x0 * (x2 * z1 - x1 * z2)^2 mod n
*/

void
n_factor_ecm_add(mp_limb_t *x, mp_limb_t *z, mp_limb_t x1, mp_limb_t z1,
                 mp_limb_t x2, mp_limb_t z2, mp_limb_t x0, mp_limb_t z0,
                 mp_limb_t n, n_ecm_t n_ecm_inf)
{
    mp_limb_t u, v, w;

    if (z1 == 0)
    {
        *x = x2;
        *z = z2;
        return;
    }

    if (z2 == 0)
    {
        *x = x1;
        *z = z1;
        return;
    }

    if (z0 == 0)
    {
        n_factor_ecm_double(x, z, x1, z1, n, n_ecm_inf);
        return;
    }

    u = n_submod(x2, z2, n);        /* u = (x2 - z2) */
    v = n_addmod(x1, z1, n);        /* v = (x1 + z1) */

    /* u = (x2 - z2) * (x1 + z1) */
    u = n_mulmod_preinv(u, v, n, n_ecm_inf->ninv, n_ecm_inf->normbits);

    v = n_submod(x1, z1, n);        /* v = (x1 - z1) */
    w = n_addmod(x2, z2, n);        /* w = (x2 + z2) */

    /* v = (x1 - z1) * (x2 + z2) */
    v = n_mulmod_preinv(v, w, n, n_ecm_inf->ninv, n_ecm_inf->normbits);

    w = n_addmod(u, v, n);          /* w = 2 * (x1 * x2 - z1 * z2) */
    v = n_submod(v, u, n);          /* v = 2 * (x2 * z1 - x1 * z2) */

    /* w = 4 * (x1 * x2 - z1 * z2)^2 */
    w = n_mulmod_preinv(w, w, n, n_ecm_inf->ninv, n_ecm_inf->normbits);

    /* v = 4 * (x2 * z1 - x1 * z2)^2 */
    v = n_mulmod_preinv(v, v, n, n_ecm_inf->ninv, n_ecm_inf->normbits);

    /* x = 4 * z0 * (x1 * x2 - z1 * z2)^2 */
    *x = n_mulmod_preinv(z0, w, n, n_ecm_inf->ninv, n_ecm_inf->normbits);

    /* z = 4 * x0 * (x2 * z1 - x1 * z2)^2 */
    *z = n_mulmod_preinv(x0, v, n, n_ecm_inf->ninv, n_ecm_inf->normbits);
}

/* P (x : z) = 2 * P1 (x0 : z0)  */

/*
    Coordinates of P :

        x = (x0 + z0)^2 * (x0 - z0)^2 mod n
        z = 4 * x0 * z0 * ((x0 - z0)^2 + a24 * 4 * x0 * z0) mod n
*/

/* a24 = (a + 2) / 4 mod n */

void
n_factor_ecm_double(mp_limb_t *x, mp_limb_t *z, mp_limb_t x0, mp_limb_t z0,
                    mp_limb_t n, n_ecm_t n_ecm_inf)
{
    mp_limb_t u, v, w;

    if (z0 == 0)
    {
        *x = x0;
        *z = 0;
        return;
    }

    u = n_addmod(x0, z0, n);
    u = n_mulmod_preinv(u, u, n, n_ecm_inf->ninv, n_ecm_inf->normbits);
    v = n_submod(x0, z0, n);
    v = n_mulmod_preinv(v, v, n, n_ecm_inf->ninv, n_ecm_inf->normbits);
    *x = n_mulmod_preinv(u, v, n, n_ecm_inf->ninv, n_ecm_inf->normbits);
    w = n_submod(u, v, n);
    u = n_mulmod_preinv(w, n_ecm_inf->a24, n, n_ecm_inf->ninv, n_ecm_inf->normbits);
    u = n_addmod(u, v, n);
    *z = n_mulmod_preinv(w, u, n, n_ecm_inf->ninv, n_ecm_inf->normbits);
}

/* P (x0 : z0) <- kP using Montgomery ladder algorithm */

void
n_factor_ecm_mul_montgomery_ladder(mp_limb_t *x, mp_limb_t *z, mp_limb_t x0, mp_limb_t z0,
                                   mp_limb_t k, mp_limb_t n, n_ecm_t n_ecm_inf)
{
    mp_limb_t x1, z1, x2, z2, len;      /* Q (x1 : z1), P (x2 : z2) */

    if (k == 0)
    {
        *x = 0;
        *z = 0;
        return;
    }

    if (k == 1)
    {
        *x = x0;
        *z = z0;
        return;
    }

    x1 = x0;    /* Q <- P0 */
    z1 = z0;
    x2 = 0;
    z2 = 0;

    n_factor_ecm_double(&x2, &z2, x0, z0, n, n_ecm_inf);    /* P <- 2P0 */

    len = n_sizeinbase(k, 2) - 2;

    while (1)
    {
        if (((UWORD(1) << len) & k) != 0)       /* ith bit is 1 */
        {
            /* Q <- P + Q */
            n_factor_ecm_add(&x1, &z1, x1, z1, x2, z2, x0, z0, n, n_ecm_inf);
            /* P <- 2 * P */
            n_factor_ecm_double(&x2, &z2, x2, z2, n, n_ecm_inf);
        }
        else
        {
            /* P <- P + Q */
            n_factor_ecm_add(&x2, &z2, x1, z1, x2, z2, x0, z0, n, n_ecm_inf);
            /* Q <- 2 * Q */
            n_factor_ecm_double(&x1, &z1, x1, z1, n, n_ecm_inf);
        }


        if (len == 0)
            break;
        else
            len -= 1;
    }

    *x = x1;
    *z = z1;
}

int
n_factor_ecm_select_curve(mp_limb_t *f, mp_limb_t sig, mp_limb_t n, n_ecm_t n_ecm_inf)
{
    mp_limb_t u, v, w, t, hi, lo;
    mp_ptr a;
    int ret = 0;
    TMP_INIT;

    TMP_START;
    a = TMP_ALLOC(2 * sizeof(mp_limb_t));

    u = sig;

    /* v = sig * 4 */
    v = n_mulmod_preinv(u, UWORD(4) << n_ecm_inf->normbits, n, n_ecm_inf->ninv,
                        n_ecm_inf->normbits);

    /* w = sig ^ 2 */
    w = n_mulmod_preinv(u, u, n, n_ecm_inf->ninv, n_ecm_inf->normbits);

    u = n_submod(w, UWORD(5) << n_ecm_inf->normbits, n);  /* u = sig^2 - 5 */

    /* w = u * u */
    w = n_mulmod_preinv(u, u, n, n_ecm_inf->ninv, n_ecm_inf->normbits);

    /* x = u * u * u */
    n_ecm_inf->x = n_mulmod_preinv(w, u, n, n_ecm_inf->ninv, n_ecm_inf->normbits);

    /* w = v * v */
    w = n_mulmod_preinv(v, v, n, n_ecm_inf->ninv, n_ecm_inf->normbits);

    /* z = v * v * v */
    n_ecm_inf->z = n_mulmod_preinv(w, v, n, n_ecm_inf->ninv, n_ecm_inf->normbits);

    w = n_mulmod_preinv(n_ecm_inf->x, v, n, n_ecm_inf->ninv, n_ecm_inf->normbits);
    t = n_mulmod_preinv(w, UWORD(4) << n_ecm_inf->normbits, n, n_ecm_inf->ninv,
                        n_ecm_inf->normbits);
    w = n_mulmod_preinv(u, UWORD(3) << n_ecm_inf->normbits, n, n_ecm_inf->ninv,
                        n_ecm_inf->normbits);

    u = n_submod(v, u, n);
    v = n_addmod(v, w, n);
    w = n_mulmod_preinv(u, u, n, n_ecm_inf->ninv, n_ecm_inf->normbits);
    u = n_mulmod_preinv(u, w, n, n_ecm_inf->ninv, n_ecm_inf->normbits);

    n_ecm_inf->a24 = n_mulmod_preinv(u, v, n, n_ecm_inf->ninv, n_ecm_inf->normbits);
    v = n_mulmod_preinv(t, n_ecm_inf->z, n, n_ecm_inf->ninv, n_ecm_inf->normbits);

    *f = n_gcdinv(&u, v, n);

    if (*f == n)
    {
        goto cleanup;
    }
    else if (*f > n_ecm_inf->one)
    {
        ret = 1;
        goto cleanup;
    }

    a[1] = UWORD(0);
    a[0] = u;
    mpn_lshift(a, a, 2, n_ecm_inf->normbits);       /* shifting back right */
    hi = a[1];
    lo = a[0];
    u = n_ll_mod_preinv(hi, lo, n, n_ecm_inf->ninv);

    v = n_mulmod_preinv(u, t, n, n_ecm_inf->ninv, n_ecm_inf->normbits);
    n_ecm_inf->x = n_mulmod_preinv(n_ecm_inf->x, v, n, n_ecm_inf->ninv,
                                   n_ecm_inf->normbits);

    v = n_mulmod_preinv(u, n_ecm_inf->z, n, n_ecm_inf->ninv, n_ecm_inf->normbits);
    w = n_mulmod_preinv(n_ecm_inf->a24, v, n, n_ecm_inf->ninv, n_ecm_inf->normbits);
    n_ecm_inf->a24 = n_addmod(w, UWORD(2) << n_ecm_inf->normbits, n);
    n_ecm_inf->a24 >>= 2;
    n_ecm_inf->a24 >>= n_ecm_inf->normbits;
    n_ecm_inf->a24 <<= n_ecm_inf->normbits;
    n_ecm_inf->z = n_ecm_inf->one;

cleanup:
    TMP_END;

    return ret;
}

int
n_factor_ecm_stage_I(mp_limb_t *f, const mp_limb_t *prime_array, mp_limb_t num,
                     mp_limb_t B1, mp_limb_t n, n_ecm_t n_ecm_inf)
{
    mp_limb_t times;
    int i, j, p;

    for (i = 0; i < num; i++)
    {
        p = n_flog(B1, prime_array[i]);
        times = prime_array[i];

        for (j = 1; j <= p; j ++)
        {
            n_factor_ecm_mul_montgomery_ladder(&(n_ecm_inf->x), &(n_ecm_inf->z),
                                               n_ecm_inf->x, n_ecm_inf->z,
                                               times, n, n_ecm_inf);
        }

        *f = n_gcd(n_ecm_inf->z, n);

        if ((*f > n_ecm_inf->one) && (*f < n))
        {
            /* Found factor in stage I */
            return 1;
        }
    }

    return 0;
}

int
n_factor_ecm_stage_II(mp_limb_t *f, mp_limb_t B1, mp_limb_t B2, mp_limb_t P,
                      mp_limb_t n, n_ecm_t n_ecm_inf)
{

    mp_limb_t g, Qx, Qz, Rx, Rz, Qdx, Qdz, a, b;
    mp_limb_t mmin, mmax, maxj, Q0x2, Q0z2;
    int i, j, ret;
    mp_ptr arrx, arrz;

    mmin = (B1 + (P/2)) / P;
    mmax = ((B2 - P/2) + P - 1)/P;      /* ceil */
    maxj = (P + 1)/2;

    g = n_ecm_inf->one;

    arrx = _nmod_vec_init((maxj >> 1) + 1);
    arrz = _nmod_vec_init((maxj >> 1) + 1);

    ret = 0;

    /* arr[0] = Q0 */
    arrx[0] = n_ecm_inf->x;
    arrz[0] = n_ecm_inf->z;

    /* Q0x2, Q0z2 = 2Q0 */
    n_factor_ecm_double(&Q0x2, &Q0z2, arrx[0], arrz[0], n, n_ecm_inf);

    /* arr[1] = 3Q0 */
    n_factor_ecm_add(arrx + 1, arrz + 1, Q0x2, Q0z2, arrx[0], arrz[0],
                     arrx[0], arrz[0], n, n_ecm_inf);

    /* For each odd j (j > 3) , compute j * Q0 [x0 :: z0] */
    /* jth stored in arr[j/2] */

    /* We are adding 2Q0 every time. Need to calculate all j's
       as (j - 2)Q0 is required for (j + 2)Q0 */

    for (j = 2; j <= (maxj >> 1); j += 1)
    {
        /* jQ0 = (j - 2)Q0 + 2Q0
           Difference is (j - 4)Q0 */

        n_factor_ecm_add(arrx + j, arrz + j, arrx[j - 1], arrz[j - 1], Q0x2,
                         Q0z2, arrx[j - 2], arrz[j - 2], n, n_ecm_inf);
    }

    /* Q = D * Q_0 */
    n_factor_ecm_mul_montgomery_ladder(&Qx, &Qz, n_ecm_inf->x, n_ecm_inf->z, P, n, n_ecm_inf);

    /* R = mmin * Q */
    n_factor_ecm_mul_montgomery_ladder(&Rx, &Rz, Qx, Qz, mmin, n, n_ecm_inf);

    /* Qd = (mmin - 1) * Q */
    n_factor_ecm_mul_montgomery_ladder(&Qdx, &Qdz, Qx, Qz, mmin - 1, n, n_ecm_inf);

    /* main stage II step */

    for (i = mmin; i <= mmax; i ++)
    {
        for (j = 1; j <= maxj; j += 2)
        {
            if (n_ecm_inf->prime_table[i - mmin][j] == 1)
            {
                a = n_mulmod_preinv(Rx, arrz[j >> 1], n, n_ecm_inf->ninv, n_ecm_inf->normbits);
                b = n_mulmod_preinv(Rz, arrx[j >> 1], n, n_ecm_inf->ninv, n_ecm_inf->normbits);
                a = n_submod(a, b, n);
                g = n_mulmod_preinv(g, a, n, n_ecm_inf->ninv, n_ecm_inf->normbits);
            }
        }

        a = Rx;
        b = Rz;

        /* R = R + Q
           difference is stored in Qd, initially (Mmin - 1)Q */

        n_factor_ecm_add(&Rx, &Rz, Rx, Rz, Qx, Qz, Qdx, Qdz, n, n_ecm_inf);

        Qdx = a;
        Qdz = b;
    }

    *f = n_gcd(g, n);

    if ((*f > n_ecm_inf->one) && (*f < n))
    {
        /* Found factor in stage I */
        ret = 1;
    }

    _nmod_vec_clear(arrx);
    _nmod_vec_clear(arrz);

    return ret;
}
