/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "ulong_extras.h"

#define FLINT_MU_LOOKUP_CUTOFF 1024

#if FLINT64
const mp_limb_t FLINT_MOEBIUS_ODD[] = 
{
    UWORD(0x4289108a05208102), UWORD(0x19988004a8a12422), UWORD(0x1a8245028906a062),
    UWORD(0x229428012aa26a00), UWORD(0x8422a98980440a18), 0x224925084068929aUL,
    UWORD(0xa1200942a8980a84), UWORD(0x8622622216a00428), UWORD(0x6060829286a590a9),
    UWORD(0x5a2190081420a1a8), UWORD(0x8a92a284a8018200), UWORD(0x980a2602491a2248),
    UWORD(0x8106a08982a26848), UWORD(0xa60085a6004a5919), UWORD(0x88a188245a221228),
    UWORD(0x0108884a22186025)
};
#else
const mp_limb_t FLINT_MOEBIUS_ODD[] =
{
    UWORD(0x05208102), 0x4289108aUL, UWORD(0xa8a12422), UWORD(0x19988004), UWORD(0x8906a062),
    UWORD(0x1a824502), UWORD(0x2aa26a00), UWORD(0x22942801), UWORD(0x80440a18), UWORD(0x8422a989),
    0x4068929aUL, UWORD(0x22492508), UWORD(0xa8980a84), UWORD(0xa1200942), UWORD(0x16a00428),
    UWORD(0x86226222), UWORD(0x86a590a9), UWORD(0x60608292), UWORD(0x1420a1a8), UWORD(0x5a219008),
    UWORD(0xa8018200), UWORD(0x8a92a284), UWORD(0x491a2248), UWORD(0x980a2602), UWORD(0x82a26848),
    UWORD(0x8106a089), UWORD(0x004a5919), UWORD(0xa60085a6), UWORD(0x5a221228), UWORD(0x88a18824),
    UWORD(0x22186025), 0x0108884aUL
};
#endif

void n_moebius_mu_vec(int * mu, ulong len)
{
    slong k;
    ulong pi;
    const mp_limb_t * primes;
    mp_limb_t p, q;

    pi = n_prime_pi(len);
    primes = n_primes_arr_readonly(pi);

    if (len)
        mu[0] = 0;
    for (k = 1; k < len; k++)
        mu[k] = 1;
    for (k = 0; k < pi; k++)
    {
        p = primes[k];
        for (q = p; q < len; q += p)
            mu[q] = -mu[q];
        p = p * p;
        for (q = p; q < len; q += p)
            mu[q] = 0;
    }
}

int n_moebius_mu(mp_limb_t n)
{
    int i;
    n_factor_t fac;

    if (n % 2 == 0)
    {
        if (n % 4 == 0)
            return 0;
        return -n_moebius_mu(n / 2);
    }

    if (n < FLINT_MU_LOOKUP_CUTOFF)
    {
        mp_limb_t m;
        n -= 1;
        m = FLINT_MOEBIUS_ODD[n / FLINT_BITS];
        m &= (UWORD(3) << (n % FLINT_BITS));
        m >>= (n % FLINT_BITS);
        return ((int) m) - 1;
    }

    /* Weed out just a few more common squares first */
    if (n % 9 == 0 || n % 25 == 0)
        return 0;

    n_factor_init(&fac);
    n_factor(&fac, n, 1);
    for (i = 0; i < fac.num; i++)
    {
        if (fac.exp[i] != 1)
            return 0;
    }

    if (fac.num % 2)
        return -1;
    else
        return 1;
}
