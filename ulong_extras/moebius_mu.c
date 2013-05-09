/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include "flint.h"
#include "ulong_extras.h"

#define FLINT_MU_LOOKUP_CUTOFF 1024

#if FLINT64
const mp_limb_t FLINT_MOEBIUS_ODD[] = 
{
    0x4289108a05208102UL, 0x19988004a8a12422UL, 0x1a8245028906a062UL,
    0x229428012aa26a00UL, 0x8422a98980440a18UL, 0x224925084068929aUL,
    0xa1200942a8980a84UL, 0x8622622216a00428UL, 0x6060829286a590a9UL,
    0x5a2190081420a1a8UL, 0x8a92a284a8018200UL, 0x980a2602491a2248UL,
    0x8106a08982a26848UL, 0xa60085a6004a5919UL, 0x88a188245a221228UL,
    0x0108884a22186025UL
};
#else
const mp_limb_t FLINT_MOEBIUS_ODD[] =
{
    0x05208102UL, 0x4289108aUL, 0xa8a12422UL, 0x19988004UL, 0x8906a062UL,
    0x1a824502UL, 0x2aa26a00UL, 0x22942801UL, 0x80440a18UL, 0x8422a989UL,
    0x4068929aUL, 0x22492508UL, 0xa8980a84UL, 0xa1200942UL, 0x16a00428UL,
    0x86226222UL, 0x86a590a9UL, 0x60608292UL, 0x1420a1a8UL, 0x5a219008UL,
    0xa8018200UL, 0x8a92a284UL, 0x491a2248UL, 0x980a2602UL, 0x82a26848UL,
    0x8106a089UL, 0x004a5919UL, 0xa60085a6UL, 0x5a221228UL, 0x88a18824UL,
    0x22186025UL, 0x0108884aUL
};
#endif

void n_moebius_mu_vec(int * mu, ulong len)
{
    len_t k;
    ulong pi;

    mp_limb_t p, q;

    pi = n_prime_pi(len);
    n_compute_primes(pi);

    if (len)
        mu[0] = 0;
    for (k = 1; k < len; k++)
        mu[k] = 1;
    for (k = 0; k < pi; k++)
    {
        p = flint_primes[k];
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
        m &= (3UL << (n % FLINT_BITS));
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
