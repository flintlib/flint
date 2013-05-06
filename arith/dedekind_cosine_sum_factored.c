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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include "flint.h"
#include "arith.h"
#include "ulong_extras.h"

static const int mod4_tab[8] = { 2, 1, 3, 0, 0, 3, 1, 2 };

static const int gcd24_tab[24] = {
    24, 1, 2, 3, 4, 1, 6, 1, 8, 3, 2, 1,
    12, 1, 2, 3, 8, 1, 6, 1, 4, 3, 2, 1
};

static mp_limb_t
n_sqrtmod_2exp(mp_limb_t a, int k)
{
    mp_limb_t x;
    int i;

    if (a == 0 || k == 0)
        return 0;

    if (k == 1)
        return 1;

    if (k == 2)
    {
        if (a == 1)
            return 1;
        return 0;
    }

    x = 1;
    for (i = 3; i < k; i++)
        x += (a - x * x) / 2;

    if (k < FLINT_BITS)
        x &= ((1UL << k) - 1);

    return x;
}

static mp_limb_t
n_sqrtmod_ppow(mp_limb_t a, mp_limb_t p, int k, mp_limb_t pk, mp_limb_t pkinv)
{
    mp_limb_t r, t;
    int i;

    r = n_sqrtmod(a, p);
    if (r == 0)
        return r;

    i = 1;
    while (i < k)
    {
        t = n_mulmod2_preinv(r, r, pk, pkinv);
        t = n_submod(t, a, pk);
        t = n_mulmod2_preinv(t, n_invmod(n_addmod(r, r, pk), pk), pk, pkinv);
        r = n_submod(r, t, pk);
        i *= 2;
    }

    return r;
}

void
trigprod_mul_prime_power(trig_prod_t prod, mp_limb_t k, mp_limb_t n,
                                                mp_limb_t p, int exp)
{
    mp_limb_t m, mod, inv;

    if (k <= 3)
    {
        if (k == 0)
        {
            prod->prefactor = 0;
        }
        else if (k == 2 && (n % 2 == 1))
        {
            prod->prefactor *= -1;
        }
        else if (k == 3)
        {
            switch (n % 3)
            {
                case 0:
                    prod->prefactor *= 2;
                    prod->cos_p[prod->n] = 1;
                    prod->cos_q[prod->n] = 18;
                    break;
                case 1:
                    prod->prefactor *= -2;
                    prod->cos_p[prod->n] = 7;
                    prod->cos_q[prod->n] = 18;
                    break;
                case 2:
                    prod->prefactor *= -2;
                    prod->cos_p[prod->n] = 5;
                    prod->cos_q[prod->n] = 18;
                    break;
            }
            prod->n++;
        }
        return;
    }

    /* Power of 2 */
    if (p == 2)
    {
        mod = 8 * k;
        inv = n_preinvert_limb(mod);

        m = n_submod(1, n_mod2_preinv(24 * n, mod, inv), mod);
        m = n_sqrtmod_2exp(m, exp + 3);
        m = n_mulmod2_preinv(m, n_invmod(3, mod), mod, inv);

        prod->prefactor *= n_jacobi(-1, m);
        if (exp % 2 == 1)
            prod->prefactor *= -1;
        prod->sqrt_p *= k;
        prod->cos_p[prod->n] = (mp_limb_signed_t)(k - m);
        prod->cos_q[prod->n] = 2 * k;
        prod->n++;
        return;
    }

    /* Power of 3 */
    if (p == 3)
    {
        mod = 3 * k;
        inv = n_preinvert_limb(mod);

        m = n_submod(1, n_mod2_preinv(24 * n, mod, inv), mod);
        m = n_sqrtmod_ppow(m, p, exp + 1, mod, inv);
        m = n_mulmod2_preinv(m, n_invmod(8, mod), mod, inv);

        prod->prefactor *= (2 * n_jacobi_unsigned(m, 3));
        if (exp % 2 == 0)
            prod->prefactor *= -1;
        prod->sqrt_p *= k;
        prod->sqrt_q *= 3;
        prod->cos_p[prod->n] = (mp_limb_signed_t)(3 * k - 8 * m);
        prod->cos_q[prod->n] = 6 * k;
        prod->n++;
        return;
    }

    /* Power of prime greater than 3 */
    inv = n_preinvert_limb(k);
    m = n_submod(1, n_mod2_preinv(24 * n, k, inv), k);

    if (m % p == 0)
    {
        if (exp == 1)
        {
            prod->prefactor *= n_jacobi(3, k);
            prod->sqrt_p *= k;
        }
        else
            prod->prefactor = 0;
        return;
    }

    m = n_sqrtmod_ppow(m, p, exp, k, inv);

    if (m == 0)
    {
        prod->prefactor = 0;
        return;
    }

    prod->prefactor *= 2;
    prod->prefactor *= n_jacobi(3, k);
    prod->sqrt_p *= k;
    prod->cos_p[prod->n] = 4 * n_mulmod2_preinv(m, n_invmod(24, k), k, inv);
    prod->cos_q[prod->n] = k;
    prod->n++;
}

/*
Solve (k2^2 * d2 * e) * n1 = (d2 * e * n + (k2^2 - 1) / d1)   mod k2

TODO: test this on 32 bit
*/
static mp_limb_t
solve_n1(mp_limb_t n, mp_limb_t k1, mp_limb_t k2,
        mp_limb_t d1, mp_limb_t d2, mp_limb_t e)
{
    mp_limb_t inv, n1, u, t[2];

    inv = n_preinvert_limb(k1);

    umul_ppmm(t[1], t[0], k2, k2);
    sub_ddmmss(t[1], t[0], t[1], t[0], 0UL, 1UL);
    mpn_divrem_1(t, 0, t, 2, d1);

    n1 = n_ll_mod_preinv(t[1], t[0], k1, inv);
    n1 = n_mod2_preinv(n1 + d2*e*n, k1, inv);

    u = n_mulmod2_preinv(k2, k2, k1, inv);
    u = n_invmod(u * d2 * e, k1);
    n1 = n_mulmod2_preinv(n1, u, k1, inv);

    return n1;
}


void
arith_hrr_expsum_factored(trig_prod_t prod, mp_limb_t k, mp_limb_t n)
{
    n_factor_t fac;
    int i;

    if (k <= 1)
    {
        prod->prefactor = k;
        return;
    }

    n_factor_init(&fac);
    n_factor(&fac, k, 0);

    /* Repeatedly factor A_k(n) into A_k1(n1)*A_k2(n2) with k1, k2 coprime */
    for (i = 0; i + 1 < fac.num && prod->prefactor != 0; i++)
    {
        mp_limb_t p, k1, k2, inv, n1, n2;

        p = fac.p[i];

        /* k = 2 * k1 with k1 odd */
        if (p == 2UL && fac.exp[i] == 1)
        {
            k2 = k / 2;
            inv = n_preinvert_limb(k2);

            n2 = n_invmod(32, k2);
            n2 = n_mulmod2_preinv(n2,
                    n_mod2_preinv(8*n + 1, k2, inv), k2, inv);
            n1 = ((k2 % 8 == 3) || (k2 % 8 == 5)) ^ (n & 1);

            trigprod_mul_prime_power(prod, 2, n1, 2, 1);
            k = k2;
            n = n2;
        }
        /* k = 4 * k1 with k1 odd */
        else if (p == 2UL && fac.exp[i] == 2)
        {
            k2 = k / 4;
            inv = n_preinvert_limb(k2);

            n2 = n_invmod(128, k2);
            n2 = n_mulmod2_preinv(n2,
                n_mod2_preinv(8*n + 5, k2, inv), k2, inv);
            n1 = (n + mod4_tab[(k2 / 2) % 8]) % 4;

            trigprod_mul_prime_power(prod, 4, n1, 2, 2);
            prod->prefactor *= -1;
            k = k2;
            n = n2;
        }
        /* k = k1 * k2 with k1 odd or divisible by 8 */
        else
        {
            mp_limb_t d1, d2, e;

            k1 = n_pow(fac.p[i], fac.exp[i]);
            k2 = k / k1;

            d1 = gcd24_tab[k1 % 24];
            d2 = gcd24_tab[k2 % 24];
            e = 24 / (d1 * d2);

            n1 = solve_n1(n, k1, k2, d1, d2, e);
            n2 = solve_n1(n, k2, k1, d2, d1, e);

            trigprod_mul_prime_power(prod, k1, n1, fac.p[i], fac.exp[i]);
            k = k2;
            n = n2;
        }
    }

    if (fac.num != 0 && prod->prefactor != 0)
        trigprod_mul_prime_power(prod, k, n,
            fac.p[fac.num - 1], fac.exp[fac.num - 1]);

}
