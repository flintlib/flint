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

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    flint_rand_t state;
    flint_randinit(state);

    printf("nmod....");
    fflush(stdout);

    /* nmod_add */
    for (i = 0; i < 10000; i++)
    {
        nmod_t mod;
        mp_limb_t m, a, b, c;
        mpz_t x, y, z;

        m = n_randtest_not_zero(state);

        nmod_init(&mod, m);
        a = n_randlimb(state) % m;
        b = n_randlimb(state) % m;

        c = nmod_add(a, b, mod);

        mpz_init(x);
        mpz_init(y);
        mpz_init(z);

        mpz_set_ui(x, a);
        mpz_set_ui(y, b);
        mpz_add(z, x, y);
        mpz_mod_ui(z, z, m);

        if (mpz_cmp_ui(z, c) != 0)
        {
            printf("FAIL (add):\n");
            printf("m = %lu\n", m);
            abort();
        }

        mpz_clear(x);
        mpz_clear(y);
        mpz_clear(z);
    }

    /* nmod_sub */
    for (i = 0; i < 10000; i++)
    {
        nmod_t mod;
        mp_limb_t m, a, b, c;
        mpz_t x, y, z;

        m = n_randtest_not_zero(state);

        nmod_init(&mod, m);
        a = n_randlimb(state) % m;
        b = n_randlimb(state) % m;

        c = nmod_sub(a, b, mod);

        mpz_init(x);
        mpz_init(y);
        mpz_init(z);

        mpz_set_ui(x, a);
        mpz_set_ui(y, b);
        mpz_sub(z, x, y);
        mpz_mod_ui(z, z, m);

        if (mpz_cmp_ui(z, c) != 0)
        {
            printf("FAIL (sub):\n");
            printf("m = %lu\n", m);
            abort();
        }

        mpz_clear(x);
        mpz_clear(y);
        mpz_clear(z);
    }

    /* nmod_mul */
    for (i = 0; i < 10000; i++)
    {
        nmod_t mod;
        mp_limb_t m, a, b, c;
        mpz_t x, y, z;

        m = n_randtest_not_zero(state);

        nmod_init(&mod, m);
        a = n_randlimb(state) % m;
        b = n_randlimb(state) % m;

        c = nmod_mul(a, b, mod);

        mpz_init(x);
        mpz_init(y);
        mpz_init(z);

        mpz_set_ui(x, a);
        mpz_set_ui(y, b);
        mpz_mul(z, x, y);
        mpz_mod_ui(z, z, m);

        if (mpz_cmp_ui(z, c) != 0)
        {
            printf("FAIL (mul):\n");
            printf("m = %lu\n", m);
            abort();
        }

        mpz_clear(x);
        mpz_clear(y);
        mpz_clear(z);
    }

    /* nmod_div */
    for (i = 0; i < 10000; i++)
    {
        nmod_t mod;
        mp_limb_t m, a, b, c;
        mpz_t x, y, z;

        m = n_randtest_prime(state, 0);

        nmod_init(&mod, m);

        a = n_randlimb(state) % m;
        do { b = n_randlimb(state) % m; } while (b == 0);

        c = nmod_div(a, b, mod);

        mpz_init(x);
        mpz_init(y);
        mpz_init(z);

        mpz_set_ui(x, a);
        mpz_set_ui(y, b);
        mpz_set_ui(z, m);
        mpz_invert(z, y, z);
        mpz_mul(z, x, z);
        mpz_mod_ui(z, z, m);

        if (mpz_cmp_ui(z, c) != 0)
        {
            printf("FAIL (div):\n");
            printf("m = %lu\n", m);
            abort();
        }

        mpz_clear(x);
        mpz_clear(y);
        mpz_clear(z);
    }

    /* nmod_inv */
    for (i = 0; i < 10000; i++)
    {
        nmod_t mod;
        mp_limb_t m, b, c;
        mpz_t y, z;

        m = n_randtest_prime(state, 0);

        nmod_init(&mod, m);

        do { b = n_randlimb(state) % m; } while (b == 0);

        c = nmod_inv(b, mod);

        mpz_init(y);
        mpz_init(z);

        mpz_set_ui(y, b);
        mpz_set_ui(z, m);
        mpz_invert(z, y, z);

        if (mpz_cmp_ui(z, c) != 0)
        {
            printf("FAIL (div):\n");
            printf("m = %lu\n", m);
            abort();
        }

        mpz_clear(y);
        mpz_clear(z);
    }

    /* nmod_pow_ui */
    for (i = 0; i < 10000; i++)
    {
        nmod_t mod;
        mp_limb_t m, b, c;
        mpz_t y, z;
        ulong exp;

        m = n_randtest_prime(state, 0);
        exp = n_randtest(state);

        nmod_init(&mod, m);

        b = n_randlimb(state) % m;

        c = nmod_pow_ui(b, exp, mod);

        mpz_init(y);
        mpz_init(z);

        mpz_set_ui(y, b);
        mpz_set_ui(z, m);
        mpz_powm_ui(z, y, exp, z);

        if (mpz_cmp_ui(z, c) != 0)
        {
            printf("FAIL (pow):\n");
            printf("m = %lu\n", m);
            abort();
        }

        mpz_clear(y);
        mpz_clear(z);
    }


    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
