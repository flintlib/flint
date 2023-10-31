/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gmpcompat.h"
#include "nmod.h"
#include "nmod_vec.h"

TEST_FUNCTION_START(nmod_vec_nmod, state)
{
    int i;

    /* nmod_add */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
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

        flint_mpz_set_ui(x, a);
        flint_mpz_set_ui(y, b);
        mpz_add(z, x, y);
        flint_mpz_mod_ui(z, z, m);

        if (flint_mpz_cmp_ui(z, c) != 0)
        {
            flint_printf("FAIL (add):\n");
            flint_printf("m = %wu\n", m);
            fflush(stdout);
            flint_abort();
        }

        mpz_clear(x);
        mpz_clear(y);
        mpz_clear(z);
    }

    /* nmod_sub */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
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

        flint_mpz_set_ui(x, a);
        flint_mpz_set_ui(y, b);
        mpz_sub(z, x, y);
        flint_mpz_mod_ui(z, z, m);

        if (flint_mpz_cmp_ui(z, c) != 0)
        {
            flint_printf("FAIL (sub):\n");
            flint_printf("m = %wu\n", m);
            fflush(stdout);
            flint_abort();
        }

        mpz_clear(x);
        mpz_clear(y);
        mpz_clear(z);
    }

    /* nmod_mul */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
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

        flint_mpz_set_ui(x, a);
        flint_mpz_set_ui(y, b);
        mpz_mul(z, x, y);
        flint_mpz_mod_ui(z, z, m);

        if (flint_mpz_cmp_ui(z, c) != 0)
        {
            flint_printf("FAIL (mul):\n");
            flint_printf("m = %wu\n", m);
            fflush(stdout);
            flint_abort();
        }

        mpz_clear(x);
        mpz_clear(y);
        mpz_clear(z);
    }

    /* nmod_div */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
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

        flint_mpz_set_ui(x, a);
        flint_mpz_set_ui(y, b);
        flint_mpz_set_ui(z, m);
        mpz_invert(z, y, z);
        mpz_mul(z, x, z);
        flint_mpz_mod_ui(z, z, m);

        if (flint_mpz_cmp_ui(z, c) != 0)
        {
            flint_printf("FAIL (div):\n");
            flint_printf("m = %wu\n", m);
            fflush(stdout);
            flint_abort();
        }

        mpz_clear(x);
        mpz_clear(y);
        mpz_clear(z);
    }

    /* nmod_inv */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
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

        flint_mpz_set_ui(y, b);
        flint_mpz_set_ui(z, m);
        mpz_invert(z, y, z);

        if (flint_mpz_cmp_ui(z, c) != 0)
        {
            flint_printf("FAIL (div):\n");
            flint_printf("m = %wu\n", m);
            fflush(stdout);
            flint_abort();
        }

        mpz_clear(y);
        mpz_clear(z);
    }

    /* nmod_pow_ui */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
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

        flint_mpz_set_ui(y, b);
        flint_mpz_set_ui(z, m);
        flint_mpz_powm_ui(z, y, exp, z);

        if (flint_mpz_cmp_ui(z, c) != 0)
        {
            flint_printf("FAIL (pow):\n");
            flint_printf("m = %wu\n", m);
            fflush(stdout);
            flint_abort();
        }

        mpz_clear(y);
        mpz_clear(z);
    }

    TEST_FUNCTION_END(state);
}
