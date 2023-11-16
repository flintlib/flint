/*
    Copyright (C) 2008, 2009 David Harvey
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz_vec.h"
#include "fmpz_extras.h"
#include "arith.h"
#include "bernoulli.h"

/* test this internal function which should really be in FLINT */
ulong _bernoulli_n_muldivrem_precomp(ulong * q, ulong a, ulong b, ulong n, double bnpre);
ulong _bernoulli_mod_p_harvey_powg(ulong p, ulong pinv, ulong k);
ulong _bernoulli_mod_p_harvey_pow2(ulong p, ulong pinv, ulong k);

void test_bern_modp_pow2(ulong p, ulong k)
{
    ulong x, y;
    ulong pinv = n_preinvert_limb(p);

    if (n_powmod2_preinv(2, k, p, pinv) == 1)
        return;

    x = _bernoulli_mod_p_harvey_powg(p, pinv, k);
    y = _bernoulli_mod_p_harvey_pow2(p, pinv, k);

    if (x != y)
    {
        flint_printf("FAIL\n");
        flint_printf("p = %wu\n", p);
        flint_printf("k = %wu\n", k);
        flint_printf("x = %wu\n", x);
        flint_printf("y = %wu\n", y);
        flint_abort();
    }
}

TEST_FUNCTION_START(bernoulli_mod_p_harvey, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ulong a, b, n, q1, r1, r2;
        mp_limb_t q2[2];
        double bnpre;

        a = n_randtest_bits(state, FLINT_D_BITS);
        b = n_randtest_bits(state, FLINT_D_BITS);
        do {
            n = n_randtest_bits(state, FLINT_D_BITS);
        } while (n == 0);

        a = a % n;
        b = b % n;

        bnpre = (double) b / (double) n;
        r1 = _bernoulli_n_muldivrem_precomp(&q1, a, b, n, bnpre);

        umul_ppmm(q2[1], q2[0], a, b);
        r2 = mpn_divrem_1(q2, 0, q2, 2, n);

        if (q2[0] != q1 || r1 != r2)
        {
            flint_printf("FAIL\n");
            flint_printf("a = %wu\n", a);
            flint_printf("b = %wu\n", b);
            flint_printf("n = %wu\n", n);
            flint_printf("q1 = %wu\n", q1);
            flint_printf("r1 = %wu\n", r1);
            flint_printf("q2 = %wu\n", q2);
            flint_printf("r2 = %wu\n", r2);
            flint_abort();
        }
    }

    {
        slong n, N, iter;
        ulong x, y, z, p, pinv;
        fmpz * num;
        fmpz * den;

        N = 300;

        num = _fmpz_vec_init(N);
        den = _fmpz_vec_init(N);
        _arith_bernoulli_number_vec_recursive(num, den, N);

        for (n = 2; n < N; n += 2)
        {
            for (iter = 0; iter < 10 * 0.1 * flint_test_multiplier(); iter++)
            {
                if (n_randint(state, 100) == 0)
                {
                    p = 1000003;
                }
                else
                {
                    p = n + 2 + n_randint(state, 2 * N);
                    p = n_nextprime(p, 1);
                }

                pinv = n_preinvert_limb(p);

                x = _bernoulli_mod_p_harvey_powg(p, pinv, n);

                y = fmpz_fdiv_ui(num + n, p);
                z = fmpz_fdiv_ui(den + n, p);

                if (y != n_mulmod2_preinv(z, n_mulmod2_preinv(x, n, p, pinv), p, pinv))
                {
                    flint_printf("FAIL\n");
                    flint_printf("n = %wu\n", n);
                    flint_printf("x = %wu\n", x);
                    flint_printf("y = %wu\n", y);
                    flint_printf("z = %wu\n", z);
                    flint_abort();
                }
            }
        }

        _fmpz_vec_clear(num, N);
        _fmpz_vec_clear(den, N);
    }

    {
        ulong p, k;

        /* exhaustive comparison over some small p and k */
        for (p = 5; p < 1000; p = n_nextprime(p, 1))
        {
            for (k = 2; k <= p - 3; k += 2)
                test_bern_modp_pow2(p, k);
        }

        /* a few larger values of p */
        for (p = n_nextprime(1000000, 1);
            p < 1000000 + 1000 * FLINT_MIN(10, 0.1 * flint_test_multiplier()); p = n_nextprime(p, 1))
        {
            k = 2 * (n_randlimb(state) % ((p-3)/2)) + 2;
            test_bern_modp_pow2(p, k);
        }

        /* these are slow */
        if (FLINT_BITS == 64 && 0.1 * flint_test_multiplier() >= 10)
        {
            test_bern_modp_pow2(2147483647, 10);
            test_bern_modp_pow2(2147483629, 10);
            test_bern_modp_pow2(2147483659, 10);
        }

        for (p = n_nextprime((1 << 15) - 1000, 1); p < (1 << 15); p = n_nextprime(p, 1))
        {
            for (iter = 0; iter < 10 * 0.1 * flint_test_multiplier(); iter++)
            {
                test_bern_modp_pow2(p, 2 * (n_randlimb(state) % ((p - 3)/2)) + 2);
            }
        }
    }

    TEST_FUNCTION_END(state);
}
