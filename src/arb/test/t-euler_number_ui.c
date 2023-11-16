/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_vec.h"
#include "ulong_extras.h"
#include "arith.h"
#include "arb.h"

ulong euler_mod_p_powsum_noredc(ulong n, ulong p, const unsigned int * divtab);
ulong euler_mod_p_powsum(ulong n, ulong p, const unsigned int * divtab);
ulong euler_mod_p_powsum_1(ulong n, ulong p);

static void
divisor_table_odd(unsigned int * tab, slong len)
{
    slong i, j;

    tab[0] = 0;

    for (i = 1; i < len; i += 2)
    {
        tab[i] = 1;
        tab[i + 1] = i;
    }

    for (i = 3; i < len; i += 2)
    {
        for (j = 3; j <= i && i * j < len; j += 2)
        {
            tab[i * j]     = j;
            tab[i * j + 1] = i;
        }
    }
}

TEST_FUNCTION_START(arb_euler_number_ui, state)
{
    slong iter;

    {
        slong nmax;
        unsigned int * divtab_odd;
        fmpz * En;

        nmax = 1000 * FLINT_MIN(1.0, 0.1 * flint_test_multiplier());

        En = _fmpz_vec_init(nmax);
        arith_euler_number_vec(En, nmax);

        {
            fmpz_t E;
            slong n;
            double alpha;

            fmpz_init(E);

            for (n = 0; n < nmax; n++)
            {
                if (n_randint(state, 2))
                    alpha = -1.0;
                else
                    alpha = n_randint(state, 11) / (double) 10;

                arb_fmpz_euler_number_ui_multi_mod(E, n, alpha);

                if (!fmpz_equal(En + n, E))
                {
                    flint_printf("FAIL: n = %wd\n", n);
                    flint_printf("vec:    "); fmpz_print(En + n); flint_printf("\n");
                    flint_printf("single: "); fmpz_print(E); flint_printf("\n");
                    flint_abort();
                }

                arb_fmpz_euler_number_ui(E, n);

                if (!fmpz_equal(En + n, E))
                {
                    flint_printf("FAIL (2): n = %wd\n", n);
                    flint_printf("vec:    "); fmpz_print(En + n); flint_printf("\n");
                    flint_printf("single: "); fmpz_print(E); flint_printf("\n");
                    flint_abort();
                }
            }

            fmpz_clear(E);
        }

        /* test the mod p code */
        for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
        {
            ulong p, n, m, m1, m5;

            n = n_randint(state, nmax);
            p = 4 + n_randint(state, 10000);
            p = n_nextprime(p, 1);

            divtab_odd = flint_malloc(sizeof(unsigned int) * (p / 4 + 2));
            divisor_table_odd(divtab_odd, p / 4 + 1);

            m = fmpz_fdiv_ui(En + n, p);

            if (n_randint(state, 2))
                m5 = euler_mod_p_powsum(n, p, divtab_odd);
            else
                m5 = euler_mod_p_powsum_noredc(n, p, divtab_odd);

            if (n_randint(state, 30) == 0)
            {
                m1 = euler_mod_p_powsum_1(n, p);

                if (m1 != UWORD_MAX && m != m1)
                {
                    flint_printf("FAIL\n\n");
                    flint_printf("n = %wu, p = %wu, m = %wu, m1 = %wu\n\n", n, p, m, m1);
                    flint_abort();
                }
            }

            if (m5 != UWORD_MAX && m != m5)
            {
                flint_printf("FAIL\n\n");
                flint_printf("n = %wu, p = %wu, m = %wu, m5 = %wu\n\n", n, p, m, m5);
                flint_abort();
            }

            flint_free(divtab_odd);
        }

        _fmpz_vec_clear(En, nmax);
    }

    TEST_FUNCTION_END(state);
}
