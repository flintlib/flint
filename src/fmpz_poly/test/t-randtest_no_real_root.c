/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_poly.h"
TEST_FUNCTION_START(fmpz_poly_randtest_no_real_root, state)
{
    int iter;

    for (iter = 0; iter < 5000 * flint_test_multiplier(); iter++)
    {
        fmpz_poly_t p;
        slong len, i;
        flint_bitcnt_t bits;

        bits = 1 + n_randint(state, 30);
        len = 1 + n_randint(state, 20);

        fmpz_poly_init(p);
        fmpz_poly_randtest(p, state, 2*len, 2*bits);
        fmpz_poly_randtest_no_real_root(p, state, len, bits);

        if (fmpz_poly_length(p) > len)
        {
            printf("ERROR:\n");
            flint_printf("got length (= %wd) above the requested limit %wd\n",
                    fmpz_poly_length(p), len);
            printf("p = "); fmpz_poly_print(p); printf("\n");
            fflush(stdout);
            flint_abort();
        }

        /* check bit size */
        for (i = 0; i < fmpz_poly_length(p); i++)
        {
            if (fmpz_bits(p->coeffs + i) > bits)
            {
                printf("ERROR:\n");
                flint_printf("%wd-th coefficient exceed requested bit size\n", i);
                printf("p = "); fmpz_poly_print(p); printf("\n");
                flint_printf("bits = %wu\n", bits);
                fflush(stdout);
                flint_abort();
            }
        }

        /* check real roots */
        if (fmpz_poly_num_real_roots_sturm(p))
        {
            printf("ERROR:\n");
            flint_printf("polynomial has real root\n");
            printf("p = ");
            fmpz_poly_print(p);
            printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(p);
    }

    TEST_FUNCTION_END(state);
}
