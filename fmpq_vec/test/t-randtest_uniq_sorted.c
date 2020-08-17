/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpq_vec.h"

int main()
{
    int iter;
    FLINT_TEST_INIT(state);

    printf("randtest_uniq_sorted....");

    for (iter = 0; iter < 500 * flint_test_multiplier(); iter++)
    {
        slong i;
        slong n;
        flint_bitcnt_t bits;
        fmpq * vec;

        n = n_randint(state, 20);
        vec = _fmpq_vec_init(n);
        bits = 21 + n_randint(state, 200);
        _fmpq_vec_randtest_uniq_sorted(vec, state, n, bits);

        /* check bit size */
        for (i = 0; i < n; i++)
        {
            if ((fmpz_bits(fmpq_numref(vec + i)) > bits) ||
                (fmpz_bits(fmpq_denref(vec + i)) > bits))
            {
                printf("ERROR\n");
                flint_printf("num size: %wu\n",
                        fmpz_bits(fmpq_numref(vec + i)));
                flint_printf("den size: %wu\n",
                        fmpz_bits(fmpq_denref(vec + i)));
                flint_printf("bits    : %wu\n", bits);
                abort();
            }
        }

        /* check uniqueness and ordering */
        for (i = 0; i < n-1; i++)
        {
            if (fmpq_cmp(vec + i + 1, vec + i) <= 0)
            {
                printf("ERROR:\n");
                flint_printf("n = %wu\n", n);
                flint_printf("bits = %wu\n", bits);
                flint_printf("got vec[%wd] = ", i); fmpq_print(vec + i);
                flint_printf(" and vec[%wd] = ", i+1); fmpq_print(vec + i + 1);
                printf("\n");
                abort();
            }
        }

        _fmpq_vec_clear(vec, n);
    }

    FLINT_TEST_CLEANUP(state);

    printf("PASS\n");
    return 0;
}
