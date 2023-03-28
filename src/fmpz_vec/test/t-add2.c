/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"

#if defined(__AVX2__) || defined(__AVX512F__)
int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("add2....");
    fflush(stdout);

    /* Check that result is equal to _fmpz_vec_add */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz * ip1, * ip2, * res, * res2;
        slong len = n_randint(state, 100);

        ip1 = _fmpz_vec_init(len);
        ip2 = _fmpz_vec_init(len);
        res = _fmpz_vec_init(len);
        res2 = _fmpz_vec_init(len);
        _fmpz_vec_randtest(ip1, state, len, 200);
        _fmpz_vec_randtest(ip2, state, len, 200);

        _fmpz_vec_add(res, ip1, ip2, len);
        _fmpz_vec_set(res2, ip1, len);
        _fmpz_vec_add2(res2, ip2, len);

        result = _fmpz_vec_equal(res, res2, len);
        if (!result)
        {
            printf("i = %d\n", i);

            for (i = 0; i < len; i++)
                if (!fmpz_equal(res + i, res2 + i))
                    break;

            flint_printf("FAIL:\n");
            flint_printf("ip1[%d]  = ", i); fmpz_print(ip1 + i); flint_printf("\n");
            flint_printf("ip2[%d]  = ", i); fmpz_print(ip2 + i); flint_printf("\n");
            flint_printf("res[%d]  = ", i); fmpz_print(res + i); flint_printf("\n");
            flint_printf("res2[%d] = ", i); fmpz_print(res2 + i); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        _fmpz_vec_clear(ip1, len);
        _fmpz_vec_clear(ip2, len);
        _fmpz_vec_clear(res, len);
        _fmpz_vec_clear(res2, len);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
#else
int
main(void)
{
    flint_printf("add2....SKIPPED\n");
    return 0;
}
#endif
