/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "aprcl.h"

TEST_FUNCTION_START(aprcl_is_prime_gauss, state)
{
    int i;

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        int pbprime, cycloprime;
        fmpz_t n;
        fmpz_init(n);

        fmpz_randtest_unsigned(n, state, 50);
        while (fmpz_cmp_ui(n, 100) <= 0)
            fmpz_randtest_unsigned(n, state, 50);

        pbprime = fmpz_is_probabprime(n);
        cycloprime = aprcl_is_prime_gauss(n);

        if (pbprime != cycloprime)
        {
            flint_printf("FAIL\n");
            flint_printf("Testing number = ");
            fmpz_print(n);
            flint_printf("\nis_probabprime = %i, aprcl_is_prime_gauss = %i\n", pbprime, cycloprime);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(n);
    }

    {
        int result;
        fmpz_t n;
        fmpz_init(n);
        result = 1;

        /* Test big primes. */
        fmpz_set_str(n, "40206835204840513073", 10);
        if (aprcl_is_prime_gauss(n) == 0)
            result = 0;

        /* 521419622856657689423872613771 % 4 == 3 */
        fmpz_set_str(n, "521419622856657689423872613771", 10);
        if (aprcl_is_prime_gauss(n) == 0)
            result = 0;

        /* Very slow. */
        if (flint_test_multiplier() > 10)
        {
            /* 5991810554633396517767024967580894321153 % 4 == 1 */
            fmpz_set_str(n, "5991810554633396517767024967580894321153", 10);
            if (aprcl_is_prime_gauss(n) == 0)
                result = 0;
        }

        /* Test big composite. */
        /* 1500450271 * 5915587277 */
        fmpz_set_str(n, "8876044532898802067", 10);
        if (aprcl_is_prime_gauss(n) == 1)
            result = 0;

        /* 5915587277 * 54673257461630679457 */
        fmpz_set_str(n, "323424426232167763068694468589", 10);
        if (aprcl_is_prime_gauss(n) == 1)
            result = 0;

        /* 48112959837082048697 * 66405897020462343733 */
        fmpz_set_str(n, "3194984256290911228520362769161858765901", 10);
        if (aprcl_is_prime_gauss(n) == 1)
            result = 0;

        if (result == 0)
        {
            flint_printf("FAIL\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(n);
    }

    TEST_FUNCTION_END(state);
}
