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

TEST_FUNCTION_START(aprcl_config_gauss, state)
{
    int i;

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_t n, s2;
        aprcl_config conf;

        fmpz_init(s2);
        fmpz_init(n);
        fmpz_randtest_unsigned(n, state, 200);

        aprcl_config_gauss_init(conf, n);

        fmpz_mul(s2, conf->s, conf->s);
        if (fmpz_cmp(s2, n) <= 0)
        {
            flint_printf("FAIL:\n");
            flint_printf("s^2 = ");
            fmpz_print(s2);
            flint_printf(" <= ");
            fmpz_print(n);
            fflush(stdout);
            flint_abort();
        }

        aprcl_config_gauss_clear(conf);
        fmpz_clear(n);
        fmpz_clear(s2);
    }

    TEST_FUNCTION_END(state);
}
