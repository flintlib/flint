/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmpcompat.h"
#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(n_ll_is_prime, state)
{
#if FLINT_BITS == 64

    int i, res, expect;
    ulong hi, lo;

    fmpz_t f;
    fmpz_init(f);

    /* NB: this function is also tested indirectly via fmpz_is_prime */

    static const char * testcases[18] = {
        /* Composites that pass many SPRP tests, from Sorenson-Webster */
        "3404730287403079539471001",  /* oversize */
        "3317044064679887385961981",  /* oversize */
        "3110269097300703345712981",
        "360681321802296925566181",
        "164280218643672633986221",
        "318665857834031151167461",
        "7395010240794120709381",
        "164280218643672633986221",
        "318665857834031151167461",
        "2995741773170734841812261",
        "667636712015520329618581",
        "3110269097300703345712981",
        "552727880697763694556181",
        "360681321802296925566181",
        "7395010240794120709381",
        "164280218643672633986221",
        /* Some other easy composites */
        "1180591620717411303428",
        "1180591620717411303431",
    };

    for (i = 0; i < 18; i++)
    {
        if (i <= 1)
            expect = -1;
        else
            expect = 0;

        fmpz_set_str(f, testcases[i], 10);

        fmpz_get_uiui(&hi, &lo, f);

        if ((res = n_ll_is_prime(hi, lo)) != expect)
        {
            flint_printf("FAIL\n");
            flint_printf("%s,  %d,  %d\n", testcases[i], res, expect);
        }
    }

    fmpz_clear(f);

#endif

    TEST_FUNCTION_END(state);
}
