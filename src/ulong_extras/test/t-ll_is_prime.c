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

    int i, res;
    ulong hi, lo;

    fmpz_t f;
    fmpz_init(f);

    static const char * testcases[] = {
        /* Composites that pass many SPRP tests, from Sorenson-Webster */
        "-" "3404730287403079539471001",  /* oversize */
        "-" "3317044064679887385961981",  /* oversize */
        "0" "3110269097300703345712981",
        "0" "360681321802296925566181",
        "0" "164280218643672633986221",
        "0" "318665857834031151167461",
        "0" "7395010240794120709381",
        "0" "164280218643672633986221",
        "0" "318665857834031151167461",
        "0" "2995741773170734841812261",
        "0" "667636712015520329618581",
        "0" "3110269097300703345712981",
        "0" "552727880697763694556181",
        "0" "360681321802296925566181",
        "0" "7395010240794120709381",
        "0" "164280218643672633986221",
        /* largest certain prime for Sorenson-Webster */
        "1" "3317044064679887385961813",
        /* Some other easy composites */
        "0" "1180591620717411303428",
        "0" "1180591620717411303431",
        /* close to the boundary */
        "0" "340282366920938463463374607431768211451",
        "0" "340282366920938463463374607431768211441",
        "0" "340282366920938463463374607431768211453",
        "-" "340282366920938463463374607431768211297",
        "0" "340282366920938457671096968286969004029",
        "-" "340282366920938457671096968286969003909",
        "0" "340282366920938457744883944581807210495",
        "-" "340282366920938457744883944581807210471",
    };

    for (i = 0; i < sizeof(testcases) / sizeof(testcases[0]); i++)
    {
        fmpz_set_str(f, testcases[i] + 1, 10);
        fmpz_get_uiui(&hi, &lo, f);

        int expect = (testcases[i][0] == '0') ? 0 : ((testcases[i][0] == '1') ? 1 : -1);

        if ((res = n_ll_is_prime(hi, lo)) != expect)
        {
            flint_printf("FAIL\n");
            flint_printf("%s,  %d,  %d\n", testcases[i] + 1, res, expect);
            flint_abort();
        }
    }

    fmpz_clear(f);

#endif

    TEST_FUNCTION_END(state);
}
