/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include <string.h>
#include "arb.h"

TEST_FUNCTION_START(arb_get_str, state)
{
    slong iter;

    /* just test no crashing... */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x;
        char * s;
        slong n;

        arb_init(x);

        arb_randtest_special(x, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));

        n = 1 + n_randint(state, 300);

        s = arb_get_str(x, n, (n_randint(state, 2) * ARB_STR_MORE)
                            |  (n_randint(state, 2) * ARB_STR_NO_RADIUS)
                            | (ARB_STR_CONDENSE * n_randint(state, 50)));

        flint_free(s);
        arb_clear(x);
    }

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x, y;
        char * s;
        slong n, prec;
        int conversion_error;

        arb_init(x);
        arb_init(y);

        arb_randtest_special(x, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        arb_randtest_special(y, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));

        n = 1 + n_randint(state, 300);
        prec = 2 + n_randint(state, 1000);

        s = arb_get_str(x, n, n_randint(state, 2) * ARB_STR_MORE);
        conversion_error = arb_set_str(y, s, prec);

        if (conversion_error || !arb_contains(y, x))
        {
            flint_printf("FAIL (roundtrip)  iter = %wd\n", iter);
            flint_printf("x = "); arb_printd(x, 50); flint_printf("\n\n");
            flint_printf("s = %s", s); flint_printf("\n\n");
            flint_printf("y = "); arb_printd(y, 50); flint_printf("\n\n");
            flint_abort();
        }

        flint_free(s);
        arb_clear(x);
        arb_clear(y);
    }

    /* test ARB_STR_NO_RADIUS */
    {
        arb_t x;
        char * s;

        arb_init(x);

        arb_set_str(x, "3.1415926535897932", 53);
        s = arb_get_str(x, 10, ARB_STR_NO_RADIUS);
        if (strcmp(s, "3.141592654"))
        {
            flint_printf("FAIL (ARB_STR_NO_RADIUS)\n");
            flint_printf("%s\n", s);
            flint_abort();
        }
        flint_free(s);

        arb_set_str(x, "+/- 3.45e-10", 53);
        s = arb_get_str(x, 10, ARB_STR_NO_RADIUS);
        if (strcmp(s, "0e-9"))
        {
            flint_printf("FAIL (ARB_STR_NO_RADIUS)\n");
            flint_printf("%s\n", s);
            flint_abort();
        }
        flint_free(s);

        arb_set_str(x, "+/- 3.45e10", 53);
        s = arb_get_str(x, 10, ARB_STR_NO_RADIUS);
        if (strcmp(s, "0e+11"))
        {
            flint_printf("FAIL (ARB_STR_NO_RADIUS)\n");
            flint_printf("%s\n", s);
            flint_abort();
        }
        flint_free(s);

        arb_set_str(x, "5e10 +/- 6e10", 53);
        s = arb_get_str(x, 10, ARB_STR_NO_RADIUS);
        if (strcmp(s, "0e+12"))
        {
            flint_printf("FAIL (ARB_STR_NO_RADIUS)\n");
            flint_printf("%s\n", s);
            flint_abort();
        }
        flint_free(s);

        arb_set_str(x, "5e-100000000000000000002 +/- 5e-100000000000000000002", 53);
        s = arb_get_str(x, 10, ARB_STR_NO_RADIUS);
        if (strcmp(s, "0e-100000000000000000000"))
        {
            flint_printf("FAIL (ARB_STR_NO_RADIUS)\n");
            flint_printf("%s\n", s);
            flint_abort();
        }
        flint_free(s);

        arb_clear(x);
    }

    TEST_FUNCTION_END(state);
}
