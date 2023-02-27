/*
    Copyright (C) 2019 Julian Rüth

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "arb.h"

int main()
{
    flint_rand_t state;
    slong iter;

    flint_printf("dump_str/load_str....");
    fflush(stdout);
    flint_randinit(state);

    /* just test no crashing... */
    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        arb_t x;
        char * s;

        arb_init(x);

        arb_randtest_special(x, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));

        s = arb_dump_str(x);

        flint_free(s);
        arb_clear(x);
    }

    for (iter = 0; iter < 100000 * arb_test_multiplier(); iter++)
    {
        arb_t x, y;
        char * s;
        int conversion_error;

        arb_init(x);
        arb_init(y);

        arb_randtest_special(x, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        arb_randtest_special(y, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));

        s = arb_dump_str(x);
        conversion_error = arb_load_str(y, s);

        if (conversion_error || !arb_equal(x, y))
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

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
