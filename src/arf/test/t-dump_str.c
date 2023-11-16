/*
    Copyright (C) 2019 Julian Rüth

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include <string.h>
#include "arf.h"

TEST_FUNCTION_START(arf_dump_str, state)
{
    slong iter;

    /* just test no crashing... */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t x;
        char * s;

        arf_init(x);

        arf_randtest_special(x, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));

        s = arf_dump_str(x);

        flint_free(s);
        arf_clear(x);
    }

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t x, y;
        char * s;
        int conversion_error;

        arf_init(x);
        arf_init(y);

        arf_randtest_special(x, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        arf_randtest_special(y, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));

        s = arf_dump_str(x);
        conversion_error = arf_load_str(y, s);

        if (conversion_error || !arf_equal(x, y))
        {
            flint_printf("FAIL (roundtrip)  iter = %wd\n", iter);
            flint_printf("x = "); arf_printd(x, 50); flint_printf("\n\n");
            flint_printf("s = %s", s); flint_printf("\n\n");
            flint_printf("y = "); arf_printd(y, 50); flint_printf("\n\n");
            flint_abort();
        }

        flint_free(s);
        arf_clear(x);
        arf_clear(y);
    }

    TEST_FUNCTION_END(state);
}
