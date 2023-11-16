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
#include "arb.h"

TEST_FUNCTION_START(arb_dump_file, state)
{
    slong iter;

/* assume tmpfile() is broken on windows */
#if !defined(_MSC_VER) && !defined(__MINGW32__)

    /* just test no crashing... */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x;
        FILE* tmp;

        arb_init(x);

        arb_randtest_special(x, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));

        tmp = tmpfile();
        if (tmp == NULL) {
            flint_printf("FAIL (creating temporary file)  iter = %wd\n\n", iter);
            flint_abort();
        }
        arb_dump_file(tmp, x);
        fflush(tmp);
        rewind(tmp);
        arb_load_file(x, tmp);
        fclose(tmp);

        arb_clear(x);
    }

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x, y, z;
        int conversion_error;
        FILE* tmp;

        arb_init(x);
        arb_init(y);
        arb_init(z);

        arb_randtest_special(x, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        arb_randtest_special(y, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        arb_randtest_special(z, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));

        tmp = tmpfile();
        arb_dump_file(tmp, x);
        fputc(' ', tmp);
        arb_dump_file(tmp, y);
        fflush(tmp);
        rewind(tmp);

        conversion_error = arb_load_file(z, tmp);
        if (conversion_error || !arb_equal(x, z))
        {
            flint_printf("FAIL (roundtrip)  iter = %wd\n", iter);
            flint_printf("x = "); arb_printd(x, 50); flint_printf("\n\n");
            flint_printf("z = "); arb_printd(z, 50); flint_printf("\n\n");
            flint_abort();
        }

        conversion_error = arb_load_file(z, tmp);
        if (conversion_error || !arb_equal(y, z))
        {
            flint_printf("FAIL (roundtrip)  iter = %wd\n", iter);
            flint_printf("y = "); arb_printd(y, 50); flint_printf("\n\n");
            flint_printf("z = "); arb_printd(z, 50); flint_printf("\n\n");
            flint_abort();
        }

        fclose(tmp);
        arb_clear(x);
        arb_clear(y);
        arb_clear(z);
    }

#endif

    TEST_FUNCTION_END(state);
}
