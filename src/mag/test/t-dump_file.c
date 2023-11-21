/*
    Copyright (C) 2019 Julian Rüth

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

TEST_FUNCTION_START(mag_dump_file, state)
{
    slong iter;

/* assume tmpfile() is broken on windows */
#if !defined(_MSC_VER) && !defined(__MINGW32__)
    /* just test no crashing... */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        mag_t x;
        FILE* tmp;

        mag_init(x);

        mag_randtest_special(x, state, 1 + n_randint(state, 100));

        tmp = tmpfile();
        mag_dump_file(tmp, x);
        fflush(tmp);
        rewind(tmp);
        mag_load_file(x, tmp);
        fclose(tmp);

        mag_clear(x);
    }

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        mag_t x, y, z;
        int conversion_error;
        FILE* tmp;

        mag_init(x);
        mag_init(y);
        mag_init(z);

        mag_randtest_special(x, state, 1 + n_randint(state, 100));
        mag_randtest_special(y, state, 1 + n_randint(state, 100));
        mag_randtest_special(z, state, 1 + n_randint(state, 100));

        tmp = tmpfile();
        mag_dump_file(tmp, x);
        fputc(' ', tmp);
        mag_dump_file(tmp, y);
        fflush(tmp);
        rewind(tmp);

        conversion_error = mag_load_file(z, tmp);
        if (conversion_error || !mag_equal(x, z))
        {
            flint_printf("FAIL (roundtrip)  iter = %wd\n", iter);
            flint_printf("x = "); mag_printd(x, 50); flint_printf("\n\n");
            flint_printf("z = "); mag_printd(z, 50); flint_printf("\n\n");
            flint_abort();
        }

        conversion_error = mag_load_file(z, tmp);
        if (conversion_error || !mag_equal(y, z))
        {
            flint_printf("FAIL (roundtrip)  iter = %wd\n", iter);
            flint_printf("y = "); mag_printd(y, 50); flint_printf("\n\n");
            flint_printf("z = "); mag_printd(z, 50); flint_printf("\n\n");
            flint_abort();
        }

        fclose(tmp);
        mag_clear(x);
        mag_clear(y);
        mag_clear(z);
    }
#endif

    TEST_FUNCTION_END(state);
}
