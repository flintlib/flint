/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acf.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("init_clear....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * arb_test_multiplier(); iter++)
    {
        acf_t x, y;

        acf_init(x);
        acf_init(y);

        arf_randtest(acf_realref(x), state, 200, 100);
        arf_randtest(acf_imagref(x), state, 200, 100);
        acf_set(y, x);

        acf_clear(x);
        acf_clear(y);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
