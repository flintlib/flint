/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_extras.h"
#include "fmpzi.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("remove_one_plus_i....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * arb_test_multiplier(); iter++)
    {
        fmpzi_t x, y, z, s;
        slong v, w;

        fmpzi_init(x);
        fmpzi_init(y);
        fmpzi_init(z);
        fmpzi_init(s);

        do {
            fmpzi_randtest(x, state, 100);
        } while (fmpzi_is_zero(x));

        fmpzi_remove_one_plus_i(x, x);

        v = n_randint(state, 10);
        fmpz_one(fmpzi_realref(y));
        fmpz_one(fmpzi_imagref(y));
        fmpzi_pow_ui(y, y, v);

        fmpzi_mul(z, x, y);

        fmpzi_randtest(s, state, 200);
        w = fmpzi_remove_one_plus_i(s, z);

        if (v != w || !fmpzi_equal(s, x))
        {
            flint_printf("FAIL\n");
            flint_printf("x = "); fmpzi_print(x); printf("\n");
            flint_printf("v = %wd\n", v);
            flint_printf("y = "); fmpzi_print(y); printf("\n");
            flint_printf("z = "); fmpzi_print(z); printf("\n");
            flint_printf("s = "); fmpzi_print(s); printf("\n");
            flint_printf("w = %wd\n", w);
            flint_abort();
        }

        fmpzi_clear(x);
        fmpzi_clear(y);
        fmpzi_clear(z);
        fmpzi_clear(s);
    }

    flint_randclear(state);
    flint_cleanup_master();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
