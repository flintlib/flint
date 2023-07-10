/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int main(void)
{
    slong iter;
    flint_rand_t state;

    flint_printf("sp2gz_inv....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 10);
        slong bits = n_randint(state, 10);
        fmpz_mat_t m1, m2;
        fmpz_t den;

        fmpz_mat_init(m1, 2 * g, 2 * g);
        fmpz_mat_init(m2, 2 * g, 2 * g);
        fmpz_init(den);

        sp2gz_randtest(m1, state, bits);
        sp2gz_inv(m2, m1);
        fmpz_mat_inv(m1, den, m1);

        if (!fmpz_mat_equal(m1, m2) || !fmpz_is_one(den))
        {
            flint_printf("FAIL\n\n");
            flint_abort();            
        }

        fmpz_mat_clear(m1);
        fmpz_mat_clear(m2);
        fmpz_clear(den);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
