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

TEST_FUNCTION_START(arb_digits_round_inplace, state)
{

    {
        char s[30];
        slong i, j, len, n;
        flint_bitcnt_t shift;
        fmpz_t inp, out, err, t;
        arf_rnd_t rnd;

        fmpz_init(inp);
        fmpz_init(out);
        fmpz_init(err);
        fmpz_init(t);

        for (i = 0; i < 100000 * 0.1 * flint_test_multiplier(); i++)
        {
            len = 1 + n_randint(state, 20);
            n = 1 + n_randint(state, 20);

            s[0] = (n_randint(state, 9) + '1');

            for (j = 1; j < len; j++)
                s[j] = (n_randint(state, 10) + '0');

            s[len] = '\0';

            fmpz_set_str(inp, s, 10);

            switch (n_randint(state, 3))
            {
                case 0:
                    rnd = ARF_RND_DOWN;
                    break;
                case 1:
                    rnd = ARF_RND_UP;
                    break;
                default:
                    rnd = ARF_RND_NEAR;
                    break;
            }

            _arb_digits_round_inplace(s, &shift, err, n, rnd);

            fmpz_set_str(out, s, 10);
            fmpz_set_ui(t, 10);
            fmpz_pow_ui(t, t, shift);
            fmpz_mul(t, t, out);
            fmpz_add(t, t, err);

            if (!fmpz_equal(t, inp) || (rnd == ARF_RND_UP && fmpz_sgn(err) > 0))
            {
                flint_printf("FAIL!\n");
                flint_printf("inp = "); fmpz_print(inp); flint_printf("\n\n");
                flint_printf("shift = %wd\n\n", shift);
                flint_printf("err = "); fmpz_print(err); flint_printf("\n\n");
                flint_printf("out = "); fmpz_print(out); flint_printf("\n\n");
                flint_printf(" t  = "); fmpz_print(t); flint_printf("\n\n");
                flint_abort();
            }
        }

        fmpz_clear(inp);
        fmpz_clear(out);
        fmpz_clear(err);
        fmpz_clear(t);
    }

    TEST_FUNCTION_END(state);
}
