/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_poly.h"
#include "nf_elem.h"

TEST_FUNCTION_START(nf_elem_is_gen, state)
{
    int i;

    /* is_gen(b) == equal(b, gen) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nf_t nf;
        nf_elem_t a, b;
        int r1, r2;

        nf_init_randtest(nf, state, 2 + n_randint(state, 5), 50);

        nf_elem_init(a, nf);
        nf_elem_init(b, nf);

        nf_elem_gen(a, nf);

        switch (n_randint(state, 3))
        {
            case 0:
                nf_elem_set(b, a, nf);
                break;
            case 1:
                nf_elem_add_si(b, a, 1 + n_randint(state, 3), nf);
                break;
            default:
                nf_elem_randtest(b, state, 10, nf);
        }

        r1 = nf_elem_is_gen(b, nf);
        r2 = nf_elem_equal(b, a, nf);

        if (r1 != r2)
        {
           flint_printf("FAIL:\n");
           flint_printf("K = "); nf_print(nf); flint_printf("\n");
           flint_printf("b = "); nf_elem_print_pretty(b, nf, "x"); flint_printf("\n");
           flint_printf("r1 = %d, r2 = %d\n", r1, r2);
           flint_abort();
        }

        nf_elem_clear(a, nf);
        nf_elem_clear(b, nf);

        nf_clear(nf);
    }

    /* Linear fields a1*x + a0 with non-unit leading coefficient */
    {
        slong a0, a1;

        for (a1 = -4; a1 <= 4; a1++)
        {
            for (a0 = -5; a0 <= 5; a0++)
            {
                fmpq_poly_t pol;
                nf_t nf;
                nf_elem_t a;

                if (a1 == 0)
                    continue;

                fmpq_poly_init(pol);
                fmpq_poly_set_coeff_si(pol, 0, a0);
                fmpq_poly_set_coeff_si(pol, 1, a1);
                nf_init(nf, pol);
                nf_elem_init(a, nf);

                nf_elem_gen(a, nf);
                if (!nf_elem_is_gen(a, nf))
                {
                    flint_printf("FAIL (linear field %wd*x + %wd)\n", a1, a0);
                    flint_abort();
                }

                nf_elem_add_si(a, a, 1, nf);
                if (nf_elem_is_gen(a, nf))
                {
                    flint_printf("FAIL (linear field %wd*x + %wd, gen + 1)\n", a1, a0);
                    flint_abort();
                }

                nf_elem_clear(a, nf);
                nf_clear(nf);
                fmpq_poly_clear(pol);
            }
        }
    }

    TEST_FUNCTION_END(state);
}
