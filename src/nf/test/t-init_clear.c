/*
    Copyright (C) 2013 William Hart
                  2019 Vincent Delecroix

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nf.h"

TEST_FUNCTION_START(nf_init_clear, state)
{
    int i;

    /* not necessarily monic */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        nf_t nf;
        fmpq_poly_t pol;

        fmpq_poly_init(pol);
        do {
           fmpq_poly_randtest_not_zero(pol, state,
                   2 + n_randint(state, 40),
                   10 + n_randint(state, 200));
        } while (fmpq_poly_degree(pol) < 1);

        nf_init(nf, pol);
        nf_clear(nf);

        nf_init(nf, pol);
        fmpq_poly_clear(pol);
        nf_clear(nf);
    }

    /* monic */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t pol;
        nf_t nf;

        fmpq_poly_init(pol);
        do {
           fmpq_poly_randtest_not_zero(pol, state,
                   2 + n_randint(state, 40),
                   10 + n_randint(state, 200));
        } while (fmpq_poly_degree(pol) < 1);
        fmpz_one(fmpq_poly_denref(pol));
        fmpz_one(pol->coeffs + fmpq_poly_degree(pol));

        nf_init(nf, pol);
        nf_clear(nf);

        nf_init(nf, pol);
        fmpq_poly_clear(pol);
        nf_clear(nf);
    }

    /* random */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        nf_t nf;

        nf_init_randtest(nf, state,
                2 + n_randint(state, 50),
                1 + n_randint(state, 200));
        nf_clear(nf);

        nf_init_randtest(nf, state,
                2 + n_randint(state, 50),
                1 + n_randint(state, 200));
        nf_clear(nf);
    }

    TEST_FUNCTION_END(state);
}
