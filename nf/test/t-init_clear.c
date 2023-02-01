/*=============================================================================

    This file is part of Antic.

    Antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version. See <http://www.gnu.org/licenses/>.

=============================================================================*/
/******************************************************************************

    Copyright (C) 2013 William Hart
                  2019 Vincent Delecroix

******************************************************************************/

#include <stdio.h>
#include "nf.h"

int
main(void)
{
    int i;
    flint_rand_t state;

    flint_printf("init/clear....");
    fflush(stdout);

    flint_randinit(state);

    /* not necessarily monic */
    for (i = 0; i < 500 * antic_test_multiplier(); i++)
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
    for (i = 0; i < 500 * antic_test_multiplier(); i++)
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
    for (i = 0; i < 500 * antic_test_multiplier(); i++)
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

    flint_randclear(state);
    flint_printf("PASS\n");
    return 0;
}
