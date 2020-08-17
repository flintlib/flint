/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "fmpz_poly.h"

int main()
{
    int iter;

    FLINT_TEST_INIT(state);

    printf("num_real_roots....");

    /* call with random nonzero polynomials */
    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        slong k;
        fmpz_poly_t p;

        fmpz_poly_init(p);
        fmpz_poly_randtest_not_zero(p, state, 20, 10 + n_randint(state, 100));
        k = fmpz_poly_num_real_roots(p);
        if (k < 0 || k > fmpz_poly_degree(p))
        {
            printf("ERROR:\n");
            flint_printf("got k in wrong range k = %wd\n", k);
            printf("p = "); fmpz_poly_print(p); printf("\n");
            abort();
        }

        fmpz_poly_clear(p);
    }

    for (iter = 0; iter < 5000 * flint_test_multiplier(); iter++)
    {
        slong k1, k2;
        fmpz_poly_t p;

        fmpz_poly_init(p);
        /* currently the code of num_real_roots only has a special branch */
        /* for length <= 5. We only test these cases.                     */
        do
        {
            fmpz_poly_randtest_not_zero(p, state, 1 + n_randint(state, 5), 100);
        } while (!fmpz_poly_is_squarefree(p));

        k1 = fmpz_poly_num_real_roots_sturm(p);
        k2 = fmpz_poly_num_real_roots(p);
        if (k1 != k2)
        {
            printf("ERROR:\n");
            flint_printf("found k1=%wd and k2=%wd\n", k1, k2);
            printf("p = "); fmpz_poly_print(p); printf("\n");
            abort();
        }

        fmpz_poly_clear(p);
    }

    FLINT_TEST_CLEANUP(state);

    printf("PASS\n");
    return 0;
}
