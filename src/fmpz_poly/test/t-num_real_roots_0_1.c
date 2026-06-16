/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_poly.h"

TEST_FUNCTION_START(fmpz_poly_num_real_roots_0_1, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        fmpz_poly_t p;
        slong k1, k2;

        fmpz_poly_init(p);

        if (iter == 0)
        {
          // e-antic issue #301
          fmpz_poly_set_coeff_si(p, 0, 52);
          fmpz_poly_set_coeff_si(p, 3, -304);
          fmpz_poly_set_coeff_si(p, 8, -23);
        }
        else
        {
          do{
            fmpz_poly_randtest(p, state, 10, 10);
          } while (fmpz_poly_is_zero(p) || !fmpz_poly_is_squarefree(p));
        }

        k1 = fmpz_poly_num_real_roots_0_1_vca(p);
        k2 = fmpz_poly_num_real_roots_0_1_sturm(p);
        if (k1 != k2)
        {
            flint_printf("ERROR:\n");
            flint_printf("vca and Sturm disagree\n");
            flint_printf("p = "); fmpz_poly_print(p); flint_printf("\n");
            flint_printf("(vca) k1 = %wd  (Sturm) k2 = %wd\n", k1, k2);
            flint_abort();
        }

        fmpz_poly_clear(p);
    }

    TEST_FUNCTION_END(state);
}
