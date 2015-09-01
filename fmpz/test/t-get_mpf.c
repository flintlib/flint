/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("get/set_mpf....");
    fflush(stdout);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t x, y;
        mpz_t z;
        mpf_t a, b, tmp;

        fmpz_init(x);
        fmpz_init(y);
        mpz_init(z);
        mpf_inits(a, b, tmp, NULL);

        fmpz_randtest(x, state, 200);
        fmpz_get_mpz(z, x);

        fmpz_get_mpf(a, x);
        mpf_set_z(b, z);

        result = (mpf_cmp(a, b) == 0);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("x = "), fmpz_print(x), flint_printf("\n");
            flint_printf("a = "), mpf_out_str(stdout, 10, 0, a),
                flint_printf("\n");
            flint_printf("b = "), mpf_out_str(stdout, 10, 0, b),
                flint_printf("\n");
            abort();
        }

        mpf_set_d(tmp,
                  (n_randtest(state) / (double) n_randtest_not_zero(state)));
        mpf_mul(a, a, tmp);

        fmpz_set_mpf(x, a);
        mpz_set_f(z, a);

        fmpz_set_mpz(y, z);
        result = fmpz_equal(x, y);

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("x = "), fmpz_print(x), flint_printf("\n");
            flint_printf("y = "), fmpz_print(y), flint_printf("\n");
            flint_printf("a = "), mpf_out_str(stdout, 10, 0, a),
                flint_printf("\n");
            abort();
        }

        fmpz_clear(x);
        fmpz_clear(y);
        mpz_clear(z);
        mpf_clears(a, b, tmp, NULL);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
