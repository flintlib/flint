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

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

int main(void)
{
    slong iter;
    FLINT_TEST_INIT(state);

    flint_printf("chebyshev_u....");
    fflush(stdout);

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        ulong n;
        fmpz_t x, y, z;
        fmpz_poly_t p;

        fmpz_init(x);
        fmpz_init(y);
        fmpz_init(z);
        fmpz_poly_init(p);

        n = n_randint(state, 40);
        fmpz_randtest(x, state, 100);

        fmpz_chebyshev_u(y, n, x);

        fmpz_poly_chebyshev_u(p, n);
        fmpz_poly_evaluate_fmpz(z, p, x);

        if (!fmpz_equal(y, z))
        {
            flint_printf("FAIL: n = %wu\n\n", n);
            fmpz_print(x); printf("\n\n");
            fmpz_print(y); printf("\n\n");
            fmpz_print(z); printf("\n\n");
            abort();
        }

        fmpz_chebyshev_u(x, n, x);

        if (!fmpz_equal(x, z))
        {
            flint_printf("FAIL (aliasing): n = %wu\n\n", n);
            fmpz_print(x); printf("\n\n");
            fmpz_print(z); printf("\n\n");
            abort();
        }

        fmpz_clear(x);
        fmpz_clear(y);
        fmpz_clear(z);
        fmpz_poly_clear(p);
    }

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}

