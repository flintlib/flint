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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    flint_rand_t state;

    printf("newton_to_monomial....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t f, g;
        fmpz * r;
        long k, n;

        fmpz_poly_init(f);
        fmpz_poly_init(g);

        fmpz_poly_randtest(f, state, 1 + n_randint(state, 20),
                                     1 + n_randint(state, 200));

        n = fmpz_poly_length(f);
        r = _fmpz_vec_init(n);

        for (k = 0; k < n; k++)
            fmpz_randtest(r + k, state, n_randint(state, 200));

        fmpz_poly_set(g, f);

        _fmpz_poly_newton_to_monomial(g->coeffs, r, n);
        _fmpz_poly_monomial_to_newton(g->coeffs, r, n);

        if (!fmpz_poly_equal(f, g))
        {
            printf("FAIL: roundtrip\n");
            fmpz_poly_print(f); printf("\n");
            fmpz_poly_print(g); printf("\n");
            abort();
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        _fmpz_vec_clear(r, n);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
