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
    int i, result;
    flint_rand_t state;

    printf("interpolate_fmpz_vec....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t P, Q;
        fmpz *x, *y;
        long j, n, npoints, bits;

        npoints = n_randint(state, 50);
        n = n_randint(state, npoints + 1);
        bits = n_randint(state, 100);

        x = _fmpz_vec_init(npoints);
        y = _fmpz_vec_init(npoints);

        fmpz_poly_init(P);
        fmpz_poly_init(Q);

        fmpz_poly_randtest(P, state, n, bits);

        for (j = 0; j < npoints; j++)
            fmpz_set_si(x + j, -npoints/2 + j);

        fmpz_poly_evaluate_fmpz_vec(y, P, x, npoints);
        fmpz_poly_interpolate_fmpz_vec(Q, x, y, npoints);

        result = (fmpz_poly_equal(P, Q));
        if (!result)
        {
            printf("FAIL (P != Q):\n");
            fmpz_poly_print(P), printf("\n\n");
            fmpz_poly_print(Q), printf("\n\n");
            abort();
        }

        fmpz_poly_clear(P);
        fmpz_poly_clear(Q);
        _fmpz_vec_clear(x, npoints);
        _fmpz_vec_clear(y, npoints);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
