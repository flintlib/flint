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

    printf("product_roots_fmpz_vec....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t P, Q, tmp;
        fmpz * x;
        len_t j, n, bits;

        n = n_randint(state, 100);
        bits = n_randint(state, 10);

        x = _fmpz_vec_init(n);
        _fmpz_vec_randtest(x, state, n, bits);

        fmpz_poly_init(P);
        fmpz_poly_init(Q);
        fmpz_poly_init(tmp);

        fmpz_poly_product_roots_fmpz_vec(P, x, n);

        fmpz_poly_set_ui(Q, 1UL);
        for (j = 0; j < n; j++)
        {
            fmpz_poly_zero(tmp);
            fmpz_poly_set_coeff_si(tmp, 1, -1L);
            fmpz_poly_set_coeff_fmpz(tmp, 0, x + j);
            fmpz_poly_neg(tmp, tmp);
            fmpz_poly_mul(Q, Q, tmp);
        }

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
        fmpz_poly_clear(tmp);
        _fmpz_vec_clear(x, n);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
