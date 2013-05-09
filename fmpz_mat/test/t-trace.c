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
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "ulong_extras.h"

int
main(void)
{
    len_t i;
    flint_rand_t state;

    printf("trace....");
    fflush(stdout);

    flint_randinit(state);

    /* Test trace(AB) = trace(BA) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mat_t A, B, AB, BA;
        fmpz_t trab, trba;
        len_t m, n;

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(B, n, m);
        fmpz_mat_init(AB, m, m);
        fmpz_mat_init(BA, n, n);

        fmpz_init(trab);
        fmpz_init(trba);

        fmpz_mat_randtest(A, state, 1 + n_randint(state, 100));
        fmpz_mat_randtest(B, state, 1 + n_randint(state, 100));

        fmpz_mat_mul(AB, A, B);
        fmpz_mat_mul(BA, B, A);

        fmpz_mat_trace(trab, AB);
        fmpz_mat_trace(trba, BA);

        if (!fmpz_equal(trab, trba))
        {
            printf("FAIL:\n");
            fmpz_mat_print_pretty(A), printf("\n");
            fmpz_mat_print_pretty(B), printf("\n");
            fmpz_mat_print_pretty(AB), printf("\n");
            fmpz_mat_print_pretty(BA), printf("\n");
            printf("tr(AB): "),  fmpz_print(trab),    printf("\n");
            printf("tr(BA): "),  fmpz_print(trba),    printf("\n");
            abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(AB);
        fmpz_mat_clear(BA);
        fmpz_clear(trab);
        fmpz_clear(trba);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
