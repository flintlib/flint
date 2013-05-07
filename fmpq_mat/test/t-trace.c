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
#include "fmpq.h"
#include "fmpq_mat.h"
#include "ulong_extras.h"

int
main(void)
{
    long i;
    flint_rand_t state;

    printf("trace....");
    fflush(stdout);

    flint_randinit(state);

    /* Test trace(AB) = trace(BA) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mat_t A, B, AB, BA;
        fmpq_t trab, trba;
        long m, n;

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        fmpq_mat_init(A, m, n);
        fmpq_mat_init(B, n, m);
        fmpq_mat_init(AB, m, m);
        fmpq_mat_init(BA, n, n);

        fmpq_init(trab);
        fmpq_init(trba);

        fmpq_mat_randtest(A, state, 1 + n_randint(state, 100));
        fmpq_mat_randtest(B, state, 1 + n_randint(state, 100));

        fmpq_mat_mul(AB, A, B);
        fmpq_mat_mul(BA, B, A);

        fmpq_mat_trace(trab, AB);
        fmpq_mat_trace(trba, BA);

        if (!fmpq_equal(trab, trba))
        {
            printf("FAIL:\n");
            fmpq_mat_print(A), printf("\n");
            fmpq_mat_print(B), printf("\n");
            fmpq_mat_print(AB), printf("\n");
            fmpq_mat_print(BA), printf("\n");
            printf("tr(AB): "),  fmpq_print(trab),    printf("\n");
            printf("tr(BA): "),  fmpq_print(trba),    printf("\n");
            abort();
        }

        fmpq_mat_clear(A);
        fmpq_mat_clear(B);
        fmpq_mat_clear(AB);
        fmpq_mat_clear(BA);
        fmpq_clear(trab);
        fmpq_clear(trba);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
