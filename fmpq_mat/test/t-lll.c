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

    Copyright (C) 2014 Abhinav Baid

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpq.h"
#include "fmpq_mat.h"
#include "fmpq_vec.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);


    flint_printf("lll....");
    fflush(stdout);

    /* check rref and |det| of LLL reduced matrix are equal to those of
       original matrix (randajtai used) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mat_t A, B, C;
        fmpq_mat_t Bq, Cq;
        fmpq_lll_t fl;
        fmpz_t temp1, temp2;

        slong m;

        m = n_randint(state, 10);

        fmpz_mat_init(A, m, m);
        fmpz_mat_init(B, m, m);
        fmpz_mat_init(C, m, m);
        fmpq_mat_init(Bq, m, m);
        fmpq_mat_init(Cq, m, m);
        fmpq_lll_context_init(fl, 3, 4, 1, 2);
        fmpz_init(temp1);
        fmpz_init(temp2);

        fmpz_mat_randajtai(A, state, 0.5);

        fmpz_mat_rref(B, temp1, A);
        fmpq_mat_set_fmpz_mat_div_fmpz(Bq, B, temp1);
        fmpz_abs(temp1, temp1);

        fmpq_mat_lll(A, A, fl);

        fmpz_mat_rref(C, temp2, A);
        fmpq_mat_set_fmpz_mat_div_fmpz(Cq, C, temp2);
        fmpz_abs(temp2, temp2);

        result = fmpq_mat_equal(Bq, Cq) && fmpz_equal(temp1, temp2);

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("B:\n");
            fmpz_mat_print_pretty(B);
            flint_printf("C:\n");
            fmpz_mat_print_pretty(C);
            fmpz_print(temp1);
            flint_printf("\n");
            fmpz_print(temp2);
            flint_printf("\n");
            abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        fmpq_mat_clear(Bq);
        fmpq_mat_clear(Cq);
        fmpq_lll_context_clear(fl);
        fmpz_clear(temp1);
        fmpz_clear(temp2);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
