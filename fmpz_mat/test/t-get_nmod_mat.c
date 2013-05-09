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
#include "fmpz_mat.h"
#include "nmod_mat.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    flint_rand_t state;
    flint_randinit(state);

    printf("get/set_nmod_mat....");
    fflush(stdout);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mat_t A;
        nmod_mat_t M, M2;
        len_t rows, cols;
        mp_limb_t mod;

        rows = n_randint(state, 50);
        cols = n_randint(state, 50);

        mod = n_randtest_prime(state, 0);

        nmod_mat_init(M, rows, cols, mod);
        nmod_mat_init(M2, rows, cols, mod);
        fmpz_mat_init(A, rows, cols);

        nmod_mat_randtest(M, state);

        if (i % 2 == 0)
            fmpz_mat_set_nmod_mat(A, M);
        else
            fmpz_mat_set_nmod_mat_unsigned(A, M);

        fmpz_mat_scalar_mul_ui(A, A, 2UL);
        nmod_mat_add(M, M, M);
        fmpz_mat_get_nmod_mat(M2, A);

        if (!nmod_mat_equal(M, M2))
        {
            printf("FAIL!\n");
            abort();
        }

        fmpz_mat_clear(A);
        nmod_mat_clear(M);
        nmod_mat_clear(M2);
    }

    flint_randclear(state);

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
