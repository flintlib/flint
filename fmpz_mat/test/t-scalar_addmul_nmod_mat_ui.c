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
#include "long_extras.h"

int
main(void)
{
    int i;
    flint_rand_t state;
    flint_randinit(state);

    printf("scalar_add/submul_nmod_mat_ui....");
    fflush(stdout);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mat_t A, B, C;
        nmod_mat_t M;
        len_t rows, cols;
        ulong c;
        ulong mod;

        rows = n_randint(state, 10);
        cols = n_randint(state, 10);
        mod = n_randtest_prime(state, 0);

        fmpz_mat_init(A, rows, cols);
        fmpz_mat_init(B, rows, cols);
        fmpz_mat_init(C, rows, cols);
        nmod_mat_init(M, rows, cols, mod);
        c = n_randtest(state);

        nmod_mat_randtest(M, state);
        fmpz_mat_set_nmod_mat_unsigned(A, M);

        fmpz_mat_randtest(B, state, 100);
        fmpz_mat_set(C, B);

        fmpz_mat_scalar_addmul_nmod_mat_ui(B, M, c);
        fmpz_mat_scalar_addmul_ui(C, A, c);

        if (!fmpz_mat_equal(B, C))
        {
            printf("FAIL!\n");
            abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        nmod_mat_clear(M);
    }

    flint_randclear(state);

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
