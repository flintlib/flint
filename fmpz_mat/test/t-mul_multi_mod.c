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

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "ulong_extras.h"

int main(void)
{
    fmpz_mat_t A, B, C, D;
    long i;
    fmpz_randstate_t rnd;

    printf("mul_multi_mod....");
    fflush(stdout);

    fmpz_randinit(rnd);

    for (i = 0; i < 1000; i++)
    {
        long m, n, k;

        m = n_randint(50);
        n = n_randint(50);
        k = n_randint(50);

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(B, n, k);
        fmpz_mat_init(C, m, k);
        fmpz_mat_init(D, m, k);

        fmpz_mat_randtest(A, rnd, n_randint(200) + 1);
        fmpz_mat_randtest(B, rnd, n_randint(200) + 1);

        fmpz_mat_mul(C, A, B);
        fmpz_mat_mul_multi_mod(D, A, B);

        if (!fmpz_mat_equal(C, D))
        {
            printf("FAIL: results not equal\n");
            abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        fmpz_mat_clear(D);
    }

    fmpz_randclear(rnd);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
