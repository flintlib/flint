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

const long test_A[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
const long test_B[] = { 10, 11, 12, 13, 14, 15, 16, 17, 18 };
const long test_C[] = { 84, 90, 96, 201, 216, 231, 318, 342, 366 };

int main(void)
{
    fmpz_mat_t A, B, C;
    long i;

    printf("mul....");
    fflush(stdout);

    fmpz_mat_init(A, 3, 3);
    fmpz_mat_init(B, 3, 3);
    fmpz_mat_init(C, 3, 3);

    for (i = 0; i < 9; i++)
    {
        fmpz_set_ui(A->entries+i, test_A[i]);
        fmpz_set_ui(B->entries+i, test_B[i]);
    }

    fmpz_mat_mul(C, A, B);

    for (i = 0; i < 9; i++)
    {
        if (fmpz_get_ui(C->entries+i) != test_C[i])
        {
            _fmpz_vec_print(C->entries, 9);
            abort();
        }
    }

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
