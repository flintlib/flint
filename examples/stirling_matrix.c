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

/*
    Demo FLINT program for generating Stirling number matrices
    and inverting them.
*/

#include <stdlib.h>
#include <stdio.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "arith.h"

int main(int argc, char* argv[])
{
    len_t n;
    fmpz_mat_t S1, S2, P;

    if (argc != 2)
    {
        printf("Syntax: stirling_matrix <integer>\n");
        return EXIT_FAILURE;
    }

    n = atoi(argv[1]);

    fmpz_mat_init(S1, n, n);
    fmpz_mat_init(S2, n, n);
    fmpz_mat_init(P, n, n);

    arith_stirling_matrix_1(S1);
    arith_stirling_matrix_2(S2);
    fmpz_mat_mul(P, S1, S2);

    printf("S1 [Stirling numbers of 1st kind]:\n");
    fmpz_mat_print_pretty(S1);
    printf("\n\n");

    printf("S2 [Stirling numbers of 2nd kind]:\n");
    fmpz_mat_print_pretty(S2);
    printf("\n\n");

    printf("S1 * S2:\n");
    fmpz_mat_print_pretty(P);
    printf("\n\n");

    fmpz_mat_clear(S1);
    fmpz_mat_clear(S2);
    fmpz_mat_clear(P);

    return EXIT_SUCCESS;
}
