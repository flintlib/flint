/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

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
    slong n;
    fmpz_mat_t S1, S2, P;

    if (argc != 2)
    {
        flint_printf("Syntax: stirling_matrix <integer>\n");
        return EXIT_FAILURE;
    }

    n = atoi(argv[1]);

    fmpz_mat_init(S1, n, n);
    fmpz_mat_init(S2, n, n);
    fmpz_mat_init(P, n, n);

    arith_stirling_matrix_1(S1);
    arith_stirling_matrix_2(S2);
    fmpz_mat_mul(P, S1, S2);

    flint_printf("S1 [Stirling numbers of 1st kind]:\n");
    fmpz_mat_print_pretty(S1);
    flint_printf("\n\n");

    flint_printf("S2 [Stirling numbers of 2nd kind]:\n");
    fmpz_mat_print_pretty(S2);
    flint_printf("\n\n");

    flint_printf("S1 * S2:\n");
    fmpz_mat_print_pretty(P);
    flint_printf("\n\n");

    fmpz_mat_clear(S1);
    fmpz_mat_clear(S2);
    fmpz_mat_clear(P);

    return EXIT_SUCCESS;
}
