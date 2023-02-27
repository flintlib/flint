/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013 Tom Bachmann (C++ adaptation)

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

#include <cstdio>
#include "fmpz_matxx.h"
#include "fmpzxx.h"
#include "arithxx.h"

using namespace std;
using namespace flint;

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        flint_printf("Syntax: stirling_matrix <integer>\n");
        return 1;
    }

    slong n = atoi(argv[1]);

    fmpz_matxx S1(stirling_matrix_1(n, n)), S2(stirling_matrix_2(n, n));

    flint_printf("S1 [Stirling numbers of 1st kind]:\n");
    print_pretty(S1);
    flint_printf("\n\n");

    flint_printf("S2 [Stirling numbers of 2nd kind]:\n");
    print_pretty(S2);
    flint_printf("\n\n");

    flint_printf("S1 * S2:\n");
    print_pretty(S1*S2);
    flint_printf("\n\n");

    return 0;
}
