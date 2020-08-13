/*
    Copyright (C) 2007 David Harvey, William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Demo FLINT program for computing the q-expansion of the delta function.
*/

#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>

#include "flint.h"
#include "fmpz.h"
#include "arith.h"

int main(int argc, char* argv[])
{
    fmpz_t c, n;
    slong N = 0;

    if (argc == 2)
        N = atol(argv[1]);

    if (argc != 2 || N < 1)
    {
        flint_printf("Syntax: delta_qexp <integer>\n");
        flint_printf("where <integer> is the (positive) number of terms to compute\n");
        return EXIT_FAILURE;
    }

    fmpz_init(c);
    fmpz_init(n);

    fmpz_set_si(n, N);
    arith_ramanujan_tau(c, n);

    flint_printf("Coefficient of q^%wd is ", N);
    fmpz_print(c);
    flint_printf("\n");

    fmpz_clear(c);
    fmpz_clear(n);

    return EXIT_SUCCESS;
}

