/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    FLINT program for demonstrating the Integer Partition function.
*/

#include <stdio.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "arith.h"

int
main(int argc, char * argv[])
{
    fmpz_t x;
    ulong n;

    if (argc != 2)
    {
        flint_printf("usage: partitions n\n");
        return 1;
    }

    flint_sscanf(argv[1], "%wu", &n);

    flint_printf("p(%wu) = \n", n);

    fmpz_init(x);
    arith_number_of_partitions(x, n);
    fmpz_print(x);
    flint_printf("\n");
    fmpz_clear(x);

    return 0;
}
