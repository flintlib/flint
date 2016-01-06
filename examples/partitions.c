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
