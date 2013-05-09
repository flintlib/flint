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

    Copyright (C) 2007 David Harvey, William Hart
    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

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
    len_t N = 0;

    if (argc == 2)
        N = atol(argv[1]);

    if (argc != 2 || N < 1)
    {
        printf("Syntax: delta_qexp <integer>\n");
        printf("where <integer> is the (positive) number of terms to compute\n");
        return EXIT_FAILURE;
    }

    fmpz_init(c);
    fmpz_init(n);

    fmpz_set_si(n, N);
    arith_ramanujan_tau(c, n);

    printf("Coefficient of q^%ld is ", N);
    fmpz_print(c);
    printf("\n");

    fmpz_clear(c);
    fmpz_clear(n);

    return EXIT_SUCCESS;
}

