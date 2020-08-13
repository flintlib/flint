/*
    Copyright (C) 2011 Andy Novocin
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Example program demonstrating the Zassenhaus factoring algorithm.
*/

#include <stdio.h>

#include "flint.h"
#include "fmpz_poly.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

int main(void)
{
    fmpz_poly_t f;
    fmpz_poly_factor_t facs;

    fmpz_poly_init(f);
    fmpz_poly_factor_init(facs);

    if (0)
    {
        FILE *polyfile;
        polyfile = fopen("examples/fmpz_poly_hensel_P1", "r");

        if (!polyfile)
        {
            flint_printf("Error.  Could not read P1 from file.\n");
            flint_abort();
        }
        fmpz_poly_fread(polyfile, f);
    }

    fmpz_poly_set_str(f, "63  1 1 1 -4 -7 -2 -6 -3 -7 18 7 25 -11 95 36 21 16 69 56 35 36 32 33 26 -26 -15 -14 -53 -96 67 72 -67 40 -79 -116 -452 -312 -260 -29 -1393 327 69 -28 -241 230 -54 -309 -125 -74 -450 -69 -3 66 -27 73 68 50 -63 -1290 372 31 -16 2");

    fmpz_poly_factor_zassenhaus(facs, f);

    flint_printf("Polynomial:\n");
    fmpz_poly_print(f);
    flint_printf("\nFactorisation:\n");
    fmpz_poly_factor_print(facs);

    fmpz_poly_clear(f);
    fmpz_poly_factor_clear(facs);

    return EXIT_SUCCESS;
}

