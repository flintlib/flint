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

    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

static void check(fmpz_t x, int expected)
{
    if (fmpz_fits_si(x) != expected)
    {
        printf("FAIL:\n\n");
        printf("x = "), fmpz_print(x), printf("\n");
        printf("fmpz_fits_si(x) = %d\n", fmpz_fits_si(x));
        printf("LONG_MIN = %ld\n", LONG_MIN);
        abort();
    }
}

int
main(void)
{
    long i;
    fmpz_t x;

    printf("fits_si....");
    fflush(stdout);

    fmpz_init(x);

    fmpz_set_si(x, COEFF_MIN);
    check(x, 1);

    fmpz_set_si(x, COEFF_MAX);
    check(x, 1);

    fmpz_set_si(x, LONG_MAX);
    check(x, 1);

    fmpz_set_si(x, LONG_MIN);
    check(x, 1);

    fmpz_set_ui(x, ULONG_MAX);
    check(x, 0);

    fmpz_set_ui(x, ULONG_MAX);
    fmpz_neg(x, x);
    check(x, 0);

    fmpz_set_si(x, LONG_MAX);
    fmpz_add_ui(x, x, 1);
    check(x, 0);

    fmpz_set_si(x, LONG_MIN);
    fmpz_sub_ui(x, x, 1);
    check(x, 0);

    for (i = 0; i < 1000; i++)
    {
        fmpz_set_ui(x, 1);
        fmpz_mul_2exp(x, x, i);
        check(x, i < FLINT_BITS - 1);
        fmpz_neg(x, x);
        check(x, i < FLINT_BITS);  /* LONG_MIN fits */
    }

    fmpz_clear(x);

    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
