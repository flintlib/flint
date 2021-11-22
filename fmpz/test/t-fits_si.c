/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

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
        flint_printf("FAIL:\n\n");
        flint_printf("x = "), fmpz_print(x), flint_printf("\n");
        flint_printf("fmpz_fits_si(x) = %d\n", fmpz_fits_si(x));
        flint_printf("WORD_MIN = %wd\n", WORD_MIN);
        fflush(stdout);
        flint_abort();
    }
}

int
main(void)
{
    slong i;
    fmpz_t x;

    FLINT_TEST_INIT(state);
    
    flint_printf("fits_si....");
    fflush(stdout);

    fmpz_init(x);

    fmpz_set_si(x, COEFF_MIN);
    check(x, 1);

    fmpz_set_si(x, COEFF_MAX);
    check(x, 1);

    fmpz_set_si(x, WORD_MAX);
    check(x, 1);

    fmpz_set_si(x, WORD_MIN);
    check(x, 1);

    fmpz_set_ui(x, UWORD_MAX);
    check(x, 0);

    fmpz_set_ui(x, UWORD_MAX);
    fmpz_neg(x, x);
    check(x, 0);

    fmpz_set_si(x, WORD_MAX);
    fmpz_add_ui(x, x, 1);
    check(x, 0);

    fmpz_set_si(x, WORD_MIN);
    fmpz_sub_ui(x, x, 1);
    check(x, 0);

    for (i = 0; i < 1000; i++)
    {
        fmpz_set_ui(x, 1);
        fmpz_mul_2exp(x, x, i);
        check(x, i < FLINT_BITS - 1);
        fmpz_neg(x, x);
        check(x, i < FLINT_BITS);  /* WORD_MIN fits */
    }

    fmpz_clear(x);

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
