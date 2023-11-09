/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"

#ifndef check
#define check check
/* Defined in t-abs_fits_ui.c, t-fits_si.c and t-moebius_mu.c */
static void check(fmpz_t input, int output, int expected)
{
    if (output != expected)
    {
        printf("FAIL:\n\n"
               "input = "), fmpz_print(input), printf("\n"
               "output = %d, expected = %d\n", output, expected);
        fflush(stdout);
        flint_abort();
    }
}
#endif

TEST_FUNCTION_START(fmpz_abs_fits_ui, state)
{
    slong i;
    fmpz_t x;

    fmpz_init(x);

    fmpz_set_si(x, COEFF_MIN);
    check(x, fmpz_abs_fits_ui(x), 1);

    fmpz_set_si(x, COEFF_MAX);
    check(x, fmpz_abs_fits_ui(x), 1);

    fmpz_set_ui(x, UWORD_MAX);
    check(x, fmpz_abs_fits_ui(x), 1);

    fmpz_set_ui(x, UWORD_MAX);
    fmpz_neg(x, x);
    check(x, fmpz_abs_fits_ui(x), 1);

    fmpz_set_ui(x, UWORD_MAX);
    fmpz_add_ui(x, x, UWORD(1));
    check(x, fmpz_abs_fits_ui(x), 0);

    fmpz_neg(x, x);
    check(x, fmpz_abs_fits_ui(x), 0);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_set_ui(x, UWORD(1));
        fmpz_mul_2exp(x, x, i);
        check(x, fmpz_abs_fits_ui(x), i < FLINT_BITS);
        fmpz_neg(x, x);
        check(x, fmpz_abs_fits_ui(x), i < FLINT_BITS);
    }

    fmpz_clear(x);

    TEST_FUNCTION_END(state);
}
