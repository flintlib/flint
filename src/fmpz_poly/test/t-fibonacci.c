/*
    Copyright (C) 2016 Shivin Srivastava

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

int main(void)
{
    fmpz_poly_t Pn, Pn1, Pn2, R;

    slong n;

    FLINT_TEST_INIT(state);

    flint_printf("fibonacci....");
    fflush(stdout);

    fmpz_poly_init(Pn);
    fmpz_poly_init(Pn1);
    fmpz_poly_init(Pn2);
    fmpz_poly_init(R);

    fmpz_poly_set_ui(Pn, UWORD(0));
    fmpz_poly_set_ui(Pn1, UWORD(1));

    for (n = 0; n <= 500; n++)
    {
        fmpz_poly_fibonacci(R, n);

        if (!fmpz_poly_equal(Pn, R))
        {
            flint_printf("FAIL: n = %wd\n", n);
            flint_printf("Direct: "); fmpz_poly_print_pretty(R, "x"); flint_printf("\n");
            flint_printf("Recur.: "); fmpz_poly_print_pretty(Pn, "x"); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_shift_left(Pn2, Pn1, 1);
        fmpz_poly_add(Pn2, Pn2, Pn);

        fmpz_poly_swap(Pn, Pn1);
        fmpz_poly_swap(Pn1, Pn2);
    }

    fmpz_poly_clear(Pn);
    fmpz_poly_clear(Pn1);
    fmpz_poly_clear(Pn2);
    fmpz_poly_clear(R);

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
