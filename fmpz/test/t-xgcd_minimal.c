/*
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

int
main(void)
{
    int i, result;
    fmpz_t maxval;
    FLINT_TEST_INIT(state);
    fmpz_set_str(maxval, "3", 10);

    flint_printf("xgcd_minimal....");
    fflush(stdout);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t d1, d2, a1, a2, b1, b2, f, g, tmpa, tmpb;

        fmpz_init(d1);
        fmpz_init(a1);
        fmpz_init(b1);
        fmpz_init(d2);
        fmpz_init(a2);
        fmpz_init(b2);
        fmpz_init(f);
        fmpz_init(g);
        fmpz_init(tmpa);
        fmpz_init(tmpb);

        fmpz_randm(f, state, maxval);
        fmpz_randm(g, state, maxval);
        if (fmpz_is_zero(g))
            fmpz_one(g);
        if (n_randint(state, 2)) fmpz_neg(g, g);
        if (n_randint(state, 2)) fmpz_neg(f, f);

        fmpz_xgcd_minimal(d1, a1, b1, f, g);
        fmpz_xgcd_minimal(d2, b2, a2, g, f);

        fmpz_divexact(tmpa, g, d1);
        fmpz_divexact(tmpb, f, d1);
        result = (fmpz_equal(d1, d2)
               /* && fmpz_equal(a1, a2) */ /* Julia's gcdx does not ensure */
               /* && fmpz_equal(b1, b2) */ /* symmetry.                    */
               && (fmpz_is_zero(g) || fmpz_is_zero(f)
               || (fmpz_cmpabs(a1, tmpa) <= 0 && fmpz_cmpabs(b1, tmpb) <= 0)));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("tmpa = "), fmpz_print(tmpa), flint_printf("\n");
            flint_printf("fmpz_cmpabs(a1, tmpa) = "), flint_printf("%d", fmpz_cmpabs(a1, tmpa)), flint_printf("\n");
            flint_printf("tmpb = "), fmpz_print(tmpb), flint_printf("\n");
            flint_printf("fmpz_cmpabs(b1, tmpb) = "), flint_printf("%d", fmpz_cmpabs(b1, tmpb)), flint_printf("\n");
            flint_printf("d1 = "), fmpz_print(d1), flint_printf("\n");
            flint_printf("a1 = "), fmpz_print(a1), flint_printf("\n");
            flint_printf("b1 = "), fmpz_print(b1), flint_printf("\n");
            flint_printf("d2 = "), fmpz_print(d2), flint_printf("\n");
            flint_printf("a2 = "), fmpz_print(a2), flint_printf("\n");
            flint_printf("b2 = "), fmpz_print(b2), flint_printf("\n");
            flint_printf("f = "), fmpz_print(f), flint_printf("\n");
            flint_printf("g = "), fmpz_print(g), flint_printf("\n");
            abort();
        }

        fmpz_clear(d1);
        fmpz_clear(a1);
        fmpz_clear(b1);
        fmpz_clear(d2);
        fmpz_clear(a2);
        fmpz_clear(b2);
        fmpz_clear(f);
        fmpz_clear(g);
        fmpz_init(tmpa);
        fmpz_init(tmpb);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
