/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_poly.h"

void cyclotomic_naive(fmpz_poly_t poly, ulong n)
{
    fmpz_poly_t t;
    slong d;

    fmpz_poly_init(t);

    fmpz_poly_set_ui(poly, UWORD(1));
    for (d = 1; d <= n; d++)
    {
        if (n % d == 0)
        {
            if (n_moebius_mu(n / d) == 1)
            {
                fmpz_poly_zero(t);
                fmpz_poly_set_coeff_si(t, d, 1);
                fmpz_poly_set_coeff_si(t, 0, -1);
                fmpz_poly_mul(poly, poly, t);
            }
        }
    }

    for (d = 1; d <= n; d++)
    {
        if (n % d == 0)
        {
            if (n_moebius_mu(n / d) == -1)
            {
                fmpz_poly_zero(t);
                fmpz_poly_set_coeff_si(t, d, 1);
                fmpz_poly_set_coeff_si(t, 0, -1);
                fmpz_poly_div(poly, poly, t);
            }
        }
    }

    fmpz_poly_clear(t);
}
TEST_FUNCTION_START(fmpz_poly_cyclotomic, state)
{
    fmpz_poly_t A, B;
    slong n;

    for (n = 0; n <= 1000; n++)
    {
        fmpz_poly_init(A);
        fmpz_poly_init(B);

        fmpz_poly_cyclotomic(A, n);
        cyclotomic_naive(B, n);

        if (!fmpz_poly_equal(A, B))
        {
            flint_printf("FAIL: wrong value of Phi_%wd(x)\n", n);
            flint_printf("Computed:\n");
            fmpz_poly_print_pretty(A, "x");
            flint_printf("\n\nExpected:\n");
            fmpz_poly_print_pretty(B, "x");
            flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(A);
        fmpz_poly_clear(B);
    }

    /* We verify the first value that does not fit on 32 bits.
       This exercises the slow path at least on a 32 bit system.
       Testing the 64 bit value is a bit too much to do by default
        as it requires ~2 GB of memory and takes a few minutes. */
    {
        fmpz_t h, ref;

        const ulong nn = UWORD(10163195);
        /* const ulong nn = UWORD(169828113);  64-bit case */

        fmpz_init(h);
        fmpz_init(ref);
        fmpz_set_str(ref, "1376877780831", 10);
        /* fmpz_set_str(ref, "31484567640915734941", 10);  64-bit case */

        fmpz_poly_init(A);
        fmpz_poly_cyclotomic(A, UWORD(10163195));
        fmpz_poly_height(h, A);

        if (!fmpz_equal(h, ref))
        {
            flint_printf("Bad computation of Phi_%wd(x)\n", nn);
            flint_printf("Computed height:\n");
            fmpz_print(h);
            flint_printf("\nExpected height:\n");
            fmpz_print(ref);
            flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(A);
        fmpz_clear(h);
        fmpz_clear(ref);
    }

    TEST_FUNCTION_END(state);
}
