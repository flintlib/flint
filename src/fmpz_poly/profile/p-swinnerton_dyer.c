/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2026 Mael Hostettler
    
    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include "fmpz.h"
#include "fmpz_poly.h"
#include "arb_poly.h"
#include "profiler.h"

/*
    Compares fmpz_poly_swinnerton_dyer (iterated resultants, exact integer
    arithmetic) with the previous implementation (arb ball arithmetic,
    rounding the coefficients).

    Usage:  p-swinnerton_dyer [-n min [max]] [-maxarb n] [-nocheck]

        -n min [max]   range of n to time (default 1 .. 13)
        -maxarb n      do not time the arb version beyond this n (default 13)
        -nocheck       do not compare the two outputs
*/

/* The previous implementation: ball arithmetic with rounding. */
static void
swinnerton_dyer_arb(fmpz_poly_t poly, ulong n)
{
    slong N = (WORD(1) << n);

    fmpz_poly_fit_length(poly, N + 1);

    if (n == 0)
    {
        fmpz_zero(poly->coeffs);
        fmpz_one(poly->coeffs + 1);
    }
    else
    {
        arb_poly_t t;
        arb_poly_init(t);
        arb_poly_swinnerton_dyer_ui(t, n, 0);
        if (!_arb_vec_get_unique_fmpz_vec(poly->coeffs, t->coeffs, t->length))
            flint_throw(FLINT_ERROR, "arb reference failed for n = %wu\n", n);
        arb_poly_clear(t);
    }

    _fmpz_poly_set_length(poly, N + 1);
}

int
main(int argc, char * argv[])
{
    fmpz_poly_t S, T;
    slong n, k, nmin, nmax, maxarb;
    int check;
    double FLINT_SET_BUT_UNUSED(cpu_new), wall_new;
    double FLINT_SET_BUT_UNUSED(cpu_arb), wall_arb;

    nmin = 1;
    nmax = 15;
    maxarb = 15;
    check = 1;

    for (k = 1; k < argc; k++)
    {
        if (strcmp(argv[k], "-n") == 0 && k + 1 < argc)
        {
            nmin = nmax = atoi(argv[k + 1]);
            if (k + 2 < argc && argv[k + 2][0] != '-')
                nmax = atoi(argv[k + 2]);
        }
        else if (strcmp(argv[k], "-maxarb") == 0 && k + 1 < argc)
            maxarb = atoi(argv[k + 1]);
        else if (strcmp(argv[k], "-nocheck") == 0)
            check = 0;
    }

    fmpz_poly_init(S);
    fmpz_poly_init(T);

    flint_printf("Swinnerton-Dyer polynomials, n = %wd .. %wd\n", nmin, nmax);
    flint_printf("resultant = iterated resultants (new), arb = ball arithmetic (old)\n");
    flint_printf("times are wall-clock seconds\n\n");
    flint_printf("%3s %9s %12s %12s %9s  %s\n",
        "n", "degree", "resultant", "arb", "arb/new", "check");

    for (n = nmin; n <= nmax; n++)
    {
        flint_printf("%3wd %9wd ", n, (WORD(1) << n));
        fflush(stdout);

        TIMEIT_START
            fmpz_poly_swinnerton_dyer(S, n);
        TIMEIT_STOP_VALUES(cpu_new, wall_new);

        flint_printf("%12.6f ", wall_new);
        fflush(stdout);

        if (n <= maxarb)
        {
            TIMEIT_START
                swinnerton_dyer_arb(T, n);
            TIMEIT_STOP_VALUES(cpu_arb, wall_arb);

            flint_printf("%12.6f %8.2fx  ", wall_arb, wall_arb / wall_new);

            if (check)
                flint_printf("%s", fmpz_poly_equal(S, T) ? "ok" : "MISMATCH");
            else
                flint_printf("-");
        }
        else
        {
            flint_printf("%12s %9s  -", "-", "-");
        }

        flint_printf("\n");
        fflush(stdout);

        if (check && n <= maxarb && !fmpz_poly_equal(S, T))
            flint_abort();
    }

    fmpz_poly_clear(S);
    fmpz_poly_clear(T);
    flint_cleanup_master();
    return 0;
}
