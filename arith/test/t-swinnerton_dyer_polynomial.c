/*
    Copyright (C) 2011 Fredrik Johansson

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
#include "arith.h"
#include "profiler.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpq_poly.h"

static const mp_limb_t known_values[] =
{
    UWORD(2147483629),
    UWORD(1073742093),
    UWORD(1342248677),
    UWORD(3319936736),
    UWORD(2947821228),
    UWORD(1019513834),
    UWORD(3324951530),
    UWORD(1995039408),
    UWORD(3505683295),
    UWORD(3567639420),
    UWORD(394942914)
};

int main()
{
    fmpz_poly_t S;
    mp_limb_t r;
    slong n;

    FLINT_TEST_INIT(state);

    flint_printf("swinnerton_dyer_polynomial....");
    fflush(stdout);

    for (n = 0; n <= 10; n++)
    {
        fmpz_poly_init(S);
        arith_swinnerton_dyer_polynomial(S, n);
        r = fmpz_poly_evaluate_mod(S, UWORD(2147483629), UWORD(4294967291));

        if (r != known_values[n])
        {
            flint_printf("ERROR: wrong evaluation of S_%wd\n", n);
            abort();
        }

        fmpz_poly_clear(S);
    }

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
