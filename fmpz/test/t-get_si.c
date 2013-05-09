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

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "long_extras.h"

int
main(void)
{
    fmpz_t x;

    int i, result;
    flint_rand_t state;

    flint_randinit(state);

    printf("get/set_si....");
    fflush(stdout);

    fmpz_init(x);

    fmpz_set_si(x, COEFF_MIN);
    if (COEFF_IS_MPZ(*x) || fmpz_get_si(x) != COEFF_MIN)
    {
        printf("FAIL: COEFF_MIN");
        abort();
    }

    fmpz_set_si(x, COEFF_MAX);
    if (COEFF_IS_MPZ(*x) || fmpz_get_si(x) != COEFF_MAX)
    {
        printf("FAIL: COEFF_MIN");
        abort();
    }

    fmpz_set_si(x, LONG_MIN);
    if (!COEFF_IS_MPZ(*x) || fmpz_get_si(x) != LONG_MIN)
    {
        printf("FAIL: LONG_MIN");
        abort();
    }

    fmpz_set_si(x, LONG_MIN);
    if (!COEFF_IS_MPZ(*x) || fmpz_get_si(x) != LONG_MIN)
    {
        printf("FAIL: LONG_MAX");
        abort();
    }

    fmpz_clear(x);

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a;
        len_t b, c;

        b = z_randtest(state);

        fmpz_init(a);

        fmpz_set_si(a, b);
        c = fmpz_get_si(a);

        result = (b == c);
        if (!result)
        {
            printf("FAIL:\n");
            printf("b = %ld, c = %ld\n", b, c);
            abort();
        }

        fmpz_clear(a);
    }

    flint_randclear(state);

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
