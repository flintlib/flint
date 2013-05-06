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
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq_poly.h"
#include "fmpq_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;
    ulong cflags = 0UL;

    printf("cosh_series....");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing of a and c */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b;
        long n = n_randint(state, 50) + 1;

        fmpq_poly_init(a);
        fmpq_poly_init(b);

        fmpq_poly_randtest_not_zero(a, state, n_randint(state, 50) + 1, 50);
        fmpq_poly_set_coeff_ui(a, 0, 0UL);

        fmpq_poly_canonicalise(a);

        fmpq_poly_cosh_series(b, a, n);
        fmpq_poly_cosh_series(a, a, n);

        cflags |= fmpq_poly_is_canonical(a) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(b) ? 0 : 2;
        result = (fmpq_poly_equal(a, b) && !cflags);
        if (!result)
        {
            printf("FAIL:\n");
            fmpq_poly_debug(a), printf("\n\n");
            fmpq_poly_debug(b), printf("\n\n");
            printf("cflags = %lu\n\n", cflags);
            abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
    }

    /* Check cosh(A)^2-1 = sinh(A)^2 */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t A, coshA, sinhA, B, C, one;
        long n = n_randint(state, 80) + 1;

        fmpq_poly_init(A);
        fmpq_poly_init(coshA);
        fmpq_poly_init(sinhA);
        fmpq_poly_init(B);
        fmpq_poly_init(C);
        fmpq_poly_init(one);

        fmpq_poly_randtest_not_zero(A, state, n_randint(state, 60) + 1, 80);
        fmpq_poly_set_coeff_ui(A, 0, 0UL);

        fmpq_poly_cosh_series(coshA, A, n);
        fmpq_poly_sinh_series(sinhA, A, n);
        fmpq_poly_mullow(B, coshA, coshA, n);
        fmpq_poly_set_coeff_ui(one, 0, 1UL);
        fmpq_poly_sub(B, B, one);
        fmpq_poly_mullow(C, sinhA, sinhA, n);

        cflags |= fmpq_poly_is_canonical(coshA) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(sinhA) ? 0 : 2;
        result = (fmpq_poly_equal(B, C) && !cflags);
        if (!result)
        {
            printf("FAIL:\n");
            printf("A = "), fmpq_poly_debug(A), printf("\n\n");
            printf("cosh(A) = "), fmpq_poly_debug(coshA), printf("\n\n");
            printf("sinh(A) = "), fmpq_poly_debug(sinhA), printf("\n\n");
            printf("cflags = %lu\n\n", cflags);
            abort();
        }

        fmpq_poly_clear(A);
        fmpq_poly_clear(coshA);
        fmpq_poly_clear(sinhA);
        fmpq_poly_clear(B);
        fmpq_poly_clear(C);
        fmpq_poly_clear(one);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
