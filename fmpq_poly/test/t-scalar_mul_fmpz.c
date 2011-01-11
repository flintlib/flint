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

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq_poly.h"
#include "long_extras.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("scalar_mul_fmpz....");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing of a and b */
    for (i = 0; i < 10000; i++)
    {
        fmpq_poly_t a, b;
        fmpz_t n;

        fmpz_init(n);
        fmpz_randtest(n, state, 200);

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_randtest(a, state, n_randint(100, state), n_randint(200, state));

        fmpq_poly_scalar_mul_fmpz(b, a, n);
        fmpq_poly_scalar_mul_fmpz(a, a, n);

        result = (fmpq_poly_equal(a, b));
        if (!result)
        {
            printf("FAIL:\n");
            fmpq_poly_print(a), printf("\n\n");
            fmpq_poly_print(b), printf("\n\n");
            abort();
        }

        fmpz_clear(n);
        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
    }

    /* Check that n (a + b) == na + nb */
    for (i = 0; i < 10000; i++)
    {
        fmpq_poly_t a, b, lhs, rhs;
        fmpz_t n;

        fmpz_init(n);
        fmpz_randtest(n, state, 200);

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(lhs);
        fmpq_poly_init(rhs);
        fmpq_poly_randtest(a, state, n_randint(100, state), n_randint(200, state));
        fmpq_poly_randtest(b, state, n_randint(100, state), n_randint(200, state));

        fmpq_poly_scalar_mul_fmpz(lhs, a, n);
        fmpq_poly_scalar_mul_fmpz(rhs, b, n);
        fmpq_poly_add(rhs, lhs, rhs);
        fmpq_poly_add(lhs, a, b);
        fmpq_poly_scalar_mul_fmpz(lhs, lhs, n);

        result = (fmpq_poly_equal(lhs, rhs));
        if (!result)
        {
            printf("FAIL:\n");
            fmpq_poly_print(a), printf("\n\n");
            fmpq_poly_print(b), printf("\n\n");
            abort();
        }

        fmpz_clear(n);
        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(lhs);
        fmpq_poly_clear(rhs);
    }

    /* Compare with fmpq_poly_scalar_mul_si */
    for (i = 0; i < 10000; i++)
    {
        fmpq_poly_t a, b;
        fmpz_t n1;
        long n;

        n = z_randtest();
        fmpz_init(n1);
        fmpz_set_si(n1, n);

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_randtest(a, state, n_randint(100, state), n_randint(200, state));

        fmpq_poly_scalar_mul_fmpz(b, a, n1);
        fmpq_poly_scalar_mul_si(a, a, n);

        result = (fmpq_poly_equal(a, b));
        if (!result)
        {
            printf("FAIL:\n");
            fmpq_poly_print(a), printf("\n\n");
            fmpq_poly_print(b), printf("\n\n");
            abort();
        }

        fmpz_clear(n1);
        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
