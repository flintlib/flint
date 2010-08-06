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
#include "ulong_extras.h"

int
main(void)
{
    int i, result;

    printf("scalar_mul_ui....");
    fflush(stdout);

    fmpq_poly_randinit();

    /* Check aliasing of a and b */
    for (i = 0; i < 10000; i++)
    {
        fmpq_poly_t a, b;
        ulong n = n_randtest();

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_randtest(a, n_randint(100), n_randint(200));

        fmpq_poly_scalar_mul_ui(b, a, n);
        fmpq_poly_scalar_mul_ui(a, a, n);

        result = (fmpq_poly_equal(a, b));
        if (!result)
        {
            printf("FAIL:\n");
            fmpq_poly_print(a), printf("\n\n");
            fmpq_poly_print(b), printf("\n\n");
            abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
    }

    /* Check that (a + b) * n == a * n + b * n */
    for (i = 0; i < 10000; i++)
    {
        fmpq_poly_t a, b, lhs, rhs;
        ulong n = n_randtest();

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(lhs);
        fmpq_poly_init(rhs);
        fmpq_poly_randtest(a, n_randint(100), n_randint(200));
        fmpq_poly_randtest(b, n_randint(100), n_randint(200));

        fmpq_poly_scalar_mul_ui(lhs, a, n);
        fmpq_poly_scalar_mul_ui(rhs, b, n);
        fmpq_poly_add(rhs, lhs, rhs);
        fmpq_poly_add(lhs, a, b);
        fmpq_poly_scalar_mul_ui(lhs, lhs, n);

        result = (fmpq_poly_equal(lhs, rhs));
        if (!result)
        {
            printf("FAIL:\n");
            fmpq_poly_print(a), printf("\n\n");
            fmpq_poly_print(b), printf("\n\n");
            printf("%li\n\n", n);
            fmpq_poly_print(lhs), printf("\n\n");
            fmpq_poly_print(rhs), printf("\n\n");
            abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
    }

    /* Check (a * n1) * n2 = a * (n1 * n2) */
    for (i = 0; i < 10000; i++)
    {
        fmpq_poly_t a, b, c;
        ulong n1 = n_randbits(FLINT_BITS / 2);
        ulong n2 = n_randbits(FLINT_BITS / 2);

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(c);
        fmpq_poly_randtest(a, n_randint(100), n_randint(200));

        fmpq_poly_scalar_mul_ui(b, a, n1);
        fmpq_poly_scalar_mul_ui(c, b, n2);
        fmpq_poly_scalar_mul_ui(b, a, n1 * n2);

        result = (fmpq_poly_equal(b, c));
        if (!result)
        {
            printf("FAIL:\n");
            printf("n1 = %lu, n2 = %lu:\n\n", n1, n2);
            fmpq_poly_print(a), printf("\n\n");
            fmpq_poly_print(b), printf("\n\n");
            fmpq_poly_print(c), printf("\n\n");
            abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(c);
    }

    fmpq_poly_randclear();
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
