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
    int result;
    printf("derivative....");
    fflush(stdout);

    fmpq_poly_randinit();

    // Check aliasing
    for (ulong i = 0; i < 10000UL; i++)
    {
        fmpq_poly_t a, b;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_randtest(a, n_randint(100), n_randint(200));

        fmpq_poly_derivative(b, a);
        fmpq_poly_derivative(a, a);

        result = (fmpq_poly_equal(a, b));
        if (!result)
        {
            printf("Error:\n");
            fmpq_poly_print(a);
            printf("\n\n");
            fmpq_poly_print(b);
            printf("\n\n");
            abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
    }

    // Check constants have derivative zero
    for (ulong i = 0; i < 10000UL; i++)
    {
        fmpq_poly_t a, b;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_randtest(a, n_randint(2), n_randint(200));

        fmpq_poly_derivative(b, a);

        result = (b->length == 0UL);
        if (!result)
        {
            printf("Error:\n");
            fmpq_poly_print(a);
            printf("\n\n");
            fmpq_poly_print(b);
            printf("\n\n");
            abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
    }

    // Check (f g)' = f' g + f g'
    for (ulong i = 0; i < 10000UL; i++)
    {
        fmpq_poly_t a, b, c, d, lhs, rhs;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(c);
        fmpq_poly_init(d);
        fmpq_poly_init(lhs);
        fmpq_poly_init(rhs);
        fmpq_poly_randtest(a, n_randint(100), n_randint(200));
        fmpq_poly_randtest(b, n_randint(100), n_randint(200));

        fmpq_poly_mul(lhs, a, b);
        fmpq_poly_derivative(lhs, lhs);
        fmpq_poly_derivative(c, a);
        fmpq_poly_derivative(d, b);
        fmpq_poly_mul(c, c, b);
        fmpq_poly_mul(d, a, d);
        fmpq_poly_add(rhs, c, d);

        result = fmpq_poly_equal(lhs, rhs);
        if (!result)
        {
            printf("Error:\n");
            fmpq_poly_print(a);
            printf("\n\n");
            fmpq_poly_print(b);
            printf("\n\n");
            abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(c);
        fmpq_poly_clear(d);
        fmpq_poly_clear(lhs);
        fmpq_poly_clear(rhs);
    }

    fmpq_poly_randclear();

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
