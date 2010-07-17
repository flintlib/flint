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
#include "fmpz_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int result;
    printf("gcd_subresultant....");
    fflush(stdout);

    fmpz_poly_randinit();

    // Check aliasing of a and b
    for (ulong i = 0; i < 500UL; i++)
    {
        fmpz_poly_t a, b, c;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_randtest(b, n_randint(50), n_randint(100));
        fmpz_poly_randtest(c, n_randint(50), n_randint(100));

        fmpz_poly_gcd_subresultant(a, b, c);
        fmpz_poly_gcd_subresultant(b, b, c);

        result = (fmpz_poly_equal(a, b));
        if (!result)
        {
            printf("Error:\n");
            fmpz_poly_print(a);
            printf("\n\n");
            fmpz_poly_print(b);
            printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    // Check aliasing of a and c
    for (ulong i = 0; i < 500UL; i++)
    {
        fmpz_poly_t a, b, c;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_randtest(b, n_randint(50), n_randint(100));
        fmpz_poly_randtest(c, n_randint(50), n_randint(100));

        fmpz_poly_gcd_subresultant(a, b, c);
        fmpz_poly_gcd_subresultant(c, b, c);

        result = (fmpz_poly_equal(a, c));
        if (!result)
        {
            printf("Error:\n");
            fmpz_poly_print(a);
            printf("\n\n");
            fmpz_poly_print(c);
            printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    // Check that a divides GCD(af, ag)
    for (ulong i = 0; i < 500UL; i++)
    {
        fmpz_poly_t a, d, f, g, q, r;

        fmpz_poly_init(a);
        fmpz_poly_init(d);
        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(q);
        fmpz_poly_init(r);
        fmpz_poly_randtest_not_zero(a, n_randint(24) + 1, n_randint(24) + 1);
        fmpz_poly_randtest(f, n_randint(50), n_randint(100));
        fmpz_poly_randtest(g, n_randint(50), n_randint(100));
        
        fmpz_poly_mul(f, a, f);
        fmpz_poly_mul(g, a, g);
        fmpz_poly_gcd_subresultant(d, f, g);
        
        fmpz_poly_divrem_divconquer(q, r, d, a);
        
        result = (r->length == 0UL);
        if (!result)
        {
            printf("Error:\n");
            fmpz_poly_print(f), printf("\n");
            fmpz_poly_print(g), printf("\n");
            fmpz_poly_print(d), printf("\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(d);
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(q);
        fmpz_poly_clear(r);
    }

    fmpz_poly_randclear();

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
