
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

    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2009 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

int
main (void)
{
    int result;

    printf ("pow....");
    fflush (stdout);

    fmpz_poly_randinit ();

    // Check aliasing of a and b
    for (ulong i = 0; i < 2000UL; i++)
    {
        fmpz_poly_t a, b;

        fmpz_poly_init (a);
        fmpz_poly_init (b);
        fmpz_poly_randtest (b, n_randint (10), n_randint (100));

        ulong exp = n_randtest () % 20UL;

        fmpz_poly_pow (a, b, exp);
        fmpz_poly_pow (b, b, exp);

        result = (fmpz_poly_equal (a, b));
        if (!result)
        {
            printf ("Error:\n");
            printf ("exp = %lu\n", exp);
            printf ("a = ");
            fmpz_poly_print (a);
            printf ("\n\n");
            printf ("b = ");
            fmpz_poly_print (b);
            printf ("\n\n");
            abort ();
        }

        fmpz_poly_clear (a);
        fmpz_poly_clear (b);
    }

    // Compare with repeated multiplications by the case
    for (ulong i = 0; i < 2000UL; i++)
    {
        fmpz_poly_t a, b, c;

        fmpz_poly_init (a);
        fmpz_poly_init (b);
        fmpz_poly_init (c);
        fmpz_poly_randtest (b, n_randint (10), n_randint (100));

        ulong exp = n_randtest () % 20UL;

        fmpz_poly_pow (a, b, exp);

        if (exp == 0UL && b->length > 0)
        {
            fmpz_poly_fit_length (c, 1UL);
            fmpz_set_ui (c->coeffs, 1UL);
            _fmpz_poly_set_length (c, 1UL);
        }
        else
        {
            fmpz_poly_set (c, b);
            ulong j;

            for (j = 1; j < exp; j++)
                fmpz_poly_mul (c, c, b);
        }

        result = (fmpz_poly_equal (a, c));
        if (!result)
        {
            printf ("Error:\n");
            printf ("exp = %lu\n", exp);
            printf ("a = ");
            fmpz_poly_print (a);
            printf ("\n\n");
            printf ("c = ");
            fmpz_poly_print (c);
            printf ("\n\n");
            abort ();
        }

        fmpz_poly_clear (a);
        fmpz_poly_clear (b);
        fmpz_poly_clear (c);
    }

    fmpz_poly_randclear ();

    _fmpz_cleanup ();
    printf ("PASS\n");
    return 0;
}
