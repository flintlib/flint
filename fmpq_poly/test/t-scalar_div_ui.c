
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
#include "fmpq_poly.h"
#include "ulong_extras.h"

int
main (void)
{
    int result;

    printf ("scalar_div_ui....");
    fflush (stdout);

    fmpq_poly_randinit ();

    // Check aliasing of a and b
    for (ulong i = 0; i < 10000UL; i++)
    {
        fmpq_poly_t a, b;

        ulong n = n_randtest_not_zero ();

        fmpq_poly_init (a);
        fmpq_poly_init (b);
        fmpq_poly_randtest (a, n_randint (100), n_randint (200));

        fmpq_poly_scalar_div_ui (b, a, n);
        fmpq_poly_scalar_div_ui (a, a, n);

        result = (fmpq_poly_equal (a, b));
        if (!result)
        {
            printf ("Error:\n");
            fmpq_poly_print (a);
            printf ("\n\n");
            fmpq_poly_print (b);
            printf ("\n\n");
            abort ();
        }

        fmpq_poly_clear (a);
        fmpq_poly_clear (b);
    }

    // check (a/n1)/n2 = a/(n1*n2)
    for (ulong i = 0; i < 10000UL; i++)
    {
        fmpq_poly_t a, b, c;

        ulong n1 = n_randbits (FLINT_BITS / 2);

        ulong n2 = n_randbits (FLINT_BITS / 2);

        if (n1 == 0UL)
            n1 = 1UL;
        if (n2 == 0UL)
            n2 = 1UL;

        fmpq_poly_init (a);
        fmpq_poly_init (b);
        fmpq_poly_init (c);
        fmpq_poly_randtest (a, n_randint (100), n_randint (200));

        fmpq_poly_scalar_div_ui (b, a, n1);
        fmpq_poly_scalar_div_ui (c, b, n2);
        fmpq_poly_scalar_div_ui (b, a, n1 * n2);

        result = (fmpq_poly_equal (b, c));
        if (!result)
        {
            printf ("Error n1 = %lu, n2 = %lu:\n", n1, n2);
            fmpq_poly_print (a);
            printf ("\n\n");
            fmpq_poly_print (b);
            printf ("\n\n");
            fmpq_poly_print (c);
            printf ("\n\n");
            abort ();
        }

        fmpq_poly_clear (a);
        fmpq_poly_clear (b);
        fmpq_poly_clear (c);
    }

    fmpq_poly_randclear ();

    _fmpz_cleanup ();
    printf ("PASS\n");
    return 0;
}
