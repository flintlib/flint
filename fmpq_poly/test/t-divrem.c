
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
    int i, result;

    printf ("divrem....");
    fflush (stdout);

    fmpq_poly_randinit ();

    // Check aliasing of {q,r} and {a,b}
    for (i = 0; i < 2000; i++)
    {
        fmpq_poly_t a, b, q, r;

        fmpq_poly_init (a);
        fmpq_poly_init (b);
        fmpq_poly_init (q);
        fmpq_poly_init (r);
        fmpq_poly_randtest (a, n_randint (50), n_randint (500));
        fmpq_poly_randtest_not_zero (b, n_randint (499) + 1,
                                     n_randint (499) + 1);

        fmpq_poly_divrem (q, r, a, b);
        fmpq_poly_divrem (a, b, a, b);

        result = (fmpq_poly_equal (q, a)) && (fmpq_poly_equal (r, b));
        if (!result)
        {
            printf ("Error:\n");
            printf ("q = "), fmpq_poly_print (q), printf ("\n\n");
            printf ("r = "), fmpq_poly_print (r), printf ("\n\n");
            printf ("a = "), fmpq_poly_print (a), printf ("\n\n");
            printf ("b = "), fmpq_poly_print (b), printf ("\n\n");
            abort ();
        }

        fmpq_poly_clear (a);
        fmpq_poly_clear (b);
        fmpq_poly_clear (q);
        fmpq_poly_clear (r);
    }

    // Check aliasing of {q,r} and {b,a}
    for (i = 0; i < 2000; i++)
    {
        fmpq_poly_t a, b, q, r;

        fmpq_poly_init (a);
        fmpq_poly_init (b);
        fmpq_poly_init (q);
        fmpq_poly_init (r);
        fmpq_poly_randtest (a, n_randint (50), n_randint (500));
        fmpq_poly_randtest_not_zero (b, n_randint (499) + 1,
                                     n_randint (499) + 1);

        fmpq_poly_divrem (q, r, a, b);
        fmpq_poly_divrem (b, a, a, b);

        result = (fmpq_poly_equal (q, b)) && (fmpq_poly_equal (r, a));
        if (!result)
        {
            printf ("Error:\n");
            printf ("q = "), fmpq_poly_print (q), printf ("\n\n");
            printf ("r = "), fmpq_poly_print (r), printf ("\n\n");
            printf ("a = "), fmpq_poly_print (a), printf ("\n\n");
            printf ("b = "), fmpq_poly_print (b), printf ("\n\n");
            abort ();
        }

        fmpq_poly_clear (a);
        fmpq_poly_clear (b);
        fmpq_poly_clear (q);
        fmpq_poly_clear (r);
    }

    // check a = q b + r
    for (i = 0; i < 2000; i++)
    {
        fmpq_poly_t a, b, q, r, rhs;

        fmpq_poly_init (a);
        fmpq_poly_init (b);
        fmpq_poly_init (q);
        fmpq_poly_init (r);
        fmpq_poly_init (rhs);
        fmpq_poly_randtest (a, n_randint (50), n_randint (500));
        fmpq_poly_randtest_not_zero (b, n_randint (499) + 1,
                                     n_randint (499) + 1);

        fmpq_poly_divrem (q, r, a, b);
        fmpq_poly_mul (rhs, q, b);
        fmpq_poly_add (rhs, rhs, r);

        result = fmpq_poly_equal (a, rhs);
        if (!result)
        {
            printf ("Error:\n");
            printf ("a   = "), fmpq_poly_print (a), printf ("\n\n");
            printf ("rhs = "), fmpq_poly_print (rhs), printf ("\n\n");
            abort ();
        }

        fmpq_poly_clear (a);
        fmpq_poly_clear (b);
        fmpq_poly_clear (q);
        fmpq_poly_clear (r);
        fmpq_poly_clear (rhs);
    }

    fmpq_poly_randclear ();

    _fmpz_cleanup ();
    printf ("PASS\n");
    return 0;
}
