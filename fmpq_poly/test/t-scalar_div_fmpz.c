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
main (void)
{
    int i, result;

    printf ("scalar_div_fmpz....");
    fflush (stdout);

    fmpq_poly_randinit ();

    // Check aliasing of a and b
    for (i = 0; i < 10000; i++)
    {
        fmpq_poly_t a, b;
        fmpz_t n;

        fmpz_init (n);
        fmpz_randtest_not_zero (n, 200);

        fmpq_poly_init (a);
        fmpq_poly_init (b);
        fmpq_poly_randtest (a, n_randint (100), n_randint (200));

        fmpq_poly_scalar_div_fmpz (b, a, n);
        fmpq_poly_scalar_div_fmpz (a, a, n);

        result = (fmpq_poly_equal (a, b));
        if (!result)
        {
            printf ("FAIL (aliasing):\n");
            fmpq_poly_print (a), printf ("\n\n");
            fmpq_poly_print (b), printf ("\n\n");
            fmpz_print (n);
            abort ();
        }

        fmpz_clear (n);
        fmpq_poly_clear (a);
        fmpq_poly_clear (b);
    }

    // Compare with fmpq_poly_scalar_mul_si
    for (i = 0; i < 10000; i++)
    {
        fmpq_poly_t a, b, c;
        fmpz_t n1;

        fmpz_init (n1);
        long n = (long) n_randtest();

        if (n == 0L)
            n = 1L;
        if (n_randint (2))
            n = -n;
        fmpz_set_si (n1, n);

        fmpq_poly_init (a);
        fmpq_poly_init (b);
        fmpq_poly_init (c);
        fmpq_poly_randtest (a, n_randint (100), n_randint (200));

        fmpq_poly_scalar_div_fmpz (b, a, n1);
        fmpq_poly_scalar_div_si (c, a, n);

        result = (fmpq_poly_equal (b, c));
        if (!result)
        {
            printf ("Error (comparison with _si):\n");
            fmpq_poly_print (a), printf ("\n\n");
            fmpz_print (n1), printf("\n\n");
            fmpq_poly_print (b), printf ("\n\n");
            fmpq_poly_print (c), printf ("\n\n");
            abort ();
        }

        fmpz_clear (n1);
        fmpq_poly_clear (a);
        fmpq_poly_clear (b);
        fmpq_poly_clear (c);
    }

    // Check that (a / n1) / n2 == a / (n1 * n2)
    for (i = 0; i < 10000; i++)
    {
        fmpq_poly_t a, lhs, rhs;
        fmpz_t n1, n2, n;

        fmpz_init (n1);
        fmpz_init (n2);
        fmpz_init (n);
        
        fmpz_randtest_not_zero(n1, n_randint(100) + 1);
        fmpz_randtest_not_zero(n2, n_randint(100) + 1);
        fmpz_mul(n, n1, n2);

        fmpq_poly_init (a);
        fmpq_poly_init (lhs);
        fmpq_poly_init (rhs);
        fmpq_poly_randtest (a, n_randint (100), n_randint (200));

        fmpq_poly_scalar_div_fmpz (lhs, a, n1);
        fmpq_poly_scalar_div_fmpz (lhs, lhs, n2);
        fmpq_poly_scalar_div_fmpz (rhs, a, n);

        result = (fmpq_poly_equal (lhs, rhs));
        if (!result)
        {
            printf ("FAIL (a / n1 / n2):\n");
            fmpq_poly_print (a), printf ("\n\n");
            fmpz_print (n1), printf("\n\n");
            fmpz_print (n2), printf("\n\n");
            fmpz_print (n), printf("\n\n");
            fmpq_poly_print (lhs), printf ("\n\n");
            fmpq_poly_print (rhs), printf ("\n\n");
            abort ();
        }

        fmpz_clear (n1);
        fmpz_clear (n2);
        fmpz_clear (n);
        fmpq_poly_clear (a);
        fmpq_poly_clear (lhs);
        fmpq_poly_clear (rhs);
    }

    // Check that (a + b) / n == a/n + b/n
    for (i = 0; i < 10000; i++)
    {
        fmpq_poly_t a, b, lhs, rhs;
        fmpz_t n;

        fmpz_init (n);
        
        fmpz_randtest_not_zero(n, n_randint(100) + 1);

        fmpq_poly_init (a);
        fmpq_poly_init (b);
        fmpq_poly_init (lhs);
        fmpq_poly_init (rhs);
        fmpq_poly_randtest (a, n_randint (100), n_randint (200));
        fmpq_poly_randtest (b, n_randint (100), n_randint (200));

        fmpq_poly_scalar_div_fmpz (lhs, a, n);
        fmpq_poly_scalar_div_fmpz (rhs, b, n);
        fmpq_poly_add (rhs, lhs, rhs);
        fmpq_poly_add (lhs, a, b);
        fmpq_poly_scalar_div_fmpz (lhs, lhs, n);

        result = (fmpq_poly_equal (lhs, rhs));
        if (!result)
        {
            printf ("FAIL ((a + b) / n):\n");
            fmpq_poly_print (a), printf ("\n\n");
            fmpq_poly_print (b), printf ("\n\n");
            fmpz_print (n), printf("\n\n");
            fmpq_poly_print (lhs), printf ("\n\n");
            fmpq_poly_print (rhs), printf ("\n\n");
            abort ();
        }

        fmpz_clear (n);
        fmpq_poly_clear (a);
        fmpq_poly_clear (b);
        fmpq_poly_clear (lhs);
        fmpq_poly_clear (rhs);
    }

    fmpq_poly_randclear ();
    _fmpz_cleanup ();
    printf ("PASS\n");
    return 0;
}
