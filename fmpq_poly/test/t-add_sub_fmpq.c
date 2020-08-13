/*
     Copyright (C) 2020 Vincent Delecroix

     This file is part of FLINT.

     FLINT is free software: you can redistribute it and/or modify it under
     the terms of the GNU Lesser General Public License (LGPL) as published
     by the Free Software Foundation; either version 2.1 of the License, or
     (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "flint.h"
#include "fmpq_poly.h"

int
main(void)
{
    int i;

    FLINT_TEST_INIT(state);

    flint_printf("add_sub_fmpq.c....");
    fflush(stdout);

    /* test b + c - c = b */
    for (i = 0; i < 100; i++)
    {
        fmpq_poly_t a, b, t;
        fmpq_t c;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(t);
        fmpq_init(c);

        fmpq_poly_randtest(b, state, 10, 200);
        fmpq_randtest(c, state, 200);

        fmpq_poly_add_fmpq(t, b, c);
        fmpq_poly_sub_fmpq(a, t, c);

        if (!fmpq_poly_equal(a, b))
        {
           printf("FAIL:\n");
           printf("a = "); fmpq_poly_print(a); printf("\n");
           printf("b = "); fmpq_poly_print(b); printf("\n");
           printf("c = "); fmpq_print(c); printf("\n");
           abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(t);
        fmpq_clear(c);
    }

    /* test aliasing a and b */
    for (i = 0; i < 100; i++)
    {
        fmpq_poly_t a, b;
        fmpq_t c;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_init(c);

        fmpq_poly_randtest(b, state, 10, 200);
        fmpq_randtest(c, state, 200);

        fmpq_poly_set(a, b);
        fmpq_poly_add_fmpq(b, b, c);
        fmpq_poly_sub_fmpq(b, b, c);

        if (!fmpq_poly_equal(a, b))
        {
           printf("FAIL:\n");
           printf("(with aliasing)\n");
           printf("a = "); fmpq_poly_print(a); printf("\n");
           printf("b = "); fmpq_poly_print(b); printf("\n");
           printf("c = "); fmpq_print(c); printf("\n");
           abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_clear(c);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");

    return 0;
}
