/*=============================================================================

    This file is part of Antic.

    Antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version. See <http://www.gnu.org/licenses/>.

=============================================================================*/
/******************************************************************************

    Copyright (C) 2020 Vincent Delecroix

******************************************************************************/

#include <stdio.h>
#include "nf.h"
#include "nf_elem.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    flint_randinit(state);

    /* test b + c - c = b */
    for (i = 0; i < 100; i++)
    {
        nf_t nf;
        nf_elem_t a, b, t;
        fmpq_t c;

        nf_init_randtest(nf, state, 40, 200);

        nf_elem_init(a, nf);
        nf_elem_init(b, nf);
        nf_elem_init(t, nf);
        fmpq_init(c);

        nf_elem_randtest(b, state, 200, nf);
        fmpq_randtest_not_zero(c, state, 200);

        nf_elem_scalar_mul_fmpq(t, b, c, nf);
        nf_elem_scalar_div_fmpq(a, t, c, nf);

        result = nf_elem_equal(a, b, nf);
        if (!result)
        {
           printf("FAIL:\n");
           printf("nf = "); nf_print(nf); printf("\n");
           printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
           printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
           printf("c = "); fmpq_print(c); printf("\n");
           abort();
        }

        nf_elem_clear(a, nf);
        nf_elem_clear(b, nf);
        nf_elem_clear(t, nf);
        fmpq_clear(c);

        nf_clear(nf);
    }

    /* test aliasing a and b */
    for (i = 0; i < 100; i++)
    {
        nf_t nf;
        nf_elem_t a, b;
        fmpq_t c;

        nf_init_randtest(nf, state, 40, 200);

        nf_elem_init(a, nf);
        nf_elem_init(b, nf);
        fmpq_init(c);

        nf_elem_randtest(b, state, 200, nf);
        fmpq_randtest_not_zero(c, state, 200);

        nf_elem_set(a, b, nf);
        nf_elem_scalar_mul_fmpq(b, b, c, nf);
        nf_elem_scalar_div_fmpq(b, b, c, nf);

        result = nf_elem_equal(a, b, nf);
        if (!result)
        {
           printf("FAIL:\n");
           printf("(with aliasing)\n");
           printf("nf = "); nf_print(nf); printf("\n");
           printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
           printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
           printf("c = "); fmpq_print(c); printf("\n");
           abort();
        }

        nf_elem_clear(a, nf);
        nf_elem_clear(b, nf);
        fmpq_clear(c);

        nf_clear(nf);
    }

    flint_randclear(state);

    return 0;
}
