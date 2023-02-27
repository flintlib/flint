/*=============================================================================

    This file is part of Antic.

    Antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version. See <http://www.gnu.org/licenses/>.

=============================================================================*/
/******************************************************************************

    Copyright (C) 2013 William Hart
                  2020 Julian Rüth

******************************************************************************/

#include <stdio.h>
#include "nf.h"
#include "nf_elem.h"

int
main(void)
{
    int i, j, result;
    flint_rand_t state;

    flint_printf("pow....");
    fflush(stdout);

    flint_randinit(state);

    /* test pow(a, e) = e*e*...*e */
    for (i = 0; i < 100 * antic_test_multiplier(); i++)
    {
        nf_t nf;
        nf_elem_t a, p1, p2;
        slong exp;

        nf_init_randtest(nf, state, 40, 20);
        
        nf_elem_init(a, nf);
        nf_elem_init(p1, nf);
        nf_elem_init(p2, nf);

        nf_elem_randtest(a, state, 20, nf);
        
        exp = n_randint(state, 10);

        nf_elem_pow(p1, a, exp, nf);
        nf_elem_one(p2, nf);

        for (j = 0; j < exp; j++)
           nf_elem_mul(p2, p2, a, nf);
        
        result = (nf_elem_equal(p1, p2, nf));
        if (!result)
        {
           printf("FAIL:\n");
           printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
           printf("p1 = "); nf_elem_print_pretty(p1, nf, "x"); printf("\n");
           printf("p2 = "); nf_elem_print_pretty(p2, nf, "x"); printf("\n");
           flint_printf("exp = %w\n", exp);
           abort();
        }

        nf_elem_clear(a, nf);
        nf_elem_clear(p1, nf);
        nf_elem_clear(p2, nf);
         
        nf_clear(nf);
    }
    
    /* test aliasing a and res */
    for (i = 0; i < 100 * antic_test_multiplier(); i++)
    {
        nf_t nf;
        nf_elem_t a, p1, p2;
        slong exp;

        nf_init_randtest(nf, state, 40, 20);
        
        nf_elem_init(a, nf);
        nf_elem_init(p1, nf);
        nf_elem_init(p2, nf);

        nf_elem_randtest(a, state, 20, nf);
        
        exp = n_randint(state, 10);

        nf_elem_pow(p1, a, exp, nf);
        nf_elem_set(p2, a, nf);
        nf_elem_pow(p2, p2, exp, nf);
        
        result = (nf_elem_equal(p1, p2, nf));
        if (!result)
        {
           printf("FAIL:\n");
           printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
           printf("p1 = "); nf_elem_print_pretty(p1, nf, "x"); printf("\n");
           printf("p2 = "); nf_elem_print_pretty(p2, nf, "x"); printf("\n");
           flint_printf("exp = %w\n", exp);
           abort();
        }

        nf_elem_clear(a, nf);
        nf_elem_clear(p1, nf);
        nf_elem_clear(p2, nf);
         
        nf_clear(nf);
    }
    
    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
