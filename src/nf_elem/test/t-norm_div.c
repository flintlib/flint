/*=============================================================================

    This file is part of Antic.

    Antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version. See <http://www.gnu.org/licenses/>.

=============================================================================*/
/******************************************************************************

    Copyright (C) 2014 William Hart

******************************************************************************/

#include <stdio.h>
#include "nf.h"
#include "nf_elem.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    flint_printf("norm_div....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 10 * antic_test_multiplier(); i++)
    {
        fmpq_poly_t pol;
        fmpz_poly_t pol2;
        nf_t nf;
        nf_elem_t a, b, c;
        fmpq_t anorm, bnorm, cnorm, cnorm2;
        fmpz_t d;
        slong nbits;

        fmpq_poly_init(pol);
        fmpz_poly_init(pol2);
        fmpz_init(d);

        do {
            fmpz_poly_randtest_not_zero(pol2, state, 25, 200);
            fmpz_one(fmpz_poly_lead(pol2));
        } while (fmpz_poly_degree(pol2) < 1);

        fmpq_poly_set_fmpz_poly(pol, pol2);
        nf_init(nf, pol);

        nf_elem_init(a, nf);
        nf_elem_init(b, nf);
        nf_elem_init(c, nf);
        
        fmpq_init(anorm);
        fmpq_init(bnorm);
        fmpq_init(cnorm);
        fmpq_init(cnorm2);

        do {
            nf_elem_randtest(a, state, 200, nf);
            nf_elem_get_den(d, a, nf);
            nf_elem_scalar_mul_fmpz(a, a, d, nf);

            nf_elem_randtest(b, state, 200, nf);
            nf_elem_get_den(d, b, nf);
            nf_elem_scalar_mul_fmpz(b, b, d, nf);

            nf_elem_mul(c, a, b, nf);
            nf_elem_norm(anorm, a, nf);
            nf_elem_norm(bnorm, b, nf);
            nf_elem_norm(cnorm2, c, nf);
        } while (fmpq_is_zero(cnorm2));

        nbits = (slong)fmpz_bits(fmpq_numref(anorm)) + 1;
        nf_elem_norm_div(cnorm, c, nf, fmpq_numref(bnorm), nbits);

        result = (fmpq_equal(cnorm, anorm));
        if (!result)
        {
           printf("FAIL:\n");
           printf("nf->pol = "); fmpq_poly_print_pretty(nf->pol, "x"); printf("\n");
           printf("c = "); nf_elem_print_pretty(c, nf, "x"); printf("\n");
           printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
           printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
           printf("norm(c) = "); fmpq_print(cnorm2); printf("\n");
           printf("norm(c) = "); fmpq_print(anorm); printf(" * "); fmpq_print(bnorm); printf("\n");
           printf("supplied divisor = "); fmpq_print(bnorm); printf("\n");
           printf("divisor found = "); fmpq_print(cnorm); printf("\n");
           printf("correct result = "); fmpq_print(anorm); printf("\n");
           flint_printf("bits bound: %wd\n", nbits);
           abort();
        }

        fmpq_clear(anorm);
        fmpq_clear(bnorm);
        fmpq_clear(cnorm);
        fmpq_clear(cnorm2);

        fmpz_clear(d);

        nf_elem_clear(a, nf);
        nf_elem_clear(b, nf);
        nf_elem_clear(c, nf);
         
        nf_clear(nf);

        fmpq_poly_clear(pol);
        fmpz_poly_clear(pol2);
    }

    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
