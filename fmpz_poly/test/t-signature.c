/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, j, result;
    FLINT_TEST_INIT(state);

    flint_printf("signature....");
    fflush(stdout);

    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t poly, linear, quadratic, rem;
        fmpz_t lhs, rhs;
        slong nreal, ncomplex, nreal_max, ncomplex_max, r1, r2;
        slong len = n_randint(state, 20) + 1;
        flint_bitcnt_t bits = n_randint(state, 50) + 1;
        
        fmpz_poly_init2(poly, len);
        fmpz_poly_init2(linear, 2);
        fmpz_poly_init2(quadratic, 3);
        fmpz_poly_init2(rem, len);
        linear->length = 2;
        quadratic->length = 3;
        fmpz_init(lhs);
        fmpz_init(rhs);

        ncomplex_max = n_randint(state, len) / 2;
        nreal_max = len - 2 * ncomplex_max;
        ncomplex = 0;
        nreal = 0;

        fmpz_poly_set_coeff_si(poly, 0, 1);
        
        for (j = 0; j < ncomplex_max; j++)
        {
            fmpz * a = quadratic->coeffs + 2;
            fmpz * b = quadratic->coeffs + 1;
            fmpz * c = quadratic->coeffs;
            
            /* Form a quadratic polynomial with complex roots: b^2 < 4ac */
            fmpz_randtest_not_zero(c, state, bits);
            fmpz_randtest(b, state, bits);
            fmpz_randtest_unsigned(a, state, bits);
            
            if (fmpz_sgn(c) < 0)
            {
                fmpz_neg(c, c);
                fmpz_neg(b, b);
            }
            
            fmpz_mul_ui(rhs, c, 4);
            fmpz_mul(lhs, b, b);
            fmpz_add(lhs, lhs, rhs);
            fmpz_fdiv_q(lhs, lhs, rhs);
            fmpz_add(a, a, lhs);
            
            /* If quadratic does not divide poly over Q, set poly *= complex */
            fmpz_poly_pseudo_rem_cohen(rem, poly, quadratic);
            if (rem->length > 0)
            {
                fmpz_poly_mul(poly, poly, quadratic);
                ncomplex++;
            }
        }
        
        for (j = 0; j < nreal_max; j++)
        {
            /* Form a linear polynomial */
            fmpz_randtest(linear->coeffs, state, bits);
            fmpz_randtest_not_zero(linear->coeffs + 1, state, bits);
            
            /* If linear does not divide poly over Q, set poly *= linear */
            fmpz_poly_pseudo_rem_cohen(rem, poly, linear);
            if (rem->length > 0)
            {
                fmpz_poly_mul(poly, poly, linear);
                nreal++;
            }
        }
        
        fmpz_poly_signature(&r1, &r2, poly);

        result = ((r1 == nreal) && (r2 == ncomplex));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("poly   = "), fmpz_poly_print(poly), flint_printf("\n\n");
            flint_printf("r1 r2  = %wd %wd\n\n", r1, r2);
            abort();
        }

        fmpz_poly_clear(poly);
        fmpz_poly_clear(linear);
        fmpz_poly_clear(quadratic);
        fmpz_poly_clear(rem);
        fmpz_clear(lhs);
        fmpz_clear(rhs);
    }

    {
        fmpz_poly_t poly;
        slong r1, r2;
        fmpz_poly_init(poly);
        fmpz_poly_set_str(poly, "6  1 1 1 10 5 1");
        fmpz_poly_signature(&r1, &r2, poly);
        result = ((r1 == 1) && (r2 == 2));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("poly   = "), fmpz_poly_print(poly), flint_printf("\n\n");
            flint_printf("r1 r2  = %wd %wd\n\n", r1, r2);
            abort();
        }
        fmpz_poly_clear(poly);
    }

    for (i = 0; i < 50; i++)
    {
        slong r, s;
        fmpz_poly_t poly;

        fmpz_poly_init(poly);
 
        fmpz_poly_cyclotomic(poly, i + 3);

        fmpz_poly_signature(&r, &s, poly);

        result = (r == 0 && s == (fmpz_poly_length(poly) - 1)/2);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("Cyclotomic(%ld) has signature (%ld, %ld)\n", i + 3, r, s);
            flint_printf("Expected signature (%ld, %ld)\n", 0, (fmpz_poly_length(poly) - 1)/2);
            abort();
        }
        
        fmpz_poly_clear(poly);
    }
    
    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
