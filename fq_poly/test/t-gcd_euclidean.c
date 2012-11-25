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

    Copyright (C) 2011 William Hart
    Copyright (C) 2012 Andres Goens

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "fq_poly.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;
    flint_randinit(state);

    printf("gcd_euclidean....");
    fflush(stdout);

    /* Check that gcd(a,a) = a (made monic) */

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        long d,len;
        fq_ctx_t ctx;

        fq_poly_t a,b,g;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        len = n_randint(state, 15) + 1;
        fq_ctx_init_conway(ctx, p, d, "a");
        fq_poly_init(a);
        fq_poly_init(b);
        fq_poly_init(g);


        fq_poly_randtest_not_zero(a, state, len, ctx);
        fq_poly_make_monic(b,a,ctx);
        fq_poly_gcd_euclidean(g,a,a,ctx);      

        result = (fq_poly_equal(g, b));
        if (!result)
        {
            printf("FAIL:\n");
            printf("a = "), fq_poly_print_pretty(a, "X", ctx), printf("\n");
            printf("b = "), fq_poly_print_pretty(b, "X", ctx), printf("\n");
            printf("g = "), fq_poly_print_pretty(g, "X", ctx), printf("\n");
            abort();
        }
        
        fq_poly_clear(a);
        fq_poly_clear(b);
        fq_poly_clear(g);
    }
    

    /* 
       Find coprime polys, multiply by another poly 
       and check the GCD is that poly 
    */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        long d,len,j;
        fq_ctx_t ctx;

        fq_poly_t a,b,c,g;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        len = n_randint(state, 15) + 1;
        fq_ctx_init_conway(ctx, p, d, "a");
        fq_poly_init(a);
        fq_poly_init(b);
        fq_poly_init(c);
        fq_poly_init(g);


        fq_poly_randtest_not_zero(a, state, len, ctx);

      
        for(j=0;(j<100* flint_test_multiplier()) && !fq_poly_is_one(g);j++){
        fq_poly_randtest_not_zero(b, state, len, ctx);
        fq_poly_gcd_euclidean(g,a,b,ctx);
        

        }

        if(!fq_poly_is_one(g)){ 
            printf("FAIL:\n");
            printf("could not find coprime polynomials after %ld tries\n",j+1);
            abort();
        }

        for(j=0;(j<100* flint_test_multiplier()) && (c->length<2);j++)
                    fq_poly_randtest_not_zero(c, state, len+2, ctx);

        if(c->length<2){
            printf("FAIL:\n");
            printf("could not find non-unit polynomial after %ld tries\n",j+1);
            abort();
        }

        fq_poly_make_monic(c, c,ctx);
        
        fq_poly_mul(a, a, c,ctx);
        fq_poly_mul(b, b, c,ctx);

        fq_poly_gcd_euclidean(g, a, b,ctx);

        result = (fq_poly_equal(g, c));
        if (!result)
        {
            printf("FAIL:\n");
            printf("a = "), fq_poly_print_pretty(a, "X", ctx), printf("\n");
            printf("b = "), fq_poly_print_pretty(b, "X", ctx), printf("\n");
            printf("c = "), fq_poly_print_pretty(c, "X", ctx), printf("\n");
            printf("g = "), fq_poly_print_pretty(g, "X", ctx), printf("\n");
            abort();
        }
        
        fq_poly_clear(a);
        fq_poly_clear(b);
        fq_poly_clear(c);
        fq_poly_clear(g);
    }

    /* Check aliasing of a and g */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        long d,len;
        fq_ctx_t ctx;

        fq_poly_t a,b,g;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        len = n_randint(state, 15) + 1;
        fq_ctx_init_conway(ctx, p, d, "a");
        fq_poly_init(a);
        fq_poly_init(b);
        fq_poly_init(g);

        fq_poly_randtest_not_zero(a, state, len, ctx);
        fq_poly_randtest_not_zero(b, state, len, ctx);

        fq_poly_gcd_euclidean(g, a, b,ctx);
        fq_poly_gcd_euclidean(a, a, b,ctx);
        
        result = (fq_poly_equal(g, a));
        if (!result)
        {
            printf("FAIL:\n");
            printf("a = "), fq_poly_print_pretty(a, "X", ctx), printf("\n");
            printf("b = "), fq_poly_print_pretty(b, "X", ctx), printf("\n");
            printf("g = "), fq_poly_print_pretty(g, "X", ctx), printf("\n");
            abort();
        }
        
        fq_poly_clear(a);
        fq_poly_clear(b);
        fq_poly_clear(g);
    }


    /* Check aliasing of b and g */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        long d,len;
        fq_ctx_t ctx;

        fq_poly_t a,b,g;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        len = n_randint(state, 15) + 1;
        fq_ctx_init_conway(ctx, p, d, "a");
        fq_poly_init(a);
        fq_poly_init(b);
        fq_poly_init(g);

        fq_poly_randtest_not_zero(a, state, len, ctx);
        fq_poly_randtest_not_zero(b, state, len, ctx);

        fq_poly_gcd_euclidean(g, a, b,ctx);
        fq_poly_gcd_euclidean(b, a, b,ctx);
        
        result = (fq_poly_equal(g, b));
        if (!result)
        {
            printf("FAIL:\n");
            printf("a = "), fq_poly_print_pretty(a, "X", ctx), printf("\n");
            printf("b = "), fq_poly_print_pretty(b, "X", ctx), printf("\n");
            printf("g = "), fq_poly_print_pretty(g, "X", ctx), printf("\n");
            abort();
        }
        
        fq_poly_clear(a);
        fq_poly_clear(b);
        fq_poly_clear(g);
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
