/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"
#include "ulong_extras.h"
#include "profiler.h"

int
main(void)
{
    int i, j;

    FLINT_TEST_INIT(state);

    flint_printf("univar....");
    fflush(stdout);
/*
    const char* vars[] = {"x","y","z","a","b","c","d","e","f","g","h"};



    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, c;
        fmpz_mpoly_univar_t ax, bx, cx;
        timeit_t time;

        fmpz_mpoly_ctx_init(ctx, 3, ORD_LEX);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(c, ctx);

        fmpz_mpoly_set_str_pretty(a, "-6*(y+z)*(z*x^3-z*y^3)", vars, ctx);
        fmpz_mpoly_set_str_pretty(b, "-8*(y+z*x)*(z^2*x^2-z^2*y^2)", vars, ctx);

        printf("A: "); fmpz_mpoly_print_pretty(a, vars, ctx); printf("\n");
        printf("B: "); fmpz_mpoly_print_pretty(b, vars, ctx); printf("\n");
        fmpz_mpoly_gcd_prs(c, a, b, ctx);
        printf("C: "); fmpz_mpoly_print_pretty(c, vars, ctx); printf("\n");

        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(c, ctx);
    }
*/

    for (i = 0; i < 60 * flint_test_multiplier(); i++)
    {
        int ok;
        fmpz_t ac, bc, d;
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, c, g, t;
        ordering_t ord;
        slong nvars, len1, len2, len3, exp_bound1, exp_bound2, exp_bound3;
        slong coeff_bits, exp_bits1, exp_bits2, exp_bits3;

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 3) + 1;

        fmpz_mpoly_ctx_init(ctx, nvars, ord);

        fmpz_init(ac);
        fmpz_init(bc);
        fmpz_init(d);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(c, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(t, ctx);

        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20);
        len3 = n_randint(state, 20);
        exp_bound1 = n_randint(state, 4) + 2;
        exp_bound2 = n_randint(state, 4) + 2;
        exp_bound3 = n_randint(state, 3) + 2;
        coeff_bits = n_randint(state, 3) + 1;
        fmpz_mpoly_randtest(a, state, len1, exp_bound1, coeff_bits, ctx);
        fmpz_mpoly_randtest(b, state, len2, exp_bound2, coeff_bits, ctx);
        fmpz_mpoly_randtest(c, state, len3, exp_bound3, coeff_bits, ctx);

        fmpz_mpoly_mul_johnson(a, a, c, ctx);
        fmpz_mpoly_mul_johnson(b, b, c, ctx);
        fmpz_mpoly_gcd_prs(g, a, b, ctx);
/*
printf("*******\n");
printf("a: "); fmpz_mpoly_print_pretty(a, vars, ctx); printf("\n");
printf("b: "); fmpz_mpoly_print_pretty(b, vars, ctx); printf("\n");
printf("g: "); fmpz_mpoly_print_pretty(g, vars, ctx); printf("\n");
*/

        if (fmpz_mpoly_is_zero(g, ctx))
        {
            if (!fmpz_mpoly_is_zero(a, ctx) || !fmpz_mpoly_is_zero(b, ctx))
            {
                flint_printf("FAIL\ngcd is zero but both inputs are not\ni: %wd\n",i);
                flint_abort();
            }
            continue;
        }

        if (!fmpz_mpoly_divides_monagan_pearce(t, g, c, ctx))
        {
                flint_printf("FAIL\nc doesn't divide the gcd\ni: %wd\n",i);
                flint_abort();
        }
        ok = fmpz_mpoly_divides_monagan_pearce(a, a, g, ctx);
        ok = ok && fmpz_mpoly_divides_monagan_pearce(b, b, g, ctx);
        if (!ok)
        {
            flint_printf("FAIL\ngcd doesn't divide both inputs\ni: %wd\n",i);
            flint_abort();
        }
        /* now check that gcd(a,b)=1 */

        if (fmpz_mpoly_is_zero(a, ctx))
        {
            if (!fmpz_mpoly_equal_si(b, WORD(1), ctx) && !fmpz_mpoly_equal_si(b, -WORD(1), ctx))
            {
                flint_printf("FAIL\nnot rp\ni: %wd\n",i);
                flint_abort();
                
            }
            continue;
        }

        if (fmpz_mpoly_is_zero(b, ctx))
        {
            if (!fmpz_mpoly_equal_si(a, WORD(1), ctx) && !fmpz_mpoly_equal_si(a, -WORD(1), ctx))
            {
                flint_printf("FAIL\nnot rp\ni: %wd\n",i);
                flint_abort();
                
            }
            continue;
        }

        _fmpz_vec_content(ac, a->coeffs, a->length);
        _fmpz_vec_content(bc, b->coeffs, b->length);
        fmpz_gcd(d, ac, bc);
        if (fmpz_cmp_ui(d, WORD(1)) != 0)
        {
            flint_printf("FAIL\nnot rp\ni: %wd\n",i);
            flint_abort();
            
        }

        for (j = 0; j < nvars; j++)
        {
            fmpz_mpoly_resultant(t, a, b, j, ctx);
            if (fmpz_mpoly_is_zero(t, ctx))
            {
                flint_printf("FAIL\nresultant is zero\ni: %wd  j: %wd\n",i,j);
                flint_abort();

            }
        }

        fmpz_clear(ac);
        fmpz_clear(bc);
        fmpz_clear(d);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(c, ctx);
        fmpz_mpoly_clear(g, ctx);       
        fmpz_mpoly_clear(t, ctx);
    }



    FLINT_TEST_CLEANUP(state);

    printf("PASS\n");
    return 0;
}

