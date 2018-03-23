/*
    Copyright (C) 2017 William Hart

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

int
main(void)
{
    int i, j, k, result;
    FLINT_TEST_INIT(state);

    flint_printf("get/set_monomial_fmpz....");
    fflush(stdout);

    /* check fmpz */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f;
        ordering_t ord;
        fmpz_t c, d;
        slong nvars, len, coeff_bits, exp_bits, index;

        fmpz_init(c);
        fmpz_init(d);

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 20) + 1;

        fmpz_mpoly_ctx_init(ctx, nvars, ord);

        fmpz_mpoly_init(f, ctx);

        len = n_randint(state, 50);
        exp_bits = n_randint(state, 100) + 1;
        coeff_bits = n_randint(state, 100);

        fmpz_mpoly_randtest_bits(f, state, len, coeff_bits, exp_bits, ctx);

        for (j = 0; j < 10; j++)
        {
            fmpz ** exp1 = (fmpz **) flint_malloc(nvars*sizeof(fmpz*));
            fmpz ** exp2 = (fmpz **) flint_malloc(nvars*sizeof(fmpz*));

            for (k = 0; k < nvars; k++)
            {
                exp1[k] = (fmpz *) flint_malloc(sizeof(fmpz)); 
                exp2[k] = (fmpz *) flint_malloc(sizeof(fmpz)); 
                fmpz_init(exp1[k]);
                fmpz_init(exp2[k]);
                fmpz_randtest_unsigned(exp1[k], state, 200);
            }

            index = n_randint(state, f->length + 1);

            fmpz_mpoly_set_monomial_fmpz(f, index, exp1, ctx);

            if (!mpoly_monomials_valid_test(f->exps, f->length, f->bits, ctx->minfo))
                flint_throw(FLINT_ERROR, "Polynomial exponents invalid");

            if (mpoly_monomials_overflow_test(f->exps, f->length, f->bits, ctx->minfo))
                flint_throw(FLINT_ERROR, "Polynomial exponents overflow");

            fmpz_mpoly_get_monomial_fmpz(exp2, f, index, ctx);

            result = 1;
            for (k = 0; k < nvars; k++)
            {
                result &= fmpz_equal(exp1[k], exp2[k]);
                fmpz_clear(exp1[k]);
                fmpz_clear(exp2[k]);
                flint_free(exp1[k]); 
                flint_free(exp2[k]); 
            }

            if (!result)
            {
                printf("FAIL\ncheck fmpz, i = %d, j = %d\n", i, j);
                flint_abort();
            }

            flint_free(exp1);
            flint_free(exp2);
       }

       fmpz_mpoly_clear(f, ctx);  
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
