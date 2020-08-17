/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"
#include "fmpz_vec.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);

    flint_printf("get_mpfr....");
    fflush(stdout);

    for (i = 0; i < 10000; i++)
    {
        fmpq_t x; mpq_t y;
        mpfr_t f1, f2;
        int r1, r2;
        mpfr_rnd_t rnd;

        fmpq_init(x);
        mpq_init(y);

        mpfr_init2(f1, 2 + n_randint(state, 1000));
        mpfr_init2(f2, mpfr_get_prec(f1));

        fmpq_randtest(x, state, 2 + n_randint(state, 1000));
        fmpq_get_mpq(y, x);

        switch (n_randint(state, 5))
        {
            case 0: rnd = MPFR_RNDZ; break;
            case 1: rnd = MPFR_RNDD; break;
            case 2: rnd = MPFR_RNDU; break;
            case 3: rnd = MPFR_RNDA; break;
            default: rnd = MPFR_RNDN; break;
        }

        r1 = fmpq_get_mpfr(f1, x, rnd);
        r2 = mpfr_set_q(f2, y, rnd);

        if (r1 != r2 || !mpfr_equal_p(f1, f2))
        {
            flint_printf("FAIL\n");
            fmpq_print(x); flint_printf("\n");
            flint_printf("r1 = %d, r2 = %d\n", r1, r2);
            flint_printf("\nf1: \n");
            mpfr_out_str(stdout, 10, 0, f1, MPFR_RNDN);
            flint_printf("\n\nf1: \n");
            mpfr_out_str(stdout, 10, 0, f2, MPFR_RNDN);
            flint_printf("\n\n");
            abort();
        }

        fmpq_clear(x);
        mpq_clear(y);
        mpfr_clear(f1);
        mpfr_clear(f2);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
