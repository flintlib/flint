/*
    Copyright (C) 2010, 2011 Sebastian Pancratz

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
#include "long_extras.h"
#include "fmpz.h"
#include "fmpq.h"

int main(void)
{
    int i;
    FLINT_TEST_INIT(state);
    
    flint_printf("mpq_init_set_readonly....");
    fflush(stdout);
    
    

    /* Create some small fmpq rationals, clear the mpq_t */
    for (i = 0; i < 100000; i++)
    {
        fmpq_t f;
        mpq_t z;

        fmpq_init(f);
        fmpq_randtest(f, state, FLINT_BITS - 2);

        flint_mpq_init_set_readonly(z, f);
        flint_mpq_clear_readonly(z);
    }

    /* Create some large fmpq rationals, do not clear the mpq_t */
    for (i = 0; i < 100000; i++)
    {
        fmpq_t f;
        mpq_t z;

        fmpq_init(f);
        fmpq_randtest(f, state, 2 * FLINT_BITS);

        if (COEFF_IS_MPZ(*fmpq_numref(f)) && COEFF_IS_MPZ(*(fmpq_denref(f))))
        {
            flint_mpq_init_set_readonly(z, f);
        }

        fmpq_clear(f);
    }

    /* Create some more fmpq rationals */
    for (i = 0; i < 100000; i++)
    {
        fmpq_t f, g;
        mpq_t z;

        fmpq_init(f);
        fmpq_init(g);
        fmpq_randtest(f, state, 2 * FLINT_BITS);

        flint_mpq_init_set_readonly(z, f);
        fmpq_set_mpq(g, z);

        if (!fmpq_equal(f, g))
        {
            flint_printf("FAIL:\n\n");
            flint_printf("f = "), fmpq_print(f), flint_printf("\n");
            flint_printf("g = "), fmpq_print(g), flint_printf("\n");
            gmp_printf("z = %Qd\n", z);
        }

        flint_mpq_clear_readonly(z);
	fmpq_clear(f);
        fmpq_clear(g);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
