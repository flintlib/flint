/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "fmpz_poly.h"
#include "profiler.h"

int
main(int argc, char * argv[])
{
    fmpz_poly_t f, g;
    fmpz_t c, d;
    slong bits, len, k, nthreads;
    flint_rand_t state;

    flint_randinit(state);

    fmpz_poly_init(f);
    fmpz_poly_init(g);
    fmpz_init(c);
    fmpz_init(d);

    fmpz_one(d);

    nthreads = 1;

    for (k = 0; k < argc - 1; k++)
    {
        if (strcmp(argv[k], "-threads") == 0)
        {
            nthreads = atoi(argv[k + 1]);
        }
    }

    flint_printf("using up to %wd threads\n\n", nthreads);
    flint_set_num_threads(nthreads);

    for (len = 32; len < 100000; len = len * 1.25)
    {
        for (bits = 32; bits < 200000; bits *= 4)
        {
            fmpz_poly_zero(f);

            for (k = 0; k < len; k++)
            {
                fmpz_randbits(c, state, bits);
                fmpz_poly_set_coeff_fmpz(f, k, c);
            }

            flint_printf("%wd  %wd  h     ", len, bits);
            TIMEIT_START
            fmpz_poly_taylor_shift_horner(g, f, d);
            TIMEIT_STOP

            flint_printf("%wd  %wd  dc    ", len, bits);
            TIMEIT_START
            fmpz_poly_taylor_shift_divconquer(g, f, d);
            TIMEIT_STOP

            for (k = 1; k <= nthreads; k *= 2)
            {
                flint_printf("%wd  %wd  mm%wd   ", len, bits, k);
                flint_set_num_threads(k);
                TIMEIT_START
                fmpz_poly_taylor_shift_multi_mod(g, f, d);
                TIMEIT_STOP
                flint_set_num_threads(1);
            }

            flint_printf("\n");
        }
    }

    fmpz_poly_clear(f);
    fmpz_poly_clear(g);
    fmpz_clear(c);
    fmpz_clear(d);

    flint_randclear(state);  
    flint_cleanup();
    return 0;
}

