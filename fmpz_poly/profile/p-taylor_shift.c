/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <ctype.h>
#include "fmpz_poly.h"
#include "profiler.h"

int
main(int argc, char * argv[])
{
    fmpz_poly_t f, g, h;
    fmpz_t c, d;
    slong bits, len, k, nthreads, minlen, maxlen, minbits, maxbits;
    double incbits, inclen;
    flint_rand_t state;

    flint_randinit(state);

    fmpz_poly_init(f);
    fmpz_poly_init(g);
    fmpz_poly_init(h);
    fmpz_init(c);
    fmpz_init(d);

    fmpz_one(d);

    nthreads = 1;

    minlen = 32;
    maxlen = 100000;
    inclen = 1.5;

    minbits = 32;
    maxbits = 400000;
    incbits = 2.0;

    for (k = 0; k < argc - 1; k++)
    {
        if (strcmp(argv[k], "-len") == 0)
        {
            minlen = maxlen = atoi(argv[k + 1]);
            if (k + 2 < argc && isdigit(argv[k + 2][0])) maxlen = atoi(argv[k + 2]);
            if (k + 3 < argc && isdigit(argv[k + 3][0])) inclen = atof(argv[k + 3]);
        }

        if (strcmp(argv[k], "-bits") == 0)
        {
            minbits = maxbits = atoi(argv[k + 1]);
            if (k + 2 < argc && isdigit(argv[k + 2][0])) maxbits = atoi(argv[k + 2]);
            if (k + 3 < argc && isdigit(argv[k + 3][0])) incbits = atof(argv[k + 3]);
        }

        if (strcmp(argv[k], "-threads") == 0)
            nthreads = atoi(argv[k + 1]);
    }

    flint_printf("len sizes: ");
    for (len = minlen; len <= maxlen; len = FLINT_MAX(len + 1, len * inclen))
        flint_printf("%wd ", len);

    flint_printf("\nbit sizes: ");
    for (bits = minbits; bits <= maxbits; bits = FLINT_MAX(bits + 1, bits * incbits))
        flint_printf("%wd ", bits);

    flint_printf("\nusing up to %wd threads\n\n", nthreads);

    for (len = minlen; len <= maxlen; len = FLINT_MAX(len + 1, len * inclen))
    {
        for (bits = minbits; bits <= maxbits; bits = FLINT_MAX(bits + 1, bits * incbits))
        {
            fmpz_poly_zero(f);

            for (k = 0; k < len; k++)
            {
                fmpz_randbits(c, state, bits);
                fmpz_poly_set_coeff_fmpz(f, k, c);
            }

            flint_printf("%wd  %wd  default ", len, bits);
            TIMEIT_START
            fmpz_poly_taylor_shift(g, f, d);
            TIMEIT_STOP

            /*
            flint_printf("%wd  %wd  compose ", len, bits);
            TIMEIT_START
            fmpz_poly_one(h);
            fmpz_poly_set_coeff_si(h, 1, 1);
            fmpz_poly_compose_divconquer(h, f, h);
            TIMEIT_STOP
            */

            flint_printf("%wd  %wd  horner  ", len, bits);
            TIMEIT_START
            fmpz_poly_taylor_shift_horner(g, f, d);
            TIMEIT_STOP

            for (k = 1; k <= nthreads; k *= 2)
            {
                flint_printf("%wd  %wd  dc%wd     ", len, bits, k);
                flint_set_num_threads(k);
                TIMEIT_START
                fmpz_poly_taylor_shift_divconquer(g, f, d);
                TIMEIT_STOP
                flint_set_num_threads(1);
            }

            for (k = 1; k <= nthreads; k *= 2)
            {
                flint_printf("%wd  %wd  mm%wd     ", len, bits, k);
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
    fmpz_poly_clear(h);
    fmpz_clear(c);
    fmpz_clear(d);

    flint_randclear(state);  
    flint_cleanup();
    return 0;
}

