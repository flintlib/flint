/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"
#include "profiler.h"

#define lenlo  2
#define lenhi  2
#define lenh   1
#define bits   16
#define elo    26
#define ehi    50
#define eh     1
#define cpumin 100

int
main(void)
{
    int len, e;
    fmpz_poly_t f, g;

    FLINT_TEST_INIT(state);


    fmpz_poly_init2(f, lenhi);
    fmpz_poly_init2(g, ehi * (lenhi - 1) + 1);

    flint_printf("| len | exp | binomial |\n");

    for (len = lenlo; len <= lenhi; len += lenh)
    {
        /*
           Construct random polynomial f of length len
         */
        {
            slong k;
            for (k = 0; k < len; k++)
                fmpz_randbits(f->coeffs + k, state, bits);
            if ((f->coeffs)[len-1] == WORD(0))
                fmpz_randtest_not_zero(f->coeffs + (len-1), state, bits);
            f->length = len;
        }

        for (e = elo; e <= ehi; e += eh)
        {
            timeit_t t[1];
            int l, loops = 1, r = 0;
            slong s[1] = {0};

          loop:

            timeit_start(t[0]);
            for (l = 0; l < loops; l++)
                fmpz_poly_pow_binomial(g, f, e);
            timeit_stop(t[0]);

            if (t[0]->cpu <= cpumin)
            {
                loops *= 10;
                goto loop;
            }

            s[0] += t[0]->cpu;
            r    += loops;

            flint_printf("%d %d %lf\n", len, e, s[0] / (double) r);
        }
    }

    fmpz_poly_clear(f);
    fmpz_poly_clear(g);

    FLINT_TEST_CLEAR(state);

    return 0;
}
