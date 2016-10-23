/*
    Copyright (C) 2011 Andy Novocin
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz_poly.h"

#define LONG_FAC_TEST 0 /* run an extra long test */
#define TEST_HARD 0 /* test hard polynomials */

#if LONG_FAC_TEST
#define FAC_MULT 5
#else
#define FAC_MULT 1
#endif

#if TEST_HARD
#include <sys/time.h>

void factor_poly(const char * file_str, const char * name)
{
   FILE * file;
   fmpz_poly_t f;
   fmpz_poly_factor_t fac;
   struct timeval start, stop;
   ulong ms;

   fmpz_poly_init(f);

   file = fopen(file_str, "rw");
       
   fmpz_poly_fread(file, f);
       
   fmpz_poly_factor_init(fac);

   gettimeofday(&start, NULL);
   fmpz_poly_factor(fac, f);
   gettimeofday(&stop, NULL);
   ms = (stop.tv_sec - start.tv_sec)*1000 + (stop.tv_usec - start.tv_usec) / 1000;

   flint_printf("%s has %wd factors\n", name, fac->num);
   flint_printf("Time to factor %s: %ld ms\n\n", name, ms);

   fmpz_poly_factor_clear(fac);
   fclose(file);

   fmpz_poly_clear(f);
}
#endif

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("factor....");
    fflush(stdout);

#if TEST_HARD
    factor_poly("/home/wbhart/.julia/v0.5/Nemo/deps/flint2/fmpz_poly_factor/test/P1_flint", "P1");
    factor_poly("/home/wbhart/.julia/v0.5/Nemo/deps/flint2/fmpz_poly_factor/test/P2_flint", "P2");
    factor_poly("/home/wbhart/.julia/v0.5/Nemo/deps/flint2/fmpz_poly_factor/test/P3_flint", "P3");
    factor_poly("/home/wbhart/.julia/v0.5/Nemo/deps/flint2/fmpz_poly_factor/test/P4_flint", "P4");
    factor_poly("/home/wbhart/.julia/v0.5/Nemo/deps/flint2/fmpz_poly_factor/test/P5_flint", "P5");
    factor_poly("/home/wbhart/.julia/v0.5/Nemo/deps/flint2/fmpz_poly_factor/test/P6_flint", "P6");
    factor_poly("/home/wbhart/.julia/v0.5/Nemo/deps/flint2/fmpz_poly_factor/test/P7_flint", "P7");
    factor_poly("/home/wbhart/.julia/v0.5/Nemo/deps/flint2/fmpz_poly_factor/test/P8_flint", "P8");
    factor_poly("/home/wbhart/.julia/v0.5/Nemo/deps/flint2/fmpz_poly_factor/test/M12_5_flint", "M12_5");
    factor_poly("/home/wbhart/.julia/v0.5/Nemo/deps/flint2/fmpz_poly_factor/test/M12_6_flint", "M12_6");
    factor_poly("/home/wbhart/.julia/v0.5/Nemo/deps/flint2/fmpz_poly_factor/test/T1_flint", "T1");
    factor_poly("/home/wbhart/.julia/v0.5/Nemo/deps/flint2/fmpz_poly_factor/test/T2_flint", "T2");
    factor_poly("/home/wbhart/.julia/v0.5/Nemo/deps/flint2/fmpz_poly_factor/test/T3_flint", "T3");
    factor_poly("/home/wbhart/.julia/v0.5/Nemo/deps/flint2/fmpz_poly_factor/test/H1_flint", "H1");
    factor_poly("/home/wbhart/.julia/v0.5/Nemo/deps/flint2/fmpz_poly_factor/test/S7_flint", "S7");
    factor_poly("/home/wbhart/.julia/v0.5/Nemo/deps/flint2/fmpz_poly_factor/test/S8_flint", "S8");
    factor_poly("/home/wbhart/.julia/v0.5/Nemo/deps/flint2/fmpz_poly_factor/test/S9_flint", "S9");
    factor_poly("/home/wbhart/.julia/v0.5/Nemo/deps/flint2/fmpz_poly_factor/test/S10_flint", "S10");
    factor_poly("/home/wbhart/.julia/v0.5/Nemo/deps/flint2/fmpz_poly_factor/test/H2_flint", "H2");
#endif

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t c;
        fmpz_poly_t f, g, h, t;
        fmpz_poly_factor_t fac;
        slong j, k, n = n_randint(state, 10*FAC_MULT);
        slong facs1 = 0, facs2 = 0;

        fmpz_init(c);
        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(h);
        fmpz_poly_init(t);
        fmpz_poly_factor_init(fac);

        fmpz_randtest_not_zero(c, state, n_randint(state, 10) + 1);
        fmpz_poly_set_fmpz(f, c);

        for (j = 0; j < n; j++)
        {
            do {
               fmpz_poly_randtest(g, state, n_randint(state, 10) + 2, n_randint(state, 100));
            } while (g->length == 0);
            k = 0;
            while (k < g->length && fmpz_is_zero(g->coeffs + k))
               k++;
            facs1 += k; /* count powers of x separately */
            if (k < g->length - 1)
               facs1++; /* rough lower bound of factors of f */
            fmpz_poly_mul(f, f, g);
        }

        fmpz_poly_factor(fac, f);

        fmpz_poly_set_fmpz(h, &fac->c);
        for (j = 0; j < fac->num; j++)
        {
            if (fac->exp[j] == 1)
                fmpz_poly_mul(h, h, fac->p + j);
            else
            {
                fmpz_poly_pow(t, fac->p + j, fac->exp[j]);
                fmpz_poly_mul(h, h, t);
            }
            facs2 += fac->exp[j];
        }

        result = fmpz_poly_equal(f, h) && facs1 <= facs2;
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("facs1 = %wd, facs2 = %wd\n", facs1, facs2);
            flint_printf("f = "), fmpz_poly_print(f), flint_printf("\n\n");
            flint_printf("h = "), fmpz_poly_print(h), flint_printf("\n\n");
            flint_printf("fac = "), fmpz_poly_factor_print(fac), flint_printf("\n\n");
            abort();
        }

        fmpz_clear(c);
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(h);
        fmpz_poly_clear(t);
        fmpz_poly_factor_clear(fac);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
