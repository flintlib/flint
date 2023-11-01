/*
    Copyright (C) 2011 Andy Novocin
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_poly_factor.h"

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

   flint_printf("%s has %wd factors: %ld ms\n", name, fac->num, ms);

   fmpz_poly_factor_clear(fac);
   fclose(file);

   fmpz_poly_clear(f);
}
#endif

TEST_FUNCTION_START(fmpz_poly_factor, state)
{
    int i, result;
    int tmul = 100;
#ifdef _WIN32
    tmul = 1;
#endif

#if TEST_HARD
#define MY_DIR "/home/wbhart/.julia/v0.5/Nemo/deps/flint2/fmpz_poly_factor/test/"
    flint_printf("\n");
    factor_poly(MY_DIR"P1_flint", "P1");
    factor_poly(MY_DIR"P2_flint", "P2");
    factor_poly(MY_DIR"P3_flint", "P3");
    factor_poly(MY_DIR"P4_flint", "P4");
    factor_poly(MY_DIR"P5_flint", "P5");
    factor_poly(MY_DIR"P6_flint", "P6");
    factor_poly(MY_DIR"P7_flint", "P7");
    factor_poly(MY_DIR"P8_flint", "P8");
    factor_poly(MY_DIR"M12_5_flint", "M12_5");
    factor_poly(MY_DIR"M12_6_flint", "M12_6");
    factor_poly(MY_DIR"T1_flint", "T1");
    factor_poly(MY_DIR"T2_flint", "T2");
    factor_poly(MY_DIR"T3_flint", "T3");
    factor_poly(MY_DIR"H1_flint", "H1");
    factor_poly(MY_DIR"S7_flint", "S7");
    factor_poly(MY_DIR"S8_flint", "S8");
    factor_poly(MY_DIR"C1_flint", "C1");
/*
    factor_poly(MY_DIR"S9_flint", "S9");
    factor_poly(MY_DIR"S10_flint", "S10");
    factor_poly(MY_DIR"H2_flint", "H2");
*/
#endif

    for (i = 0; i < tmul * flint_test_multiplier(); i++)
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
            fflush(stdout);
            flint_abort();
        }

        for (j = 0; j < fac->num; j++)
        {
            for (k = j + 1; k < fac->num; k++)
            {
                if (fmpz_poly_equal(fac->p + j, fac->p + k))
                {
                    flint_printf("FAIL (repeated factor):\n");
                    flint_printf("facs1 = %wd, facs2 = %wd\n", facs1, facs2);
                    flint_printf("f = "), fmpz_poly_print(f), flint_printf("\n\n");
                    flint_printf("h = "), fmpz_poly_print(h), flint_printf("\n\n");
                    flint_printf("fac = "), fmpz_poly_factor_print(fac), flint_printf("\n\n");
                    fflush(stdout);
                    flint_abort();
                }
            }

            fmpz_poly_content(c, fac->p + j);
            if (!fmpz_is_one(c) || fmpz_sgn((fac->p + j)->coeffs + fmpz_poly_degree(fac->p + j)) < 0)
            {
                flint_printf("FAIL (factor not reduced):\n");
                flint_printf("facs1 = %wd, facs2 = %wd\n", facs1, facs2);
                flint_printf("f = "), fmpz_poly_print(f), flint_printf("\n\n");
                flint_printf("h = "), fmpz_poly_print(h), flint_printf("\n\n");
                flint_printf("fac = "), fmpz_poly_factor_print(fac), flint_printf("\n\n");
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_clear(c);
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(h);
        fmpz_poly_clear(t);
        fmpz_poly_factor_clear(fac);
    }

    /* Regression test: #783 */
    {
        fmpz_poly_t a;
        fmpz_poly_factor_t f;
        fmpz_poly_init(a);

        fmpz_poly_set_str(a, "31  294114975 822759704 601207031 3459410600 6329635407 5561735376 2870497527 -364079376 -984488613 -2855108728 -5168185293 -2678402184 -3199593893 -2740409376 -1003657917 -549998688 -252445155 -80724312 -19101979 -28418280 -2087859 -9732528 -2615547 -159120 -148311 -4680 -1359 2504 73 0 1");

        fmpz_poly_factor_init(f);
        fmpz_poly_factor(f, a);

        fmpz_poly_factor_clear(f);
        fmpz_poly_clear(a);
    }

    TEST_FUNCTION_END(state);
}
