/*
    Copyright 2009 William Hart
    Copyright 2010, 2011, 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <stdlib.h>
#include "profiler.h"
#include "flint.h"
#include "fmpz_mat.h"
#include "fmpz.h"
#include "ulong_extras.h"

#define RANDOMIZED 1
#define TABLE 2

slong ns_tab[10];
slong bits_tab[10];

#define DEFAULT 0
#define CLASSICAL 1
#define WAKSMAN 2
#define MULTIMOD 3
#define STRASSEN 4

const char * names[] = {"DEFAULT", "CLASSICAL", "WAKSMAN", "MULTIMOD", "STRASSEN"};

double timing(flint_rand_t state, slong n, slong m, slong k, slong bits1, slong bits2, int randtest, int algorithm)
{
    double __attribute__((unused)) tcpu, twall;
    fmpz_mat_t A, B, C;

    fmpz_mat_init(A, m, n);
    fmpz_mat_init(B, n, k);
    fmpz_mat_init(C, m, k);

    if (randtest)
    {
        fmpz_mat_randtest(A, state, bits1);
        fmpz_mat_randtest(B, state, bits2);
    }
    else
    {
        fmpz_mat_randbits(A, state, bits1);
        fmpz_mat_randbits(B, state, bits2);
    }

    TIMEIT_START
    if (algorithm == DEFAULT)
        fmpz_mat_mul(C, A, B);
    else if (algorithm == CLASSICAL)
        fmpz_mat_mul_classical(C, A, B);
    else if (algorithm == WAKSMAN)
        fmpz_mat_mul_waksman(C, A, B);
    else if (algorithm == MULTIMOD)
        fmpz_mat_mul_multi_mod(C, A, B);
    else if (algorithm == STRASSEN)
	    fmpz_mat_mul_strassen(C, A, B);
    TIMEIT_STOP_VALUES(tcpu, twall)

    fmpz_mat_clear(A);
    fmpz_mat_clear(B);
    fmpz_mat_clear(C);

    return twall;
}

void
compare_two(int parameters, int alg1, int alg2)
{
    double t1, t2;
    slong m, n, k, bits1, bits2, iter;
    int randtest;

    flint_rand_t state, state1, state2;
    flint_randinit(state);
    flint_randinit(state1);
    flint_randinit(state2);

    flint_printf("A = %s\n", names[alg1]);
    flint_printf("B = %s\n\n", names[alg2]);

    flint_printf("   m    n    k   bits1 bits2   r            A           B   speedup\n");

    for (iter = 0; iter < 100; iter++)
    {
        if (parameters == RANDOMIZED)
        {
            m = 1 + n_randint(state, 60);

            if (n_randint(state, 2))
            {
                n = 1 + n_randint(state, 60);
                k = 1 + n_randint(state, 60);
            }
            else
            {
                n = m + n_randint(state, 2);
                k = m + n_randint(state, 2);
            }

            bits1 = 2 + n_randint(state, 1 << (5 + n_randint(state, 10)));

            if (n_randint(state, 2))
                bits2 = bits1;
            else
                bits2 = 2 + n_randint(state, 1 + bits1);

            randtest = (n_randint(state, 4) == 0);
        }
        else
        {
            m = n = k = ns_tab[iter / 10];
            bits1 = bits2 = bits_tab[iter % 10];
            randtest = 0;
        }

        t1 = timing(state1, m, n, k, bits1, bits2, randtest, alg1);
        t2 = timing(state2, m, n, k, bits1, bits2, randtest, alg2);

        flint_printf("%4wd %4wd %4wd   %5wd %5wd   %d   %6.8f  %6.8f     %.3f\n", m, n, k, bits1, bits2, randtest, t1, t2, t1 / t2);
    }

    flint_randclear(state);
    flint_randclear(state1);
    flint_randclear(state2);
}

int main(int argc, char * argv[])
{
    int params = RANDOMIZED;

    if (argc == 5)
    {
        slong i, a, b;

        a = atol(argv[1]);
        b = atol(argv[2]);
        for (i = 0; i < 10; i++)
        {
            ns_tab[i] = round(exp(((9-i) * log(a) + i * log(b))/9));
            flint_printf("%wd\n", ns_tab[i]);
        }

        a = atol(argv[3]);
        b = atol(argv[4]);
        for (i = 0; i < 10; i++)
        {
            bits_tab[i] = round(exp(((9-i) * log(a) + i * log(b))/9));
            flint_printf("%wd\n", bits_tab[i]);
        }

        params = TABLE;
    }

    compare_two(params, DEFAULT, WAKSMAN);

    return 0;
}
