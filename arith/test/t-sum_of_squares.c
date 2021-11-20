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
#include "flint.h"
#include "arith.h"
#include "fmpz_vec.h"

#define N 10

static const fmpz known[N][N] = {
    {1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {1, 2, 0, 0, 2, 0, 0, 0, 0, 2}, 
    {1, 4, 4, 0, 4, 8, 0, 0, 4, 4},
    {1, 6, 12, 8, 6, 24, 24, 0, 12, 30}, 
    {1, 8, 24, 32, 24, 48, 96, 64, 24, 104}, 
    {1, 10, 40, 80, 90, 112, 240, 320, 200, 250}, 
    {1, 12, 60, 160, 252, 312, 544, 960, 1020, 876}, 
    {1, 14, 84, 280, 574, 840, 1288, 2368, 3444, 3542}, 
    {1, 16, 112, 448, 1136, 2016, 3136, 5504, 9328, 12112}, 
    {1, 18, 144, 672, 2034, 4320, 7392, 12672, 22608, 34802}
};

int main(void)
{
    fmpz * r;
    fmpz_t t;
    slong i, j;

    FLINT_TEST_INIT(state);

    flint_printf("sum_of_squares....");
    fflush(stdout);

    r = _fmpz_vec_init(N);
    fmpz_init(t);

    for (i = 0; i < N; i++)
    {
        arith_sum_of_squares_vec(r, i, N);

        for (j = 0; j < N; j++)
        {
            fmpz_set_ui(t, j);
            arith_sum_of_squares(t, i, t);

            if (!fmpz_equal(t, r + j) || !fmpz_equal(t, known[i] + j))
            {
                flint_printf("FAIL:\n");
                flint_printf("i, j = %wd, %wd, r[j] = %wd, r(j) = %wd, "
                    "expected: %wd\n",
                    i, j, fmpz_get_si(r + j), fmpz_get_si(t), known[i][j]);
                fflush(stdout);
                flint_abort();
            }
        }
    }

    _fmpz_vec_clear(r, N);
    fmpz_clear(t);

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
