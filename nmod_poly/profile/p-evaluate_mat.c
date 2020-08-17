/*
    Copyright (C) 2015 Ashish Kedia

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
#include "nmod_poly.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "nmod_mat.h"
#include <time.h>

#define MIN_DIMENSION 1
#define MAX_DIMENSION 501
#define DIMENSION_SCALE 5
#define MIN_LENGTH 1
#define MAX_LENGTH 40

int
main()
{
    FLINT_TEST_INIT(state);
    slong dim, len, i;
    int result;
    nmod_mat_t A, B, C;
    nmod_poly_t poly;
    mp_limb_t n = n_randtest_not_zero(state);
    clock_t horner_begin, paterson_begin;
    double horner_time, paterson_time;

    printf
        ("#Dimension\tLength\t\tHoner's\t\t\tPaterson\t\tBetter(0 = Horner, 1 = Paterson)\n");

    for (dim = MIN_DIMENSION; dim <= MAX_DIMENSION; dim += DIMENSION_SCALE)
    {
        nmod_mat_init(A, dim, dim, n);
        nmod_mat_init(B, dim, dim, n);
        nmod_mat_init(C, dim, dim, n);

        do
        {
            nmod_mat_randtest(A, state);
        } while (nmod_mat_is_zero(A));

        for (len = MIN_LENGTH; len <= MAX_LENGTH; len++)
        {
            nmod_poly_init2(poly, n, len);
            for (i = 0; i < len; i++)
            {
                poly->coeffs[i] = n_randint(state, n);
            }
            poly->length = len;

            horner_begin = clock();
            nmod_poly_evaluate_mat_horner(B, poly, A);
            horner_time = (double) (clock() - horner_begin) / CLOCKS_PER_SEC;

            paterson_begin = clock();
            nmod_poly_evaluate_mat_paterson_stockmeyer(C, poly, A);
            paterson_time =
                (double) (clock() - paterson_begin) / CLOCKS_PER_SEC;

            result = horner_time < paterson_time ? 0 : 1;
            flint_printf("%wd\t\t%wd\t\t%lf\t\t%lf\t\t%d\n", dim, len,
                         horner_time, paterson_time, result);

            nmod_poly_clear(poly);
        }

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(C);
    }

    FLINT_TEST_CLEANUP(state);
    return 0;
}
