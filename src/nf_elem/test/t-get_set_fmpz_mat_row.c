/*
    Copyright (C) 2013 William Hart
                  2020 Julian Rüth

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mat.h"
#include "nf_elem.h"

TEST_FUNCTION_START(nf_elem_get_set_fmpz_mat_row, state)
{
    int i, result;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nf_t nf;
        nf_elem_t a, b;
        fmpz_mat_t mat;
        slong rows, j;
        fmpz_t d;

        nf_init_randtest(nf, state, 40, 200);

        rows = n_randint(state, 100) + 1;
        j = n_randint(state, rows);
        fmpz_mat_init(mat, rows, fmpq_poly_degree(nf->pol));

        nf_elem_init(a, nf);
        nf_elem_init(b, nf);

        nf_elem_randtest(a, state, 200, nf);

        fmpz_init(d);

        nf_elem_get_fmpz_mat_row(mat, j, d, a, nf);
        nf_elem_set_fmpz_mat_row(b, mat, j, d, nf);

        result = nf_elem_equal(a, b, nf);
        if (!result)
        {
           flint_printf("FAIL:\n");
           flint_printf("rows = %wd, cols = %wd, j = %wd\n", rows, fmpq_poly_degree(nf->pol), j);
           flint_printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
           flint_printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
           flint_printf("d = "); fmpz_print(d); printf("\n");
           flint_abort();
        }

        nf_elem_clear(a, nf);
        nf_elem_clear(b, nf);

        fmpz_mat_clear(mat);

        fmpz_clear(d);

        nf_clear(nf);
    }

    TEST_FUNCTION_END(state);
}
