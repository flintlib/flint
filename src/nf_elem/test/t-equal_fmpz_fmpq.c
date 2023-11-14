/*
    Copyright (C) 2018 Vincent Delecroix
                  2020 Julian Rüth

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_poly.h"
#include "nf_elem.h"

TEST_FUNCTION_START(nf_elem_equal_fmpz_fmpq, state)
{
    int i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t z;
        fmpq_t q;
        fmpq_poly_t f;

        nf_t nf;
        nf_elem_t a;

        fmpq_init(q);
        fmpz_init(z);

        nf_init_randtest(nf, state, 20, 200);

        nf_elem_init(a, nf);

        fmpq_poly_init(f);

        fmpq_poly_randtest(f, state, fmpq_poly_degree(nf->pol) - 1, 200);
        nf_elem_set_fmpq_poly(a, f, nf);

        fmpq_poly_get_coeff_fmpq(q, f, 0);
        if (nf_elem_equal_fmpq(a, q, nf) != (fmpq_poly_length(f) <= 1))
        {
                flint_printf("nf_elem_equal_fmpq wrong\n");
                flint_printf("nf = "); nf_print(nf); flint_printf("\n");
                flint_printf("f = "); fmpq_poly_print_pretty(f, "x"); flint_printf("\n");
                flint_printf("a = "); nf_elem_print_pretty(a, nf, "x"); flint_printf("\n");
                flint_printf("q = "); fmpq_print(q); flint_printf("\n");
                flint_abort();
        }

        fmpz_set(z, fmpq_numref(q));
        if (nf_elem_equal_fmpz(a, z, nf) != (fmpq_poly_length(f) <= 1 && fmpz_is_one(fmpq_denref(q))))
        {
                flint_printf("nf_elem_equal_fmpz wrong\n");
                flint_printf("nf = "); nf_print(nf); flint_printf("\n");
                flint_printf("f = "); fmpq_poly_print_pretty(f, "x"); flint_printf("\n");
                flint_printf("a = "); nf_elem_print_pretty(a, nf, "x"); flint_printf("\n");
                flint_printf("z = "); fmpz_print(z); flint_printf("\n");
                flint_abort();
        }

        fmpq_add_si(q, q, 1);
        if (nf_elem_equal_fmpq(a, q, nf))
        {
                flint_printf("nf_elem_equal_fmpq wrong\n");
                flint_printf("nf = "); nf_print(nf); flint_printf("\n");
                flint_printf("f = "); fmpq_poly_print_pretty(f, "x"); flint_printf("\n");
                flint_printf("a = "); nf_elem_print_pretty(a, nf, "x"); flint_printf("\n");
                flint_printf("q = "); fmpq_print(q); flint_printf("\n");
                flint_abort();
        }

        fmpz_add_ui(z, z, 1);
        if (nf_elem_equal_fmpz(a, z, nf))
        {
                flint_printf("nf_elem_equal_fmpz wrong\n");
                flint_printf("nf = "); nf_print(nf); flint_printf("\n");
                flint_printf("f = "); fmpq_poly_print_pretty(f, "x"); flint_printf("\n");
                flint_printf("a = "); nf_elem_print_pretty(a, nf, "x"); flint_printf("\n");
                flint_printf("z = "); fmpz_print(z); flint_printf("\n");
                flint_abort();
        }

        fmpq_poly_clear(f);
        nf_elem_clear(a, nf);
        nf_clear(nf);
    }

    TEST_FUNCTION_END(state);
}
