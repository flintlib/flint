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
#include "nf_elem.h"

TEST_FUNCTION_START(nf_elem_set_equal_si_ui, state)
{
    int i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f;
        ulong m, n;
        slong sm, sn;
        int f_is_m;

        nf_t nf;
        nf_elem_t a;
        nf_elem_t b;

        nf_init_randtest(nf, state, 20, 200);

        nf_elem_init(a, nf);
        nf_elem_init(b, nf);

        fmpq_poly_init(f);

        /* unsigned words */
        m = n_randtest(state);
        n = n_randtest(state);
        nf_elem_set_ui(a, m, nf);
        nf_elem_set_ui(b, n, nf);

        if (!nf_elem_equal_ui(a, m, nf) ||
            !nf_elem_equal_ui(b, n, nf) ||
            nf_elem_equal_ui(a, n, nf) != (m == n) ||
            nf_elem_equal_ui(b, m, nf) != (m == n))
        {
            flint_printf("set_ui/equal_ui wrong with\n");
            flint_printf("nf = "); nf_print(nf); flint_printf("\n");
            flint_printf("m = %wu\n", m);
            flint_printf("n = %wu\n", n);
            flint_printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
            flint_printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
            flint_abort();
        }

        fmpq_poly_randtest(f, state, fmpq_poly_degree(nf->pol) - 1, 200);
        nf_elem_set_fmpq_poly(a, f, nf);

        if (fmpq_poly_length(f) == 0)
            f_is_m = m == 0;
        else
            f_is_m = fmpq_poly_length(f) == 1 &&
                     fmpz_equal_ui(fmpq_poly_numref(f), m) &&
                     fmpz_is_one(fmpq_poly_denref(f));

        if (nf_elem_equal_ui(a, m, nf) != f_is_m)
        {
            flint_printf("equal_ui wrong with\n");
            flint_printf("nf = "); nf_print(nf); flint_printf("\n");
            flint_printf("m = %wu\n", m);
            flint_printf("f = "); fmpq_poly_print_pretty(f, "x"); flint_printf("\n");
            flint_printf("a = "); nf_elem_print_pretty(a, nf, "x"); flint_printf("\n");
            flint_abort();
        }

        /* with signed words */
        sm = (slong) m;
        sn = (slong) n;

        nf_elem_set_si(a, sm, nf);
        nf_elem_set_si(b, sn, nf);

        if (!nf_elem_equal_si(a, sm, nf) ||
            !nf_elem_equal_si(b, sn, nf) ||
            nf_elem_equal_si(a, sn, nf) != (m == n) ||
            nf_elem_equal_si(b, sm, nf) != (m == n))
        {
            flint_printf("set_si/equal_si wrong with\n");
            flint_printf("nf = "); nf_print(nf); flint_printf("\n");
            flint_printf("sm = %wd\n", sm);
            flint_printf("sn = %wd\n", sn);
            flint_printf("a = "); nf_elem_print_pretty(a, nf, "x"); flint_printf("\n");
            flint_printf("b = "); nf_elem_print_pretty(b, nf, "x"); flint_printf("\n");
            flint_abort();
        }

        nf_elem_set_fmpq_poly(a, f, nf);
        if (sm == 0)
            f_is_m = fmpq_poly_length(f) == 0;
        else
            f_is_m = fmpq_poly_length(f) == 1 &&
                     fmpz_equal_si(fmpq_poly_numref(f), sm) &&
                     fmpz_is_one(fmpq_poly_denref(f));

        if (nf_elem_equal_si(a, sm, nf) != f_is_m)
        {
            flint_printf("equal_si wrong with\n");
            flint_printf("nf = "); nf_print(nf); flint_printf("\n");
            flint_printf("sm = %wd\n", sm);
            flint_printf("f = "); fmpq_poly_print_pretty(f, "x"); flint_printf("\n");
            flint_printf("a = "); nf_elem_print_pretty(a, nf, "x"); flint_printf("\n");
            flint_abort();
        }

        /* cleaning */

        fmpq_poly_clear(f);
        nf_elem_clear(a, nf);
        nf_elem_clear(b, nf);
        nf_clear(nf);
    }

    TEST_FUNCTION_END(state);
}
