/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_poly.h"
#include "fmpz_poly_factor.h"
#include "qqbar.h"

TEST_FUNCTION_START(qqbar_real_roots_fmpz_poly, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpz_poly_t f, g;
        qqbar_ptr r1, r2;
        slong i, j, d, n1, n2;
        int flags;

        fmpz_poly_init(f);
        fmpz_poly_init(g);

        flags = 0;

        if (n_randint(state, 2))
        {
            do {
                fmpz_poly_randtest_irreducible(f, state, 2 + n_randint(state, 10), 1 + n_randint(state, 20));
            } while (fmpz_poly_degree(f) < 1);

            if (n_randint(state, 2))
                flags |= QQBAR_ROOTS_IRREDUCIBLE;
        }
        else
        {
            fmpz_poly_randtest(f, state, 1 + n_randint(state, 6), 1 + n_randint(state, 20));
            fmpz_poly_randtest(g, state, 1 + n_randint(state, 5), 1 + n_randint(state, 20));
            fmpz_poly_mul(f, f, g);
            if (n_randint(state, 2))
                fmpz_poly_mul(f, f, g);
            fmpz_poly_randtest(g, state, 1 + n_randint(state, 4), 1 + n_randint(state, 5));
            fmpz_poly_mul(f, f, g);
        }

        if (fmpz_poly_is_zero(f))
            fmpz_poly_set_ui(f, 1);

        if (n_randint(state, 2))
            flags |= QQBAR_ROOTS_UNSORTED;

        d = fmpz_poly_degree(f);

        r1 = _qqbar_vec_init(FLINT_MAX(d, 0));
        r2 = _qqbar_vec_init(FLINT_MAX(d, 0));

        n1 = qqbar_real_roots_fmpz_poly(r1, f, flags);

        qqbar_roots_fmpz_poly(r2, f, flags & QQBAR_ROOTS_IRREDUCIBLE);
        for (i = n2 = 0; i < d; i++)
            if (qqbar_is_real(r2 + i))
                qqbar_swap(r2 + n2++, r2 + i);

        if (n1 != n2)
        {
            flint_printf("FAIL (count)\n");
            flint_printf("f = %{fmpz_poly}\n", f);
            flint_printf("n1 = %wd, n2 = %wd\n", n1, n2);
            flint_abort();
        }

        for (i = 0; i < n1; i++)
        {
            if (!qqbar_is_real(r1 + i) || !arb_is_zero(acb_imagref(QQBAR_ENCLOSURE(r1 + i))))
            {
                flint_printf("FAIL (real)\n");
                flint_abort();
            }

            if (!(flags & QQBAR_ROOTS_UNSORTED) && i > 0 && qqbar_cmp_re(r1 + i - 1, r1 + i) < 0)
            {
                flint_printf("FAIL (sorted)\n");
                flint_printf("f = %{fmpz_poly}\n", f);
                flint_abort();
            }

            /* the complex version lists the real roots first, in the
               same order */
            if (!(flags & QQBAR_ROOTS_UNSORTED) && !qqbar_equal(r1 + i, r2 + i))
            {
                flint_printf("FAIL (equal)\n");
                flint_printf("f = %{fmpz_poly}\n", f);
                qqbar_print(r1 + i); flint_printf("\n");
                qqbar_print(r2 + i); flint_printf("\n");
                flint_abort();
            }

            if (flags & QQBAR_ROOTS_UNSORTED)
            {
                slong c1 = 0, c2 = 0;
                for (j = 0; j < n1; j++)
                {
                    c1 += qqbar_equal(r1 + i, r1 + j);
                    c2 += qqbar_equal(r1 + i, r2 + j);
                }

                if (c1 != c2)
                {
                    flint_printf("FAIL (multiplicity)\n");
                    flint_printf("f = %{fmpz_poly}\n", f);
                    flint_abort();
                }
            }
        }

        _qqbar_vec_clear(r1, FLINT_MAX(d, 0));
        _qqbar_vec_clear(r2, FLINT_MAX(d, 0));
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
    }

    TEST_FUNCTION_END(state);
}
