/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <float.h>
#include "test_helpers.h"
#include "double_extras.h"
#include "arf.h"
#include "arb.h"
#include "gr.h"
#include "dfloat.h"

/* Soundness of the predicates: a definite answer (not contains zero,
   disjoint, cmp decided, equal) must agree with exact arithmetic on
   the balls; conservative answers are allowed in borderline cases. */

/* [m +/- r] as an exact pair of arf endpoints */
static void
endpoints(arf_t lo, arf_t hi, const double * m, int n, double r)
{
    arf_t t;
    arf_init(t);
    _dfloat_get_arf(lo, m, n);
    arf_set(hi, lo);
    if (r == D_INF)
    {
        arf_neg_inf(lo);
        arf_pos_inf(hi);
    }
    else
    {
        arf_set_d(t, r);
        arf_sub(lo, lo, t, ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_add(hi, hi, t, ARF_PREC_EXACT, ARF_RND_DOWN);
    }
    arf_clear(t);
}

TEST_FUNCTION_START(predicates, state)
{
    slong iter;

    DFLOAT_TEST_SKIP_IF_UNSUPPORTED(state);

    for (iter = 0; iter < 100000 * flint_test_multiplier(); iter++)
    {
        d3b_t x, y;
        arf_t xlo, xhi, ylo, yhi;
        int c, status, r1, exact_cz, exact_disjoint;

        if (n_randint(state, 3))
        {
            d3b_randtest(x, state);
            d3b_randtest(y, state);
        }
        else
        {
            d3b_randtest_special(x, state);
            d3b_randtest_special(y, state);
        }
        /* make touching balls likely */
        if (n_randint(state, 4) == 0)
        {
            d3b_set(y, x);
            y->d[0] += y->rad * (n_randint(state, 2) ? 2.0 : -2.0);
            y->rad = d_randtest(state) * x->rad;
            if (!(fabs(y->d[0]) <= DBL_MAX)) continue;
        }
        if (n_randint(state, 8) == 0)
        {
            x->rad = fabs(x->d[0]) * (1.0 + n_randint(state, 3) * D_EPS);
            if (!(x->rad <= DBL_MAX)) x->rad = D_INF;
        }

        arf_init(xlo); arf_init(xhi); arf_init(ylo); arf_init(yhi);
        endpoints(xlo, xhi, x->d, 3, x->rad);
        endpoints(ylo, yhi, y->d, 3, y->rad);

        exact_cz = (arf_sgn(xlo) <= 0 && arf_sgn(xhi) >= 0);
        r1 = d3b_contains_zero(x);
        if (exact_cz && !r1)
        {
            flint_printf("FAIL: contains_zero unsound\n");
            d3b_print(x); flint_printf("\n");
            flint_abort();
        }

        exact_disjoint = (arf_cmp(xhi, ylo) < 0 || arf_cmp(yhi, xlo) < 0);
        r1 = d3b_overlaps(x, y);
        if (!r1 && !exact_disjoint)
        {
            flint_printf("FAIL: overlaps unsound\n");
            d3b_print(x); flint_printf("\n");
            d3b_print(y); flint_printf("\n");
            flint_abort();
        }

        status = d3b_cmp(&c, x, y);
        if (status == GR_SUCCESS)
        {
            if (c == 0)
            {
                if (!(arf_equal(xlo, xhi) && arf_equal(ylo, yhi) && arf_equal(xlo, ylo)))
                {
                    flint_printf("FAIL: cmp == 0 unsound\n");
                    flint_abort();
                }
            }
            else if ((c < 0 && !(arf_cmp(xhi, ylo) < 0)) || (c > 0 && !(arf_cmp(yhi, xlo) < 0)))
            {
                flint_printf("FAIL: cmp unsound\n");
                d3b_print(x); flint_printf("\n");
                d3b_print(y); flint_printf("\n");
                flint_abort();
            }
        }

        r1 = d3b_equal(x, y);
        if (r1 && !(arf_equal(xlo, xhi) && arf_equal(ylo, yhi) && arf_equal(xlo, ylo)))
        {
            flint_printf("FAIL: equal unsound\n");
            flint_abort();
        }

        r1 = d3b_contains(x, y);
        if (r1 && !(arf_cmp(xlo, ylo) <= 0 && arf_cmp(yhi, xhi) <= 0))
        {
            flint_printf("FAIL: contains unsound\n");
            d3b_print(x); flint_printf("\n");
            d3b_print(y); flint_printf("\n");
            flint_abort();
        }

        r1 = d3b_contains_d(x, y->d[0]);
        if (r1)
        {
            arf_t p;
            arf_init(p);
            arf_set_d(p, y->d[0]);
            if (!(arf_cmp(xlo, p) <= 0 && arf_cmp(p, xhi) <= 0))
            {
                flint_printf("FAIL: contains_d unsound\n");
                flint_abort();
            }
            arf_clear(p);
        }

        arf_clear(xlo); arf_clear(xhi); arf_clear(ylo); arf_clear(yhi);
    }

    TEST_FUNCTION_END(state);
}
