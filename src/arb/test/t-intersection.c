/*
    Copyright (C) 2016 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "arb.h"

TEST_FUNCTION_START(arb_intersection, state)
{
    slong iter;

    /* check a containment requirement */
    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x, y, z, w;
        arb_t xy, yz;
        slong pa, pb, pc;
        int alias;

        arb_init(x);
        arb_init(y);
        arb_init(z);
        arb_init(w);
        arb_init(xy);
        arb_init(yz);

        arb_randtest_special(x, state, 200, 10);
        arb_randtest_special(y, state, 200, 10);
        arb_randtest_special(z, state, 200, 10);
        arb_randtest_special(w, state, 200, 10);
        arb_randtest_special(xy, state, 200, 10);
        arb_randtest_special(yz, state, 200, 10);

        pa = 2 + n_randint(state, 200);
        pb = 2 + n_randint(state, 200);
        pc = 2 + n_randint(state, 200);

        arb_union(xy, x, y, pa);
        arb_union(yz, y, z, pb);
        arb_intersection(w, xy, yz, pc);

        if (!arb_contains(w, y))
        {
            flint_printf("FAIL (containment):\n\n");
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_printf("y = "); arb_print(y); flint_printf("\n\n");
            flint_printf("z = "); arb_print(z); flint_printf("\n\n");
            flint_printf("w = "); arb_print(w); flint_printf("\n\n");
            flint_abort();
        }

        if (n_randint(state, 2))
        {
            arb_intersection(xy, xy, yz, pc);
            alias = arb_equal(xy, w);
        }
        else
        {
            arb_intersection(yz, xy, yz, pc);
            alias = arb_equal(yz, w);
        }

        if (!alias)
        {
            flint_printf("FAIL (aliasing):\n\n");
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_printf("y = "); arb_print(y); flint_printf("\n\n");
            flint_printf("z = "); arb_print(z); flint_printf("\n\n");
            flint_printf("w = "); arb_print(w); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x);
        arb_clear(y);
        arb_clear(z);
        arb_clear(w);
        arb_clear(xy);
        arb_clear(yz);
    }

    /* require that the return value is the same as for arb_overlaps */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, y;
        fmpq_t am, ar, bm, br, t, u;
        int c1, c2, c3;
        slong prec;

        prec = 2 + n_randint(state, 200);

        arb_init(a);
        arb_init(b);
        arb_init(y);

        fmpq_init(am);
        fmpq_init(ar);
        fmpq_init(bm);
        fmpq_init(br);
        fmpq_init(t);
        fmpq_init(u);

        arb_randtest(a, state, 1 + n_randint(state, 500), 14);
        arb_randtest(b, state, 1 + n_randint(state, 500), 14);

        arf_get_fmpq(am, arb_midref(a));
        mag_get_fmpq(ar, arb_radref(a));
        arf_get_fmpq(bm, arb_midref(b));
        mag_get_fmpq(br, arb_radref(b));

        fmpq_sub(t, am, bm);
        fmpz_abs(fmpq_numref(t), fmpq_numref(t));
        fmpq_add(u, ar, br);

        c1 = arb_overlaps(a, b);

        c2 = (fmpq_cmp(t, u) <= 0);

        c3 = arb_intersection(y, a, b, prec);

        if (c1 != c2 || c1 != c3)
        {
            flint_printf("FAIL (compatibility with arb_overlaps):\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("y = "); arb_print(y); flint_printf("\n\n");
            flint_printf("am = "); fmpq_print(am); flint_printf("\n\n");
            flint_printf("ar = "); fmpq_print(ar); flint_printf("\n\n");
            flint_printf("bm = "); fmpq_print(bm); flint_printf("\n\n");
            flint_printf("br = "); fmpq_print(br); flint_printf("\n\n");
            flint_printf("t = "); fmpq_print(t); flint_printf("\n\n");
            flint_printf("u = "); fmpq_print(u); flint_printf("\n\n");
            flint_printf("c1 = %d, c2 = %d, c3 = %d\n\n", c1, c2, c3);
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(y);

        fmpq_clear(am);
        fmpq_clear(ar);
        fmpq_clear(bm);
        fmpq_clear(br);
        fmpq_clear(t);
        fmpq_clear(u);
    }

    /* check a simple hardcoded example */
    {
        slong prec;
        arb_t xy, yz, y, v, w;

        prec = 32;

        arb_init(xy);
        arb_init(yz);
        arb_init(y);
        arb_init(v);
        arb_init(w);

        arb_set_str(xy, "1 +/- 1", prec);
        arb_set_str(yz, "2 +/- 1", prec);
        arb_set_str(y, "1.5 +/- 0.6", prec);
        arb_set_str(v, "1.5 +/- 0.4", prec);

        arb_intersection(w, xy, yz, prec);

        if (!arb_contains(y, w) || !arb_contains(w, v))
        {
            flint_printf("FAIL (hardcoded example)\n\n");
            flint_printf("xy = "); arb_print(xy); flint_printf("\n\n");
            flint_printf("yx = "); arb_print(yz); flint_printf("\n\n");
            flint_printf("y = "); arb_print(y); flint_printf("\n\n");
            flint_printf("w = "); arb_print(w); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(xy);
        arb_clear(yz);
        arb_clear(y);
        arb_clear(v);
        arb_clear(w);
    }

    TEST_FUNCTION_END(state);
}
