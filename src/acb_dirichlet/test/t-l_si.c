/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_dirichlet.h"

TEST_FUNCTION_START(acb_dirichlet_l_si, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * 0.1 * flint_test_multiplier(); iter++)
    {
        dirichlet_group_t G;
        dirichlet_char_t chi;
        acb_t a, b, t;
        ulong q;
        slong s, prec;

        q = 1 + n_randint(state, n_randint(state, 4) ? 40 : 300);

        switch (n_randint(state, 4))
        {
            case 0:
                prec = 2 + n_randint(state, 64);
                break;
            case 1:
                prec = 2 + n_randint(state, 300);
                break;
            case 2:
                prec = 2 + n_randint(state, 1000);
                break;
            default:
                prec = 2 + n_randint(state, 3000);
        }

        if (q > 100)
            prec = FLINT_MIN(prec, 500);

        if (n_randint(state, 4) == 0)
            s = (slong) n_randint(state, 200) - 100;
        else
            s = (slong) n_randint(state, 30) - 12;

        dirichlet_group_init(G, q);
        dirichlet_char_init(chi, G);
        dirichlet_char_index(chi, G, n_randint(state, G->phi_q));

        acb_init(a);
        acb_init(b);
        acb_init(t);

        flint_set_num_threads(1 + n_randint(state, 3));

        acb_dirichlet_l_si(a, s, (q == 1 && n_randint(state, 2)) ? NULL : G, chi, prec);

        flint_set_num_threads(1);

        acb_set_si(t, s);
        _acb_dirichlet_l_generic(b, t, G, chi, prec + 20);

        if (!acb_overlaps(a, b) || (acb_is_finite(b) && !acb_is_finite(a)))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("q = %wu  chi = %wu  s = %wd  prec = %wd\n\n", q, chi->n, s, prec);
            flint_printf("a = "); acb_printn(a, 50, 0); flint_printf("\n\n");
            flint_printf("b = "); acb_printn(b, 50, 0); flint_printf("\n\n");
            flint_abort();
        }

        /* the closed forms and the special algorithms should be accurate */
        if (acb_is_finite(b) && !acb_contains_zero(b) && prec >= 16 &&
            acb_rel_accuracy_bits(a) < prec - 10)
        {
            flint_printf("FAIL: accuracy\n\n");
            flint_printf("q = %wu  chi = %wu  s = %wd  prec = %wd\n\n", q, chi->n, s, prec);
            flint_printf("a = "); acb_printn(a, 50, 0); flint_printf("\n\n");
            flint_printf("b = "); acb_printn(b, 50, 0); flint_printf("\n\n");
            flint_abort();
        }

        /* exact zeros */
        if (s <= 0 && dirichlet_char_is_primitive(G, chi) && q > 1 &&
            ((1 - s) % 2) != (slong) dirichlet_parity_char(G, chi) && !acb_is_zero(a))
        {
            flint_printf("FAIL: zero\n\n");
            flint_printf("q = %wu  chi = %wu  s = %wd  prec = %wd\n\n", q, chi->n, s, prec);
            flint_printf("a = "); acb_printn(a, 50, 0); flint_printf("\n\n");
            flint_abort();
        }

        /* acb_dirichlet_l with an exact integer uses acb_dirichlet_l_si */
        acb_dirichlet_l(b, t, G, chi, prec);
        if (!acb_equal(a, b) && !(acb_is_finite(a) == 0 && acb_is_finite(b) == 0))
        {
            /* the thread count can only change the radius, not the value */
            if (!acb_overlaps(a, b))
            {
                flint_printf("FAIL: acb_dirichlet_l\n\n");
                flint_printf("q = %wu  chi = %wu  s = %wd  prec = %wd\n\n", q, chi->n, s, prec);
                flint_abort();
            }
        }

        dirichlet_char_clear(chi);
        dirichlet_group_clear(G);
        acb_clear(a);
        acb_clear(b);
        acb_clear(t);
    }

    /* Borwein's algorithm, for all applicable primitive characters */
    for (iter = 0; iter < 300 * 0.1 * flint_test_multiplier(); iter++)
    {
        dirichlet_group_t G;
        dirichlet_char_t chi;
        acb_t a, b, t;
        ulong q, s;
        slong prec;

        do {
            q = 1 + n_randint(state, 100);
        } while (q % 4 == 2);

        s = 1 + n_randint(state, 12);
        prec = 2 + n_randint(state, 1500);

        dirichlet_group_init(G, q);
        dirichlet_char_init(chi, G);
        dirichlet_char_index(chi, G, n_randint(state, G->phi_q));

        acb_init(a);
        acb_init(b);
        acb_init(t);

        if (dirichlet_char_is_primitive(G, chi) && !(q == 1 && s == 1))
        {
            flint_set_num_threads(1 + n_randint(state, 3));
            acb_dirichlet_l_ui_borwein_bsplit(a, s, (q == 1) ? NULL : G, chi, prec);
            flint_set_num_threads(1);

            acb_set_ui(t, s);
            _acb_dirichlet_l_generic(b, t, G, chi, prec + 20);

            if (!acb_overlaps(a, b))
            {
                flint_printf("FAIL: borwein overlap\n\n");
                flint_printf("q = %wu  chi = %wu  s = %wu  prec = %wd\n\n", q, chi->n, s, prec);
                flint_printf("a = "); acb_printn(a, 50, 0); flint_printf("\n\n");
                flint_printf("b = "); acb_printn(b, 50, 0); flint_printf("\n\n");
                flint_abort();
            }

            /* indeterminate only if 1 - chi(2) 2^(1-s) = 0 */
            if (!(s == 1 && q % 2 == 1 && dirichlet_chi(G, chi, 2) == 0) &&
                (!acb_is_finite(a) || (prec >= 16 && acb_rel_accuracy_bits(a) < prec - 10)))
            {
                flint_printf("FAIL: borwein accuracy\n\n");
                flint_printf("q = %wu  chi = %wu  s = %wu  prec = %wd\n\n", q, chi->n, s, prec);
                flint_printf("a = "); acb_printn(a, 50, 0); flint_printf("\n\n");
                flint_abort();
            }
        }

        dirichlet_char_clear(chi);
        dirichlet_group_clear(G);
        acb_clear(a);
        acb_clear(b);
        acb_clear(t);
    }

    TEST_FUNCTION_END(state);
}
