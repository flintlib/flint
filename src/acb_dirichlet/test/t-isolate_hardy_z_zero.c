/*
    Copyright (C) 2019 D.H.J. Polymath

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_dirichlet.h"

static void
_check_interval(const arf_t a, const arf_t b, const fmpz_t n)
{
    arb_t v;
    int sa, sb;
    slong prec = arf_bits(b) + 8;

    arb_init(v);

    if (arf_cmp(a, b) >= 0)
    {
        flint_printf("FAIL: interval endpoint ordering\n\n");
        flint_printf("n = "); fmpz_print(n); flint_printf("\n\n");
        flint_printf("a = "); arf_print(a); flint_printf("\n\n");
        flint_printf("b = "); arf_print(b); flint_printf("\n\n");
        flint_abort();
    }

    sa = _acb_dirichlet_definite_hardy_z(v, a, &prec);
    sb = _acb_dirichlet_definite_hardy_z(v, b, &prec);

    if (sa == sb)
    {
        flint_printf("FAIL: interval endpoint function signs\n\n");
        flint_printf("n = "); fmpz_print(n); flint_printf("\n\n");
        flint_printf("a = "); arf_print(a); flint_printf("\n\n");
        flint_printf("b = "); arf_print(b); flint_printf("\n\n");
        flint_printf("sa = %d\n\n", sa);
        flint_printf("sb = %d\n\n", sb);
        flint_abort();
    }

    arb_clear(v);
}

TEST_FUNCTION_START(acb_dirichlet_isolate_hardy_z_zero, state)
{
    slong iter;

    for (iter = 0; iter < 126 + 20 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t a, b;
        fmpz_t n;

        arf_init(a);
        arf_init(b);
        fmpz_init(n);

        if (iter < 126)
        {
            fmpz_set_si(n, iter + 1);
        }
        else
        {
            fmpz_randtest_unsigned(n, state, 20);
            fmpz_add_ui(n, n, 127);
        }

        acb_dirichlet_isolate_hardy_z_zero(a, b, n);
        _check_interval(a, b, n);
        if (fmpz_cmp_si(n, 126) <= 0)
        {
            _acb_dirichlet_isolate_gram_hardy_z_zero(a, b, n);
            _check_interval(a, b, n);
        }
        if (fmpz_cmp_si(n, 13999526) <= 0)
        {
            _acb_dirichlet_isolate_rosser_hardy_z_zero(a, b, n);
            _check_interval(a, b, n);
        }
        if (fmpz_cmp_si(n, 2) >= 0)
        {
            _acb_dirichlet_isolate_turing_hardy_z_zero(a, b, n);
            _check_interval(a, b, n);
        }

        arf_clear(a);
        arf_clear(b);
        fmpz_clear(n);
    }

    TEST_FUNCTION_END(state);
}
