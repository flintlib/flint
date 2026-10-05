/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "test_helpers.h"
#include "acb.h"
#include "fexpr.h"

static void
check(const acb_t z, slong digits, const char * expected)
{
    fexpr_t x;
    char * s;

    fexpr_init(x);
    fexpr_set_acb_decimal(x, z, digits);
    s = fexpr_get_str(x);

    if (strcmp(s, expected) != 0)
    {
        flint_printf("FAIL\n\n");
        flint_printf("z = "); acb_printd(z, 10); flint_printf("\n");
        flint_printf("got      %s\n", s);
        flint_printf("expected %s\n", expected);
        flint_abort();
    }

    flint_free(s);
    fexpr_clear(x);
}

TEST_FUNCTION_START(fexpr_set_acb_decimal, state)
{
    acb_t z;

    acb_init(z);

    check(z, 10, "0");

    acb_set_d(z, 1.5);
    check(z, 10, "Decimal(\"1.500000000\")");

    acb_set_d_d(z, -0.25, 2.0);
    check(z, 4, "Add(Decimal(\"-0.2500\"), Mul(Decimal(\"2.000\"), NumberI))");

    acb_set_d_d(z, 0.0, -3.0);
    check(z, 3, "Mul(Decimal(\"-3.00\"), NumberI)");

    acb_indeterminate(z);
    check(z, 3, "Undefined");

    arb_pos_inf(acb_realref(z));
    arb_zero(acb_imagref(z));
    check(z, 3, "Infinity");

    acb_clear(z);

    TEST_FUNCTION_END(state);
}
