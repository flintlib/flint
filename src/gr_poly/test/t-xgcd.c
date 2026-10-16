/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "gr_poly.h"

TEST_FUNCTION_START(gr_poly_xgcd, state)
{
    _gr_poly_test_xgcd((gr_method_poly_xgcd_op) _gr_poly_xgcd,
        state,
        100 * flint_test_multiplier(),
        0,  /* default maxn */
        1,  /* expected to succeed over appropriate UFDs (not just fields) */
        NULL /* test random contexts */
    );

    TEST_FUNCTION_END(state);
}
