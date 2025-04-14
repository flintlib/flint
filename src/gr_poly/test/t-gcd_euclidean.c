/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "gr_poly.h"

FLINT_DLL extern gr_static_method_table _ca_methods;

TEST_FUNCTION_START(gr_poly_gcd_euclidean, state)
{
    _gr_poly_test_gcd_field((gr_method_poly_gcd_op) _gr_poly_gcd_euclidean, state, 1000, 6, NULL);

    TEST_FUNCTION_END(state);
}
