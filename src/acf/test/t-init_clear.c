/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acf.h"

TEST_FUNCTION_START(acf_init_clear, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acf_t x, y;

        acf_init(x);
        acf_init(y);

        arf_randtest(acf_realref(x), state, 200, 100);
        arf_randtest(acf_imagref(x), state, 200, 100);
        acf_set(y, x);

        acf_clear(x);
        acf_clear(y);
    }

    TEST_FUNCTION_END(state);
}
