/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arf.h"

TEST_FUNCTION_START(arf_sgn, state)
{
    slong ix;

    for (ix = 0; ix < 100 * flint_test_multiplier(); ix++)
    {
        arf_t x;
        int sgn;

        arf_init(x);

        arf_randtest_special(x, state, 30, 30);

        sgn = arf_sgn(x);

        switch (sgn)
        {
            case -1:
                if (arf_is_neg_inf(x) || arf_cmp_si(x, 0) < 0)
                    break;
                else
                    goto fail;
            case 0:
                if (arf_is_zero(x) || arf_is_nan(x))
                    break;
                else
                    goto fail;
            case 1:
                if (arf_is_pos_inf(x) || arf_cmp_si(x, 0) > 0)
                    break;
                else
                    goto fail;
            default:
fail:           flint_throw(FLINT_TEST_FAIL,
                        "x = %{arf}\n"
                        "sgn = %d\n",
                        x, sgn);
        }

        arf_clear(x);
    }

    TEST_FUNCTION_END(state);
}
