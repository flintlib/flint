/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"
#include "fixed.h"

/* every dynamic table entry's value limbs (all limbs above the
   guard) must equal the exact floor of the infinite-precision value
   at that precision -- hence at every shorter truncation too, floors
   being nested; the guard limb itself is one-sided but not exact
   (e.g. the multi-summation tier biases it down deliberately) */

TEST_FUNCTION_START(fixed_tab_floor, state)
{
    slong iter;

    for (iter = 0; iter < 2 + flint_test_multiplier() / 3; iter++)
    {
        slong nv = 1 + n_randint(state, 40);
        slong rc = 8 + n_randint(state, 120);
        slong i, nc;
        nn_ptr exact;

        _fixed_exp_logs_ensure(nv, rc);
        _fixed_atans_ensure(nv, rc);

        /* compare at the requested width nv through the exported
           entry accessors (the thread-local data is not
           DLL-exported): exact floors truncate to exact floors, so
           this checks the guarantee at every ensure width */
        nc = nv;
        exact = flint_malloc(nc * sizeof(ulong));

        for (i = 0; i <= rc; i++)
        {
            _fixed_tab_entry_exact(exact, 0, (ulong) i, nc);
            if (mpn_cmp(exact, _fixed_exp_logs_entry(i, nc), nc) != 0)
                TEST_FUNCTION_FAIL("logs: nv = %wd, rc = %wd, "
                    "i = %wd\n", nv, rc, i);
        }

        for (i = 0; i <= rc; i++)
        {
            _fixed_tab_entry_exact(exact, 1, (ulong) i, nc);
            if (mpn_cmp(exact, _fixed_atans_entry(i, nc), nc) != 0)
                TEST_FUNCTION_FAIL("atans: nv = %wd, rc = %wd, "
                    "i = %wd\n", nv, rc, i);
        }
        flint_free(exact);
    }

    TEST_FUNCTION_END(state);
}
