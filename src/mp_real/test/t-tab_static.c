/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mp_real.h"
#include "mp_real/impl.h"

/* the static tables hold only value limbs: the top
   MP_REAL_STATIC_TAB_N limbs of the dynamic entries built one limb
   deeper (whose own bottom limb is the generation guard and is not
   stored); enforce that identity bit for bit */

TEST_FUNCTION_START(mp_real_tab_static, state)
{
    slong i, j;

    /* go through the exported entry accessors: the thread-local
       table data itself is not DLL-exported on Windows */
    _mp_real_exp_logs_ensure(MP_REAL_STATIC_TAB_N, MP_REAL_STATIC_TAB_R);
    _mp_real_atans_ensure(MP_REAL_STATIC_TAB_N, MP_REAL_STATIC_TAB_R);

    for (i = 0; i <= MP_REAL_STATIC_TAB_R; i++)
    {
        nn_srcptr el = _mp_real_exp_logs_entry(i, MP_REAL_STATIC_TAB_N);
        nn_srcptr ea = _mp_real_atans_entry(i, MP_REAL_STATIC_TAB_N);
        for (j = 0; j < MP_REAL_STATIC_TAB_N; j++)
        {
            if (_mp_real_exp_logs_static[i * MP_REAL_STATIC_TAB_N + j]
                != el[j])
                TEST_FUNCTION_FAIL("logs entry %wd limb %wd\n", i, j);
            if (_mp_real_atans_static[i * MP_REAL_STATIC_TAB_N + j]
                != ea[j])
                TEST_FUNCTION_FAIL("atans entry %wd limb %wd\n", i, j);
        }
    }

    TEST_FUNCTION_END(state);
}
