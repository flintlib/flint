/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2025 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "ulong_extras.h"
#include "mpn_extras.h"
#include "test_helpers.h"

#define DCL static void

DCL nadd2(nn_ptr rp, nn_srcptr ap, nn_srcptr bp) { NN_ADD_2(rp, ap, bp); }
DCL nadd3(nn_ptr rp, nn_srcptr ap, nn_srcptr bp) { NN_ADD_3(rp, ap, bp); }
DCL nadd4(nn_ptr rp, nn_srcptr ap, nn_srcptr bp) { NN_ADD_4(rp, ap, bp); }
DCL nadd5(nn_ptr rp, nn_srcptr ap, nn_srcptr bp) { NN_ADD_5(rp, ap, bp); }
DCL nadd6(nn_ptr rp, nn_srcptr ap, nn_srcptr bp) { NN_ADD_6(rp, ap, bp); }
DCL nadd7(nn_ptr rp, nn_srcptr ap, nn_srcptr bp) { NN_ADD_7(rp, ap, bp); }
DCL nadd8(nn_ptr rp, nn_srcptr ap, nn_srcptr bp) { NN_ADD_8(rp, ap, bp); }
static void (*adds[])(nn_ptr, nn_srcptr, nn_srcptr)=
    {NULL, NULL, nadd2, nadd3, nadd4, nadd5, nadd6, nadd7, nadd8};

TEST_FUNCTION_START(nn_add, state)
{
    int i, j, result;

    for (i = 0; i < 50000 * flint_test_multiplier(); i++)
    {
        ulong _rp[8], sp[8], ap[8], bp[8];
        slong n = 2 + n_randint(state, 7);
        int alias = n_randint(state, 2);
        nn_ptr rp = alias ? ap : _rp;

        for (j = 0; j < 8; j++)
        {
            sp[j] = n_randtest(state);
            ap[j] = n_randtest(state);
            bp[j] = n_randtest(state);
        }

        mpn_add_n(sp, ap, bp, n);
        adds[n](rp, ap, bp);

        result = flint_mpn_equal_p(rp, sp, n);

        if (!result)
            TEST_FUNCTION_FAIL("Aliasing: %d\n"     "n = %d\n"
                               "rp = %{ulong*}\n"   "sp = %{ulong*}\n"
                               "ap = %{ulong*}\n"   "bp = %{ulong*}\n" "i=%d\n",
                               alias, n, rp, n, sp, n, ap, n, bp, n, i);
    }

    TEST_FUNCTION_END(state);
}

#undef DCL
