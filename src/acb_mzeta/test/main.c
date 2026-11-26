/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Include functions *********************************************************/

#include "t-epstein.c"
#include "t-epstein_lattice.c"
#include "t-epstein_shell.c"
#include "t-epstein_theta.c"
#include "t-g_suffixes.c"
#include "t-gamma_upper_scaled.c"
#include "t-madelung_nacl.c"
#include "t-mzv.c"
#include "t-polylog.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(acb_mzeta_epstein),
    TEST_FUNCTION(acb_mzeta_epstein_lattice),
    TEST_FUNCTION(acb_mzeta_epstein_shell),
    TEST_FUNCTION(acb_mzeta_epstein_theta),
    TEST_FUNCTION(acb_mzeta_G_suffixes),
    TEST_FUNCTION(acb_mzeta_gamma_upper_scaled),
    TEST_FUNCTION(acb_mzeta_madelung_nacl),
    TEST_FUNCTION(acb_mzeta_mzv),
    TEST_FUNCTION(acb_mzeta_polylog),
};

/* main function *************************************************************/

TEST_MAIN(tests)
