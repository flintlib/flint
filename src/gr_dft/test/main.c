/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Include functions *********************************************************/

#include "t-dft.c"
#include "t-dft_inverse.c"
#include "t-acb_dft.c"
#include "t-dirichlet_dft.c"
#include "t-nfixed.c"
#include "t-prod.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(gr_dft),
    TEST_FUNCTION(gr_dft_inverse),
    TEST_FUNCTION(gr_dft_acb),
    TEST_FUNCTION(gr_dft_dirichlet),
    TEST_FUNCTION(gr_dft_nfixed),
    TEST_FUNCTION(gr_dft_prod),
};

/* main function *************************************************************/

TEST_MAIN(tests)
