/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"

/* Include functions *********************************************************/

#include "t-series.c"
#include "t-series_acb.c"
#include "t-series_arb.c"
#include "t-series_fmpq.c"
#include "t-series_fmpz.c"
#include "t-series_mod.c"
#include "t-series_nmod8.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(gr_series),
    TEST_FUNCTION(gr_series_acb),
    TEST_FUNCTION(gr_series_arb),
    TEST_FUNCTION(gr_series_fmpq),
    TEST_FUNCTION(gr_series_fmpz),
    TEST_FUNCTION(gr_series_mod),
    TEST_FUNCTION(gr_series_nmod8),
};

/* main function *************************************************************/

TEST_MAIN(tests)
