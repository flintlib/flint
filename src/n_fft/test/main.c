/*
    Copyright (C) 2024 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Include functions *********************************************************/

#include "t-init.c"
#include "t-dft.c"
#include "t-idft.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(n_fft_ctx_init2),
    TEST_FUNCTION(n_fft_dft),
    TEST_FUNCTION(n_fft_idft),
};

/* main function *************************************************************/

TEST_MAIN(tests)
