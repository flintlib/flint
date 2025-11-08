/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Include functions *********************************************************/

#include "t-init.c"
#include "t-special_divrems.c"
#include "t-dft.c"
#include "t-idft.c"
#include "t-dft_t.c"
#include "t-idft_t.c"
#include "t-tft.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(n_fft_ctx_init2),
    TEST_FUNCTION(n_fft_special_divrems),
    TEST_FUNCTION(n_fft_dft),
    TEST_FUNCTION(n_fft_idft),
    TEST_FUNCTION(n_fft_dft_t),
    TEST_FUNCTION(n_fft_idft_t),
    TEST_FUNCTION(n_fft_tft),
    /* TEST_FUNCTION(n_fft_itft), */
    /* TEST_FUNCTION(n_fft_tft_t), */
    /* TEST_FUNCTION(n_fft_itft_t), */
};

/* main function *************************************************************/

TEST_MAIN(tests)
