/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <stdlib.h>

/* Include functions *********************************************************/

#include "t-add_sub_neg.c"
#include "t-discrete_log_pohlig_hellman.c"
#include "t-divides.c"
#include "t-inv.c"
#include "t-mul.c"
#include "t-next_smooth_prime.c"
#include "t-pow_fmpz.c"
#include "t-pow_ui.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(fmpz_mod_add_sub_neg),
    TEST_FUNCTION(fmpz_mod_discrete_log_pohlig_hellman),
    TEST_FUNCTION(fmpz_mod_divides),
    TEST_FUNCTION(fmpz_mod_inv),
    TEST_FUNCTION(fmpz_mod_mul),
    TEST_FUNCTION(fmpz_mod_next_smooth_prime),
    TEST_FUNCTION(fmpz_mod_pow_fmpz),
    TEST_FUNCTION(fmpz_mod_pow_ui)
};

/* main function *************************************************************/

TEST_MAIN(tests)
