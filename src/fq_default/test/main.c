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

#include "t-ctx_init_modulus.c"
#include "t-ctx_init_modulus_nmod.c"
#include "t-ctx_modulus.c"
#include "t-get_set_fmpz.c"
#include "t-get_set_fmpz_mod_poly.c"
#include "t-get_set_fmpz_poly.c"
#include "t-init.c"
#include "t-inlines.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(fq_default_ctx_init_modulus),
    TEST_FUNCTION(fq_default_ctx_init_modulus_nmod),
    TEST_FUNCTION(fq_default_ctx_modulus),
    TEST_FUNCTION(fq_default_get_set_fmpz),
    TEST_FUNCTION(fq_default_get_set_fmpz_mod_poly),
    TEST_FUNCTION(fq_default_get_set_fmpz_poly),
    TEST_FUNCTION(fq_default_init),
    TEST_FUNCTION(fq_default_inlines)
};

/* main function *************************************************************/

TEST_MAIN(tests)
