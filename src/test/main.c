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

#include "t-add_ssaaaa.c"
#include "t-add_sssaaaaaa.c"
#include "t-add_ssssaaaaaaaa.c"
#include "t-byte_swap.c"
#include "t-flint_clz.c"
#include "t-flint_ctz.c"
#include "t-io.c"
#include "t-sdiv_qrnnd.c"
#include "t-smul_ppmm.c"
#include "t-sub_dddmmmsss.c"
#include "t-sub_ddmmss.c"
#include "t-udiv_qrnnd.c"
#include "t-udiv_qrnnd_preinv.c"
#include "t-umul_ppmm.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(add_ssaaaa),
    TEST_FUNCTION(add_sssaaaaaa),
    TEST_FUNCTION(add_ssssaaaaaaaa),
    TEST_FUNCTION(byte_swap),
    TEST_FUNCTION(flint_clz),
    TEST_FUNCTION(flint_ctz),
    TEST_FUNCTION(flint_fprintf),
    TEST_FUNCTION(flint_printf),
    TEST_FUNCTION(sdiv_qrnnd),
    TEST_FUNCTION(smul_ppmm),
    TEST_FUNCTION(sub_dddmmmsss),
    TEST_FUNCTION(sub_ddmmss),
    TEST_FUNCTION(udiv_qrnnd),
    TEST_FUNCTION(udiv_qrnnd_preinv),
    TEST_FUNCTION(umul_ppmm)
};

/* main function *************************************************************/

TEST_MAIN(tests)
