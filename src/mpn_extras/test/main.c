/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <stdlib.h>

/* Include functions *********************************************************/

#include "t-2add_n_inplace.c"
#include "t-divides.c"
#include "t-divrem_preinv1.c"
#include "t-divrem_preinvn.c"
#include "t-fmms1.c"
#include "t-gcd_full.c"
#include "t-mod_preinvn.c"
#include "t-mul.c"
#include "t-mul_n.c"
#include "t-mul_toom22.c"
#include "t-mullow_n.c"
#include "t-mulhigh_n.c"
#include "t-mulhigh_normalised.c"
#include "t-mulmod_2expp1.c"
#include "t-mulmod_preinv1.c"
#include "t-mulmod_preinvn.c"
#include "t-remove_2exp.c"
#include "t-remove_power.c"
#include "t-sqr.c"
#include "t-sqrhigh.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(flint_mpn_2add_n_inplace),
    TEST_FUNCTION(flint_mpn_divides),
    TEST_FUNCTION(flint_mpn_divrem_preinv1),
    TEST_FUNCTION(flint_mpn_divrem_preinvn),
    TEST_FUNCTION(flint_mpn_fmms1),
    TEST_FUNCTION(flint_mpn_gcd_full),
    TEST_FUNCTION(flint_mpn_mod_preinvn),
    TEST_FUNCTION(flint_mpn_mul),
    TEST_FUNCTION(flint_mpn_mul_n),
    TEST_FUNCTION(flint_mpn_mul_toom22),
    TEST_FUNCTION(flint_mpn_mullow_n),
    TEST_FUNCTION(flint_mpn_mulhigh_n),
    TEST_FUNCTION(flint_mpn_mulhigh_normalised),
    TEST_FUNCTION(flint_mpn_mulmod_2expp1),
    TEST_FUNCTION(flint_mpn_mulmod_preinv1),
    TEST_FUNCTION(flint_mpn_mulmod_preinvn),
    TEST_FUNCTION(flint_mpn_remove_2exp),
    TEST_FUNCTION(flint_mpn_remove_power),
    TEST_FUNCTION(flint_mpn_sqr),
    TEST_FUNCTION(flint_mpn_sqrhigh)
};

/* main function *************************************************************/

TEST_MAIN(tests)
