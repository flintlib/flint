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

#include "t-e_inc.c"
#include "t-f.c"
#include "t-invariants.c"
#include "t-inv_p.c"
#include "t-pi.c"
#include "t-pi_inc.c"
#include "t-p_p_prime.c"
#include "t-rc1.c"
#include "t-rf.c"
#include "t-rg.c"
#include "t-rj.c"
#include "t-sigma.c"
#include "t-zeta.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(acb_elliptic_e_inc),
    TEST_FUNCTION(acb_elliptic_f),
    TEST_FUNCTION(acb_elliptic_invariants),
    TEST_FUNCTION(acb_elliptic_inv_p),
    TEST_FUNCTION(acb_elliptic_pi),
    TEST_FUNCTION(acb_elliptic_pi_inc),
    TEST_FUNCTION(acb_elliptic_p_p_prime),
    TEST_FUNCTION(acb_elliptic_rc1),
    TEST_FUNCTION(acb_elliptic_rf),
    TEST_FUNCTION(acb_elliptic_rg),
    TEST_FUNCTION(acb_elliptic_rj),
    TEST_FUNCTION(acb_elliptic_sigma),
    TEST_FUNCTION(acb_elliptic_zeta)
};

/* main function *************************************************************/

TEST_MAIN(tests)
