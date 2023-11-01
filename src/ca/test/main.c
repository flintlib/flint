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

#include "t-acos.c"
#include "t-add.c"
#include "t-asin.c"
#include "t-atan.c"
#include "t-conj.c"
#include "t-ctx_init_clear.c"
#include "t-div.c"
#include "t-erf.c"
#include "t-exp.c"
#include "t-field_init_clear.c"
#include "t-fmpz_mpoly_evaluate.c"
#include "t-gamma.c"
#include "t-get_fexpr.c"
#include "t-get_str.c"
#include "t-init_clear.c"
#include "t-log.c"
#include "t-log_identities.c"
#include "t-mul.c"
#include "t-neg.c"
#include "t-pow.c"
#include "t-pow_si_arithmetic.c"
#include "t-properties.c"
#include "t-re_im.c"
#include "t-sin_cos.c"
#include "t-sqrt.c"
#include "t-sqrt_factor.c"
#include "t-sub.c"
#include "t-tan.c"
#include "t-transfer.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(ca_acos),
    TEST_FUNCTION(ca_add),
    TEST_FUNCTION(ca_asin),
    TEST_FUNCTION(ca_atan),
    TEST_FUNCTION(ca_conj),
    TEST_FUNCTION(ca_ctx_init_clear),
    TEST_FUNCTION(ca_div),
    TEST_FUNCTION(ca_erf),
    TEST_FUNCTION(ca_exp),
    TEST_FUNCTION(ca_field_init_clear),
    TEST_FUNCTION(ca_fmpz_mpoly_evaluate),
    TEST_FUNCTION(ca_gamma),
    TEST_FUNCTION(ca_get_fexpr),
    TEST_FUNCTION(ca_get_str),
    TEST_FUNCTION(ca_init_clear),
    TEST_FUNCTION(ca_log),
    TEST_FUNCTION(ca_log_identities),
    TEST_FUNCTION(ca_mul),
    TEST_FUNCTION(ca_neg),
    TEST_FUNCTION(ca_pow),
    TEST_FUNCTION(ca_pow_si_arithmetic),
    TEST_FUNCTION(ca_properties),
    TEST_FUNCTION(ca_re_im),
    TEST_FUNCTION(ca_sin_cos),
    TEST_FUNCTION(ca_sqrt),
    TEST_FUNCTION(ca_sqrt_factor),
    TEST_FUNCTION(ca_sub),
    TEST_FUNCTION(ca_tan),
    TEST_FUNCTION(ca_transfer)
};

/* main function *************************************************************/

TEST_MAIN(tests)
