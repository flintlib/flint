/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Include functions *********************************************************/

#include "t-add_sub_neg.c"
#include "t-discrete_log_pohlig_hellman.c"
#include "t-dot_nlimbs.c"
#include "t-dot.c"
#include "t-dot_ptr.c"
#include "t-nmod.c"
#include "t-nmod_pow_fmpz.c"
#include "t-reduce.c"
#include "t-scalar_addmul_nmod.c"
#include "t-scalar_mul_nmod.c"
#include "t-scalar_mul_nmod_shoup.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(nmod_vec_add_sub_neg),
    TEST_FUNCTION(nmod_vec_discrete_log_pohlig_hellman),
    TEST_FUNCTION(_nmod_vec_dot_params),
    TEST_FUNCTION(nmod_vec_dot),
    TEST_FUNCTION(nmod_vec_dot_ptr),
    TEST_FUNCTION(nmod_vec_nmod),
    TEST_FUNCTION(nmod_vec_nmod_pow_fmpz),
    TEST_FUNCTION(nmod_vec_reduce),
    TEST_FUNCTION(nmod_vec_scalar_addmul_nmod),
    TEST_FUNCTION(nmod_vec_scalar_mul_nmod),
    TEST_FUNCTION(nmod_vec_scalar_mul_nmod_shoup)
};

/* main function *************************************************************/

TEST_MAIN(tests)
