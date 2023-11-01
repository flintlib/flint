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

#include "t-bell_number.c"
#include "t-bell_number_multi_mod.c"
#include "t-bell_number_nmod.c"
#include "t-bell_number_nmod_vec.c"
#include "t-bell_number_vec.c"
#include "t-bernoulli_number.c"
#include "t-bernoulli_number_denom.c"
#include "t-bernoulli_number_vec.c"
#include "t-bernoulli_polynomial.c"
#include "t-chebyshev_t_polynomial.c"
#include "t-chebyshev_u_polynomial.c"
#include "t-divisors.c"
#include "t-euler_number_vec.c"
#include "t-euler_number_zeta.c"
#include "t-euler_polynomial.c"
#include "t-harmonic_number.c"
#include "t-landau_function_vec.c"
#include "t-number_of_partitions_vec.c"
#include "t-ramanujan_tau.c"
#include "t-stirling.c"
#include "t-sum_of_squares.c"
#include "t-swinnerton_dyer_polynomial.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(arith_bell_number),
    TEST_FUNCTION(arith_bell_number_multi_mod),
    TEST_FUNCTION(arith_bell_number_nmod),
    TEST_FUNCTION(arith_bell_number_nmod_vec),
    TEST_FUNCTION(arith_bell_number_vec),
    TEST_FUNCTION(arith_bernoulli_number),
    TEST_FUNCTION(arith_bernoulli_number_denom),
    TEST_FUNCTION(arith_bernoulli_number_vec),
    TEST_FUNCTION(arith_bernoulli_polynomial),
    TEST_FUNCTION(arith_chebyshev_t_polynomial),
    TEST_FUNCTION(arith_chebyshev_u_polynomial),
    TEST_FUNCTION(arith_divisors),
    TEST_FUNCTION(arith_euler_number_vec),
    TEST_FUNCTION(arith_euler_number_zeta),
    TEST_FUNCTION(arith_euler_polynomial),
    TEST_FUNCTION(arith_harmonic_number),
    TEST_FUNCTION(arith_landau_function_vec),
    TEST_FUNCTION(arith_number_of_partitions_vec),
    TEST_FUNCTION(arith_ramanujan_tau),
    TEST_FUNCTION(arith_stirling),
    TEST_FUNCTION(arith_sum_of_squares),
    TEST_FUNCTION(arith_swinnerton_dyer_polynomial)
};

/* main function *************************************************************/

TEST_MAIN(tests)
