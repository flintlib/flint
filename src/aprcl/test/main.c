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

#include "t-config_gauss.c"
#include "t-config_jacobi.c"
#include "t-f_table.c"
#include "t-is_prime.c"
#include "t-is_prime_gauss.c"
#include "t-is_prime_jacobi.c"
#include "t-unity_zp_add.c"
#include "t-unity_zp_aut_inv.c"
#include "t-unity_zp_equal.c"
#include "t-unity_zp_init.c"
#include "t-unity_zp_is_unity.c"
#include "t-unity_zp_jacobi_sum.c"
#include "t-unity_zp_mul11.c"
#include "t-unity_zp_mul2.c"
#include "t-unity_zp_mul3.c"
#include "t-unity_zp_mul5.c"
#include "t-unity_zp_mul7.c"
#include "t-unity_zp_mul.c"
#include "t-unity_zp_pow_2k.c"
#include "t-unity_zp_pow.c"
#include "t-unity_zp_pow_sliding.c"
#include "t-unity_zpq_add.c"
#include "t-unity_zpq_equal.c"
#include "t-unity_zpq_gauss_sum.c"
#include "t-unity_zpq_init.c"
#include "t-unity_zpq_mul.c"
#include "t-unity_zpq_mul_unity_p.c"
#include "t-unity_zpq_pow.c"
#include "t-unity_zp_reduce_cyclotomic.c"
#include "t-unity_zp_sqr11.c"
#include "t-unity_zp_sqr2.c"
#include "t-unity_zp_sqr3.c"
#include "t-unity_zp_sqr5.c"
#include "t-unity_zp_sqr7.c"
#include "t-unity_zp_sqr.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(aprcl_config_gauss),
    TEST_FUNCTION(aprcl_config_jacobi),
    TEST_FUNCTION(aprcl_f_table),
    TEST_FUNCTION(aprcl_is_prime),
    TEST_FUNCTION(aprcl_is_prime_gauss),
    TEST_FUNCTION(aprcl_is_prime_jacobi),
    TEST_FUNCTION(aprcl_unity_zp_add),
    TEST_FUNCTION(aprcl_unity_zp_aut_inv),
    TEST_FUNCTION(aprcl_unity_zp_equal),
    TEST_FUNCTION(aprcl_unity_zp_init),
    TEST_FUNCTION(aprcl_unity_zp_is_unity),
    TEST_FUNCTION(aprcl_unity_zp_jacobi_sum),
    TEST_FUNCTION(aprcl_unity_zp_mul11),
    TEST_FUNCTION(aprcl_unity_zp_mul2),
    TEST_FUNCTION(aprcl_unity_zp_mul3),
    TEST_FUNCTION(aprcl_unity_zp_mul5),
    TEST_FUNCTION(aprcl_unity_zp_mul7),
    TEST_FUNCTION(aprcl_unity_zp_mul),
    TEST_FUNCTION(aprcl_unity_zp_pow_2k),
    TEST_FUNCTION(aprcl_unity_zp_pow),
    TEST_FUNCTION(aprcl_unity_zp_pow_sliding),
    TEST_FUNCTION(aprcl_unity_zpq_add),
    TEST_FUNCTION(aprcl_unity_zpq_equal),
    TEST_FUNCTION(aprcl_unity_zpq_gauss_sum),
    TEST_FUNCTION(aprcl_unity_zpq_init),
    TEST_FUNCTION(aprcl_unity_zpq_mul),
    TEST_FUNCTION(aprcl_unity_zpq_mul_unity_p),
    TEST_FUNCTION(aprcl_unity_zpq_pow),
    TEST_FUNCTION(aprcl_unity_zp_reduce_cyclotomic),
    TEST_FUNCTION(aprcl_unity_zp_sqr11),
    TEST_FUNCTION(aprcl_unity_zp_sqr2),
    TEST_FUNCTION(aprcl_unity_zp_sqr3),
    TEST_FUNCTION(aprcl_unity_zp_sqr5),
    TEST_FUNCTION(aprcl_unity_zp_sqr7),
    TEST_FUNCTION(aprcl_unity_zp_sqr)
};

/* main function *************************************************************/

TEST_MAIN(tests)
