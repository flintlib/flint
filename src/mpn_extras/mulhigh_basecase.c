/*
    Copyright (C) 2024 Albin Ahlbäck
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

#if FLINT_HAVE_ASSEMBLY_x86_64_adx
mp_limb_t flint_mpn_mulhigh_1(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mulhigh_2(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mulhigh_3(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mulhigh_4(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mulhigh_5(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mulhigh_6(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mulhigh_7(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mulhigh_8(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mulhigh_9(mp_ptr, mp_srcptr, mp_srcptr);

mp_limb_pair_t flint_mpn_mulhigh_normalised_1(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_pair_t flint_mpn_mulhigh_normalised_2(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_pair_t flint_mpn_mulhigh_normalised_3(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_pair_t flint_mpn_mulhigh_normalised_4(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_pair_t flint_mpn_mulhigh_normalised_5(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_pair_t flint_mpn_mulhigh_normalised_6(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_pair_t flint_mpn_mulhigh_normalised_7(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_pair_t flint_mpn_mulhigh_normalised_8(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_pair_t flint_mpn_mulhigh_normalised_9(mp_ptr, mp_srcptr, mp_srcptr);

mp_limb_t flint_mpn_sqrhigh_1(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqrhigh_2(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqrhigh_3(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqrhigh_4(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqrhigh_5(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqrhigh_6(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqrhigh_7(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqrhigh_8(mp_ptr, mp_srcptr);

const flint_mpn_mul_func_t flint_mpn_mulhigh_func_tab[] =
{
    NULL,
    flint_mpn_mulhigh_1,
    flint_mpn_mulhigh_2,
    flint_mpn_mulhigh_3,
    flint_mpn_mulhigh_4,
    flint_mpn_mulhigh_5,
    flint_mpn_mulhigh_6,
    flint_mpn_mulhigh_7,
    flint_mpn_mulhigh_8,
    flint_mpn_mulhigh_9
};

const flint_mpn_mulhigh_normalised_func_t flint_mpn_mulhigh_normalised_func_tab[] =
{
    NULL,
    flint_mpn_mulhigh_normalised_1,
    flint_mpn_mulhigh_normalised_2,
    flint_mpn_mulhigh_normalised_3,
    flint_mpn_mulhigh_normalised_4,
    flint_mpn_mulhigh_normalised_5,
    flint_mpn_mulhigh_normalised_6,
    flint_mpn_mulhigh_normalised_7,
    flint_mpn_mulhigh_normalised_8,
    flint_mpn_mulhigh_normalised_9
};

const flint_mpn_sqr_func_t flint_mpn_sqrhigh_func_tab[] =
{
    NULL,
    flint_mpn_sqrhigh_1,
    flint_mpn_sqrhigh_2,
    flint_mpn_sqrhigh_3,
    flint_mpn_sqrhigh_4,
    flint_mpn_sqrhigh_5,
    flint_mpn_sqrhigh_6,
    flint_mpn_sqrhigh_7,
    flint_mpn_sqrhigh_8
};
#elif FLINT_HAVE_ASSEMBLY_armv8
mp_limb_t flint_mpn_mulhigh_1(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mulhigh_2(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mulhigh_3(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mulhigh_4(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mulhigh_5(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mulhigh_6(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mulhigh_7(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mulhigh_8(mp_ptr, mp_srcptr, mp_srcptr);

mp_limb_t flint_mpn_sqrhigh_1(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqrhigh_2(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqrhigh_3(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqrhigh_4(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqrhigh_5(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqrhigh_6(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqrhigh_7(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqrhigh_8(mp_ptr, mp_srcptr);

const flint_mpn_mul_func_t flint_mpn_mulhigh_func_tab[] =
{
    NULL,
    flint_mpn_mulhigh_1,
    flint_mpn_mulhigh_2,
    flint_mpn_mulhigh_3,
    flint_mpn_mulhigh_4,
    flint_mpn_mulhigh_5,
    flint_mpn_mulhigh_6,
    flint_mpn_mulhigh_7,
    flint_mpn_mulhigh_8,
};

const flint_mpn_mulhigh_normalised_func_t flint_mpn_mulhigh_normalised_func_tab[] =
{
    NULL,
};

const flint_mpn_sqr_func_t flint_mpn_sqrhigh_func_tab[] =
{
    NULL,
    flint_mpn_sqrhigh_1,
    flint_mpn_sqrhigh_2,
    flint_mpn_sqrhigh_3,
    flint_mpn_sqrhigh_4,
    flint_mpn_sqrhigh_5,
    flint_mpn_sqrhigh_6,
    flint_mpn_sqrhigh_7,
    flint_mpn_sqrhigh_8,
};
#else

/* todo: add MPFR-like basecase for use in mulders */
/* todo: squaring code */
/* todo: define the generic basecase also on x86_64_adx,
   and use to test the assembly versions */

mp_limb_t _flint_mpn_mulhigh_basecase(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n)
{
    mp_limb_t b, a, low, t[2];
    slong i;

    FLINT_ASSERT(n >= 3);

    NN_DOTREV_S3_1X1_HIGH(b, a, u, v, n - 1);
    NN_DOTREV_S3_A3_1X1(b, a, low, 0, b, a, u, v, n);
    t[0] = a;
    t[1] = b;

    umul_ppmm(res[1], res[0], u[n - 1], v[1]);
    for (i = 2; i < n; i++)
        res[i] = mpn_addmul_1(res, u + n - i, i, v[i]);

    mpn_add(res, res, n, t, 2);
    return low;
}

mp_limb_t flint_mpn_mulhigh_1(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t low;
    umul_ppmm(res[0], low, u[0], v[0]);
    return low;
}

mp_limb_t flint_mpn_mulhigh_2(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, low;
    FLINT_MPN_MUL_2X2(res[1], res[0], low, b, u[1], u[0], v[1], v[0]);
    return low;
}

mp_limb_t flint_mpn_mulhigh_3(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a, low;
    NN_DOTREV_S3_1X1_HIGH(b, a, u, v, 2);
    NN_DOTREV_S3_A3_1X1(b, a, low, 0, b, a, u, v, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[0], 0, b, a, u + 1, v + 1, 2);
    NN_ADDMUL_S2_A2_1X1(res[2], res[1], b, a, u[2], v[2]);
    return low;
}

mp_limb_t flint_mpn_mulhigh_4(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a, low;
    NN_DOTREV_S3_1X1_HIGH(b, a, u, v, 3);
    NN_DOTREV_S3_A3_1X1(b, a, low, 0, b, a, u, v, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[0], 0, b, a, u + 1, v + 1, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, b, a, u + 2, v + 2, 2);
    NN_ADDMUL_S2_A2_1X1(res[3], res[2], b, a, u[3], v[3]);
    return low;
}

mp_limb_t flint_mpn_mulhigh_5(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a, low;
    NN_DOTREV_S3_1X1_HIGH(b, a, u, v, 4);
    NN_DOTREV_S3_A3_1X1(b, a, low, 0, b, a, u, v, 5);
    NN_DOTREV_S3_A3_1X1(b, a, res[0], 0, b, a, u + 1, v + 1, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, b, a, u + 2, v + 2, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u + 3, v + 3, 2);
    NN_ADDMUL_S2_A2_1X1(res[4], res[3], b, a, u[4], v[4]);
    return low;
}

mp_limb_t flint_mpn_mulhigh_6(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a, low;
    NN_DOTREV_S3_1X1_HIGH(b, a, u, v, 5);
    NN_DOTREV_S3_A3_1X1(b, a, low, 0, b, a, u, v, 6);
    NN_DOTREV_S3_A3_1X1(b, a, res[0], 0, b, a, u + 1, v + 1, 5);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, b, a, u + 2, v + 2, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u + 3, v + 3, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[3], 0, b, a, u + 4, v + 4, 2);
    NN_ADDMUL_S2_A2_1X1(res[5], res[4], b, a, u[5], v[5]);
    return low;
}

mp_limb_t flint_mpn_mulhigh_7(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a, low;
    NN_DOTREV_S3_1X1_HIGH(b, a, u, v, 6);
    NN_DOTREV_S3_A3_1X1(b, a, low, 0, b, a, u, v, 7);
    NN_DOTREV_S3_A3_1X1(b, a, res[0], 0, b, a, u + 1, v + 1, 6);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, b, a, u + 2, v + 2, 5);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u + 3, v + 3, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[3], 0, b, a, u + 4, v + 4, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[4], 0, b, a, u + 5, v + 5, 2);
    NN_ADDMUL_S2_A2_1X1(res[6], res[5], b, a, u[6], v[6]);
    return low;
}

mp_limb_t flint_mpn_mulhigh_8(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a, low;
    NN_DOTREV_S3_1X1_HIGH(b, a, u, v, 7);
    NN_DOTREV_S3_A3_1X1(b, a, low, 0, b, a, u, v, 8);
    NN_DOTREV_S3_A3_1X1(b, a, res[0], 0, b, a, u + 1, v + 1, 7);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, b, a, u + 2, v + 2, 6);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u + 3, v + 3, 5);
    NN_DOTREV_S3_A3_1X1(b, a, res[3], 0, b, a, u + 4, v + 4, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[4], 0, b, a, u + 5, v + 5, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[5], 0, b, a, u + 6, v + 6, 2);
    NN_ADDMUL_S2_A2_1X1(res[7], res[6], b, a, u[7], v[7]);
    return low;
}

mp_limb_t flint_mpn_mulhigh_9(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a, low;
    NN_DOTREV_S3_1X1_HIGH(b, a, u, v, 8);
    NN_DOTREV_S3_A3_1X1(b, a, low, 0, b, a, u, v, 9);
    NN_DOTREV_S3_A3_1X1(b, a, res[0], 0, b, a, u + 1, v + 1, 8);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, b, a, u + 2, v + 2, 7);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u + 3, v + 3, 6);
    NN_DOTREV_S3_A3_1X1(b, a, res[3], 0, b, a, u + 4, v + 4, 5);
    NN_DOTREV_S3_A3_1X1(b, a, res[4], 0, b, a, u + 5, v + 5, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[5], 0, b, a, u + 6, v + 6, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[6], 0, b, a, u + 7, v + 7, 2);
    NN_ADDMUL_S2_A2_1X1(res[8], res[7], b, a, u[8], v[8]);
    return low;
}

mp_limb_t flint_mpn_mulhigh_10(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a, low;
    NN_DOTREV_S3_1X1_HIGH(b, a, u, v, 9);
    NN_DOTREV_S3_A3_1X1(b, a, low, 0, b, a, u, v, 10);
    NN_DOTREV_S3_A3_1X1(b, a, res[0], 0, b, a, u + 1, v + 1, 9);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, b, a, u + 2, v + 2, 8);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u + 3, v + 3, 7);
    NN_DOTREV_S3_A3_1X1(b, a, res[3], 0, b, a, u + 4, v + 4, 6);
    NN_DOTREV_S3_A3_1X1(b, a, res[4], 0, b, a, u + 5, v + 5, 5);
    NN_DOTREV_S3_A3_1X1(b, a, res[5], 0, b, a, u + 6, v + 6, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[6], 0, b, a, u + 7, v + 7, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[7], 0, b, a, u + 8, v + 8, 2);
    NN_ADDMUL_S2_A2_1X1(res[9], res[8], b, a, u[9], v[9]);
    return low;
}

mp_limb_t flint_mpn_mulhigh_11(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a, low;
    NN_DOTREV_S3_1X1_HIGH(b, a, u, v, 10);
    NN_DOTREV_S3_A3_1X1(b, a, low, 0, b, a, u, v, 11);
    NN_DOTREV_S3_A3_1X1(b, a, res[0], 0, b, a, u + 1, v + 1, 10);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, b, a, u + 2, v + 2, 9);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u + 3, v + 3, 8);
    NN_DOTREV_S3_A3_1X1(b, a, res[3], 0, b, a, u + 4, v + 4, 7);
    NN_DOTREV_S3_A3_1X1(b, a, res[4], 0, b, a, u + 5, v + 5, 6);
    NN_DOTREV_S3_A3_1X1(b, a, res[5], 0, b, a, u + 6, v + 6, 5);
    NN_DOTREV_S3_A3_1X1(b, a, res[6], 0, b, a, u + 7, v + 7, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[7], 0, b, a, u + 8, v + 8, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[8], 0, b, a, u + 9, v + 9, 2);
    NN_ADDMUL_S2_A2_1X1(res[10], res[9], b, a, u[10], v[10]);
    return low;
}

mp_limb_t flint_mpn_mulhigh_12(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a, low;
    NN_DOTREV_S3_1X1_HIGH(b, a, u, v, 11);
    NN_DOTREV_S3_A3_1X1(b, a, low, 0, b, a, u, v, 12);
    NN_DOTREV_S3_A3_1X1(b, a, res[0], 0, b, a, u + 1, v + 1, 11);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, b, a, u + 2, v + 2, 10);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u + 3, v + 3, 9);
    NN_DOTREV_S3_A3_1X1(b, a, res[3], 0, b, a, u + 4, v + 4, 8);
    NN_DOTREV_S3_A3_1X1(b, a, res[4], 0, b, a, u + 5, v + 5, 7);
    NN_DOTREV_S3_A3_1X1(b, a, res[5], 0, b, a, u + 6, v + 6, 6);
    NN_DOTREV_S3_A3_1X1(b, a, res[6], 0, b, a, u + 7, v + 7, 5);
    NN_DOTREV_S3_A3_1X1(b, a, res[7], 0, b, a, u + 8, v + 8, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[8], 0, b, a, u + 9, v + 9, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[9], 0, b, a, u + 10, v + 10, 2);
    NN_ADDMUL_S2_A2_1X1(res[11], res[10], b, a, u[11], v[11]);
    return low;
}

mp_limb_t flint_mpn_mulhigh_13(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a, low;
    NN_DOTREV_S3_1X1_HIGH(b, a, u, v, 12);
    NN_DOTREV_S3_A3_1X1(b, a, low, 0, b, a, u, v, 13);
    NN_DOTREV_S3_A3_1X1(b, a, res[0], 0, b, a, u + 1, v + 1, 12);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, b, a, u + 2, v + 2, 11);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u + 3, v + 3, 10);
    NN_DOTREV_S3_A3_1X1(b, a, res[3], 0, b, a, u + 4, v + 4, 9);
    NN_DOTREV_S3_A3_1X1(b, a, res[4], 0, b, a, u + 5, v + 5, 8);
    NN_DOTREV_S3_A3_1X1(b, a, res[5], 0, b, a, u + 6, v + 6, 7);
    NN_DOTREV_S3_A3_1X1(b, a, res[6], 0, b, a, u + 7, v + 7, 6);
    NN_DOTREV_S3_A3_1X1(b, a, res[7], 0, b, a, u + 8, v + 8, 5);
    NN_DOTREV_S3_A3_1X1(b, a, res[8], 0, b, a, u + 9, v + 9, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[9], 0, b, a, u + 10, v + 10, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[10], 0, b, a, u + 11, v + 11, 2);
    NN_ADDMUL_S2_A2_1X1(res[12], res[11], b, a, u[12], v[12]);
    return low;
}

mp_limb_t flint_mpn_mulhigh_14(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a, low;
    NN_DOTREV_S3_1X1_HIGH(b, a, u, v, 13);
    NN_DOTREV_S3_A3_1X1(b, a, low, 0, b, a, u, v, 14);
    NN_DOTREV_S3_A3_1X1(b, a, res[0], 0, b, a, u + 1, v + 1, 13);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, b, a, u + 2, v + 2, 12);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u + 3, v + 3, 11);
    NN_DOTREV_S3_A3_1X1(b, a, res[3], 0, b, a, u + 4, v + 4, 10);
    NN_DOTREV_S3_A3_1X1(b, a, res[4], 0, b, a, u + 5, v + 5, 9);
    NN_DOTREV_S3_A3_1X1(b, a, res[5], 0, b, a, u + 6, v + 6, 8);
    NN_DOTREV_S3_A3_1X1(b, a, res[6], 0, b, a, u + 7, v + 7, 7);
    NN_DOTREV_S3_A3_1X1(b, a, res[7], 0, b, a, u + 8, v + 8, 6);
    NN_DOTREV_S3_A3_1X1(b, a, res[8], 0, b, a, u + 9, v + 9, 5);
    NN_DOTREV_S3_A3_1X1(b, a, res[9], 0, b, a, u + 10, v + 10, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[10], 0, b, a, u + 11, v + 11, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[11], 0, b, a, u + 12, v + 12, 2);
    NN_ADDMUL_S2_A2_1X1(res[13], res[12], b, a, u[13], v[13]);
    return low;
}

mp_limb_t flint_mpn_mulhigh_15(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a, low;
    NN_DOTREV_S3_1X1_HIGH(b, a, u, v, 14);
    NN_DOTREV_S3_A3_1X1(b, a, low, 0, b, a, u, v, 15);
    NN_DOTREV_S3_A3_1X1(b, a, res[0], 0, b, a, u + 1, v + 1, 14);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, b, a, u + 2, v + 2, 13);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u + 3, v + 3, 12);
    NN_DOTREV_S3_A3_1X1(b, a, res[3], 0, b, a, u + 4, v + 4, 11);
    NN_DOTREV_S3_A3_1X1(b, a, res[4], 0, b, a, u + 5, v + 5, 10);
    NN_DOTREV_S3_A3_1X1(b, a, res[5], 0, b, a, u + 6, v + 6, 9);
    NN_DOTREV_S3_A3_1X1(b, a, res[6], 0, b, a, u + 7, v + 7, 8);
    NN_DOTREV_S3_A3_1X1(b, a, res[7], 0, b, a, u + 8, v + 8, 7);
    NN_DOTREV_S3_A3_1X1(b, a, res[8], 0, b, a, u + 9, v + 9, 6);
    NN_DOTREV_S3_A3_1X1(b, a, res[9], 0, b, a, u + 10, v + 10, 5);
    NN_DOTREV_S3_A3_1X1(b, a, res[10], 0, b, a, u + 11, v + 11, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[11], 0, b, a, u + 12, v + 12, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[12], 0, b, a, u + 13, v + 13, 2);
    NN_ADDMUL_S2_A2_1X1(res[14], res[13], b, a, u[14], v[14]);
    return low;
}

mp_limb_t flint_mpn_mulhigh_16(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a, low;
    NN_DOTREV_S3_1X1_HIGH(b, a, u, v, 15);
    NN_DOTREV_S3_A3_1X1(b, a, low, 0, b, a, u, v, 16);
    NN_DOTREV_S3_A3_1X1(b, a, res[0], 0, b, a, u + 1, v + 1, 15);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, b, a, u + 2, v + 2, 14);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u + 3, v + 3, 13);
    NN_DOTREV_S3_A3_1X1(b, a, res[3], 0, b, a, u + 4, v + 4, 12);
    NN_DOTREV_S3_A3_1X1(b, a, res[4], 0, b, a, u + 5, v + 5, 11);
    NN_DOTREV_S3_A3_1X1(b, a, res[5], 0, b, a, u + 6, v + 6, 10);
    NN_DOTREV_S3_A3_1X1(b, a, res[6], 0, b, a, u + 7, v + 7, 9);
    NN_DOTREV_S3_A3_1X1(b, a, res[7], 0, b, a, u + 8, v + 8, 8);
    NN_DOTREV_S3_A3_1X1(b, a, res[8], 0, b, a, u + 9, v + 9, 7);
    NN_DOTREV_S3_A3_1X1(b, a, res[9], 0, b, a, u + 10, v + 10, 6);
    NN_DOTREV_S3_A3_1X1(b, a, res[10], 0, b, a, u + 11, v + 11, 5);
    NN_DOTREV_S3_A3_1X1(b, a, res[11], 0, b, a, u + 12, v + 12, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[12], 0, b, a, u + 13, v + 13, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[13], 0, b, a, u + 14, v + 14, 2);
    NN_ADDMUL_S2_A2_1X1(res[15], res[14], b, a, u[15], v[15]);
    return low;
}

const flint_mpn_mul_func_t flint_mpn_mulhigh_func_tab[] =
{
    NULL,
    flint_mpn_mulhigh_1,
    flint_mpn_mulhigh_2,
    flint_mpn_mulhigh_3,
    flint_mpn_mulhigh_4,
    flint_mpn_mulhigh_5,
    flint_mpn_mulhigh_6,
    flint_mpn_mulhigh_7,
    flint_mpn_mulhigh_8,
    flint_mpn_mulhigh_9,
    flint_mpn_mulhigh_10,
    flint_mpn_mulhigh_11,
    flint_mpn_mulhigh_12,
    flint_mpn_mulhigh_13,
    flint_mpn_mulhigh_14,
    flint_mpn_mulhigh_15,
    flint_mpn_mulhigh_16
};

const flint_mpn_mulhigh_normalised_func_t flint_mpn_mulhigh_normalised_func_tab[] =
{
    NULL,
};

mp_limb_t flint_mpn_sqrhigh_1(mp_ptr res, mp_srcptr u)
{
    mp_limb_t low;
    umul_ppmm(res[0], low, u[0], u[0]);
    return low;
}

/* todo */
mp_limb_t flint_mpn_sqrhigh_2(mp_ptr res, mp_srcptr u)
{
    mp_limb_t b, low;
    FLINT_MPN_MUL_2X2(res[1], res[0], low, b, u[1], u[0], u[1], u[0]);
    return low;
}

/* todo: higher cases */

const flint_mpn_sqr_func_t flint_mpn_sqrhigh_func_tab[] = {
    NULL,
    flint_mpn_sqrhigh_1,
    flint_mpn_sqrhigh_2,
};

#endif
