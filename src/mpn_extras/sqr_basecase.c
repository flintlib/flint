/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

#if FLINT_HAVE_ASSEMBLY_x86_64_adx

mp_limb_t flint_mpn_sqr_1(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqr_2(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqr_3(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqr_4(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqr_5(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqr_6(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqr_7(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqr_8(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqr_9(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqr_10(mp_ptr, mp_srcptr);

mp_limb_t flint_mpn_mul_7_2(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mul_7_3(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mul_7_4(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mul_7_5(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mul_7_6(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mul_7_7(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mul_8_7(mp_ptr, mp_srcptr, mp_srcptr);

mp_limb_t flint_mpn_sqr_11(mp_ptr res, mp_srcptr x)
{
    mp_limb_t t[14];
    flint_mpn_sqr_7(res, x);
    flint_mpn_sqr_4(res + 14, x + 7);
    flint_mpn_mul_7_4(t, x, x + 7);
    t[11] = 0;
    t[12] = 0;
    t[13] = 0;
    res[21] += mpn_addmul_1(res + 7, t, 14, 2);
    return res[21];
}

mp_limb_t flint_mpn_sqr_12(mp_ptr res, mp_srcptr x)
{
    mp_limb_t t[12];
    mp_limb_t cy;
    flint_mpn_sqr_7(res, x);
    flint_mpn_sqr_5(res + 14, x + 7);
    flint_mpn_mul_7_5(t, x, x + 7);
    cy = mpn_addmul_1(res + 7, t, 12, 2);
    mpn_add_1(res + 19, res + 19, 5, cy);
    return res[23];
}

mp_limb_t flint_mpn_sqr_13(mp_ptr res, mp_srcptr x)
{
    mp_limb_t t[13];
    mp_limb_t cy;
    flint_mpn_sqr_7(res, x);
    flint_mpn_sqr_6(res + 14, x + 7);
    flint_mpn_mul_7_6(t, x, x + 7);
    cy = mpn_addmul_1(res + 7, t, 13, 2);
    mpn_add_1(res + 20, res + 20, 6, cy);
    return res[25];
}

mp_limb_t flint_mpn_sqr_14(mp_ptr res, mp_srcptr x)
{
    mp_limb_t t[14];
    mp_limb_t cy;
    flint_mpn_sqr_7(res, x);
    flint_mpn_sqr_7(res + 14, x + 7);
    flint_mpn_mul_7_7(t, x, x + 7);
    cy = mpn_addmul_1(res + 7, t, 14, 2);
    mpn_add_1(res + 21, res + 21, 7, cy);
    return res[27];
}

const flint_mpn_sqr_func_t flint_mpn_sqr_func_tab[] = {
    NULL,
    flint_mpn_sqr_1,
    flint_mpn_sqr_2,
    flint_mpn_sqr_3,
    flint_mpn_sqr_4,
    flint_mpn_sqr_5,
    flint_mpn_sqr_6,
    flint_mpn_sqr_7,
    flint_mpn_sqr_8,
    flint_mpn_sqr_9,
    flint_mpn_sqr_10,
    flint_mpn_sqr_11,
    flint_mpn_sqr_12,
    flint_mpn_sqr_13,
    flint_mpn_sqr_14
};

#elif FLINT_HAVE_ASSEMBLY_armv8
mp_limb_t flint_mpn_sqr_1(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqr_2(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqr_3(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqr_4(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqr_5(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqr_6(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqr_7(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqr_8(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqr_9(mp_ptr, mp_srcptr);

const flint_mpn_sqr_func_t flint_mpn_sqr_func_tab[] =
{
    NULL,
    flint_mpn_sqr_1,
    flint_mpn_sqr_2,
    flint_mpn_sqr_3,
    flint_mpn_sqr_4,
    flint_mpn_sqr_5,
    flint_mpn_sqr_6,
    flint_mpn_sqr_7,
    flint_mpn_sqr_8,
    flint_mpn_sqr_9
};
#else

/* Currently generic C code performs worse than GMP. */

const flint_mpn_sqr_func_t flint_mpn_sqr_func_tab[] = {
    NULL,
};

#endif
