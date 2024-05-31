/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

#if FLINT_HAVE_ASSEMBLY_x86_64_adx || FLINT_HAVE_ASSEMBLY_armv8
mp_limb_t flint_mpn_sqrhigh_1(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqrhigh_2(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqrhigh_3(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqrhigh_4(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqrhigh_5(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqrhigh_6(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqrhigh_7(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqrhigh_8(mp_ptr, mp_srcptr);

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
#else
mp_limb_t flint_mpn_sqrhigh_1(mp_ptr res, mp_srcptr u)
{
    mp_limb_t low;
    umul_ppmm(res[0], low, u[0], u[0]);
    return low;
}

mp_limb_t flint_mpn_sqrhigh_2(mp_ptr res, mp_srcptr u)
{
    mp_limb_t b, low;
    FLINT_MPN_SQR_2X2(res[1], res[0], low, b, u[1], u[0]);
    return low;
}

/* todo: higher cases */

const flint_mpn_sqr_func_t flint_mpn_sqrhigh_func_tab[] = {
    NULL,
    flint_mpn_sqrhigh_1,
    flint_mpn_sqrhigh_2,
};
#endif

#if FLINT_HAVE_ASSEMBLY_x86_64_adx
mp_limb_pair_t flint_mpn_sqrhigh_normalised_1(mp_ptr, mp_srcptr);
mp_limb_pair_t flint_mpn_sqrhigh_normalised_2(mp_ptr, mp_srcptr);
mp_limb_pair_t flint_mpn_sqrhigh_normalised_3(mp_ptr, mp_srcptr);
mp_limb_pair_t flint_mpn_sqrhigh_normalised_4(mp_ptr, mp_srcptr);
mp_limb_pair_t flint_mpn_sqrhigh_normalised_5(mp_ptr, mp_srcptr);
mp_limb_pair_t flint_mpn_sqrhigh_normalised_6(mp_ptr, mp_srcptr);
mp_limb_pair_t flint_mpn_sqrhigh_normalised_7(mp_ptr, mp_srcptr);
mp_limb_pair_t flint_mpn_sqrhigh_normalised_8(mp_ptr, mp_srcptr);

const flint_mpn_sqrhigh_normalised_func_t flint_mpn_sqrhigh_normalised_func_tab[] =
{
    NULL,
    flint_mpn_sqrhigh_normalised_1,
    flint_mpn_sqrhigh_normalised_2,
    flint_mpn_sqrhigh_normalised_3,
    flint_mpn_sqrhigh_normalised_4,
    flint_mpn_sqrhigh_normalised_5,
    flint_mpn_sqrhigh_normalised_6,
    flint_mpn_sqrhigh_normalised_7,
    flint_mpn_sqrhigh_normalised_8
};
#else
const flint_mpn_sqrhigh_normalised_func_t flint_mpn_sqrhigh_normalised_func_tab[] =
{
    NULL
};
#endif
