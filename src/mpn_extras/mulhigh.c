/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

#if FLINT_HAVE_ADX
mp_limb_t flint_mpn_mulhigh_1(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mulhigh_2(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mulhigh_3(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mulhigh_4(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mulhigh_5(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mulhigh_6(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mulhigh_7(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mulhigh_8(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mulhigh_9(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mulhigh_10(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mulhigh_11(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mulhigh_12(mp_ptr, mp_srcptr, mp_srcptr);

const flint_mpn_mul_func_t flint_mpn_mulhigh_n_func_tab[] =
{
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
    flint_mpn_mulhigh_12
};
#else
typedef int this_file_is_empty;
#endif
