/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

#if FLINT_HAVE_ASSEMBLY_x86_64_adx
mp_limb_t flint_mpn_mullow_1(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mullow_2(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mullow_3(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mullow_4(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mullow_5(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mullow_6(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mullow_7(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mullow_8(mp_ptr, mp_srcptr, mp_srcptr);

const flint_mpn_mul_func_t flint_mpn_mullow_func_tab[] =
{
    NULL,
    flint_mpn_mullow_1,
    flint_mpn_mullow_2,
    flint_mpn_mullow_3,
    flint_mpn_mullow_4,
    flint_mpn_mullow_5,
    flint_mpn_mullow_6,
    flint_mpn_mullow_7,
    flint_mpn_mullow_8
};
#else
const flint_mpn_mul_func_t flint_mpn_mullow_func_tab[] = { NULL };
#endif

#if !FLINT_HAVE_NATIVE_mpn_mullow_basecase
mp_limb_t
flint_mpn_mullow_basecase(mp_ptr rp, mp_srcptr xp, mp_srcptr yp, mp_size_t n)
{
    mp_limb_t ret;
    mp_size_t ix;

    ret = mpn_mul_1(rp, xp, n, yp[0]);

    for (ix = 1; ix < n; ix++)
    {
        ret += mpn_addmul_1(rp + ix, xp, n - ix, yp[ix]);
        ret += xp[n - ix] * yp[ix];
    }

    return ret;
}
#endif
