/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

#if FLINT_HAVE_ADX

mp_limb_t flint_mpn_sqr_1(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqr_2(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqr_3(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqr_4(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqr_5(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqr_6(mp_ptr, mp_srcptr);
mp_limb_t flint_mpn_sqr_7(mp_ptr, mp_srcptr);

static mp_limb_t (* sqr_funcs[])(mp_ptr, mp_srcptr) =
{
    flint_mpn_sqr_1,
    flint_mpn_sqr_2,
    flint_mpn_sqr_3,
    flint_mpn_sqr_4,
    flint_mpn_sqr_5,
    flint_mpn_sqr_6,
    flint_mpn_sqr_7
};

/* NOTE: Important that the inputs are in this order to avoid additional moves
 * when calling `mul_funcs'. */
mp_limb_t
flint_mpn_sqr_basecase(mp_ptr res, mp_srcptr ap, slong alen)
{
    FLINT_ASSERT(alen > 0);
    FLINT_ASSERT(alen < 8);

    return sqr_funcs[alen - 1](res, ap);
}

#else

typedef int this_file_is_empty;

#endif
