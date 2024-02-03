/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
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

const flint_mpn_sqr_func_t flint_mpn_sqr_func_tab[] = {
    NULL,
    flint_mpn_sqr_1,
    flint_mpn_sqr_2,
    flint_mpn_sqr_3,
    flint_mpn_sqr_4,
    flint_mpn_sqr_5,
    flint_mpn_sqr_6,
    flint_mpn_sqr_7
};

#else

/* Currently generic C code performs worse than GMP. */

const flint_mpn_sqr_func_t flint_mpn_sqr_func_tab[] = {
    NULL,
};

#endif

mp_limb_t
flint_mpn_sqr_basecase(mp_ptr r, mp_srcptr x, mp_size_t n)
{
    return flint_mpn_sqr_func_tab[n](r, x);
}
