/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef MPN_MOD_IMPL_H
#define MPN_MOD_IMPL_H

#include "flint.h"

#ifdef __cplusplus
extern "C" {
#endif

void _mpn_dot_rev_2x2_3(nn_ptr s, nn_srcptr a, nn_srcptr b, slong len);
void _mpn_dot_rev_2x2_4(nn_ptr s, nn_srcptr a, nn_srcptr b, slong len);
void _mpn_dot_rev_2x2_5(nn_ptr s, nn_srcptr a, nn_srcptr b, slong len);

void _mpn_dot_rev_3x3_5(nn_ptr s, nn_srcptr a, nn_srcptr b, slong len);

void _mpn_dot_rev_nxn_2n(nn_ptr res, nn_srcptr a, nn_srcptr b, slong len, slong nlimbs);

void _mpn_dot_rev_nxn_2nm1(nn_ptr res, nn_srcptr a, nn_srcptr b, slong len, slong nlimbs);

void _mpn_dot_rev_nxn_2np1(nn_ptr res, nn_srcptr a, nn_srcptr b, slong len, slong nlimbs);

#ifdef __cplusplus
}
#endif

#endif /* MPN_MOD_IMPL_H */
