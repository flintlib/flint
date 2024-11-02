/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ARB_MAT_IMPL_H
#define ARB_MAT_IMPL_H

#include "arb_types.h"

#ifdef __cplusplus
extern "C" {
#endif

slong _arb_mat_exp_choose_N(const mag_t norm, slong prec);

#ifdef __cplusplus
}
#endif

#endif /* ARB_MAT_IMPL_H */
