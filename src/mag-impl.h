/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef MAG_IMPL_H
#define MAG_IMPL_H

#include "arb_types.h"

#ifdef __cplusplus
extern "C" {
#endif

void _mag_exp_d(mag_t res, double x, int roundup);

void mag_exp_huge(mag_t res, const mag_t x);

void mag_exp_huge_lower(mag_t res, const mag_t x);

#ifdef __cplusplus
}
#endif

#endif /* MAG_IMPL_H */
