/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef QQBAR_IMPL_H
#define QQBAR_IMPL_H

#include "qqbar.h"

#ifdef __cplusplus
extern "C" {
#endif

void
qqbar_get_decimal_root_nearest(char ** re_s, char ** im_s, const qqbar_t x, slong default_digits);

ulong
qqbar_try_as_cyclotomic(qqbar_t zeta, fmpq_poly_t poly, const qqbar_t x);

void
best_rational_fast(slong * p, ulong * q, double x, slong N);

#ifdef __cplusplus
}
#endif

#endif /* QQBAR_IMPL_H */
