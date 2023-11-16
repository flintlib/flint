/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef PARTITIONS_H
#define PARTITIONS_H

#include "arb_types.h"

#ifdef __cplusplus
extern "C" {
#endif

void partitions_rademacher_bound(arf_t b, const fmpz_t n, ulong N);

void partitions_hrr_sum_arb(arb_t x, const fmpz_t n, slong N0, slong N, int use_doubles);

void partitions_fmpz_fmpz(fmpz_t p, const fmpz_t n, int use_doubles);

void partitions_fmpz_ui(fmpz_t p, ulong n);

/* deprecated */
#define partitions_fmpz_ui_using_doubles partitions_fmpz_ui

void partitions_leading_fmpz(arb_t res, const fmpz_t n, slong prec);

#ifdef __cplusplus
}
#endif

#endif
