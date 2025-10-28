/*
    Copyright (C) 2025 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPQ_POLY_IMPL_H
#define FMPQ_POLY_IMPL_H

#include "fmpq_types.h"

void _fmpq_poly_sin_cos_series_basecase_can(fmpz * S, fmpz_t Sden, fmpz * C, fmpz_t Cden, const fmpz * A, const fmpz_t Aden, slong Alen, slong n, int can);

#endif
