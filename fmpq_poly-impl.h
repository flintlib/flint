/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPQ_POLY_IMPL_H
#define FMPQ_POLY_IMPL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include "fmpq_mat.h"
#include "fmpq_poly.h"

/* For _fmpz_gcd_ui */
#include "fmpq-impl.h"

/* defined in revert_series_lagrange.c */
void
_set_vec(fmpz * rnum, fmpz_t den, const fmpz * xnum, const fmpz * xden, slong len);

#endif
