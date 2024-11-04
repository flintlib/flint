/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_POLY_IMPL_H
#define FMPZ_POLY_IMPL_H

#include "fmpz_types.h"

#ifdef __cplusplus
extern "C" {
#endif

double _fmpz_poly_evaluate_horner_d_2exp2_precomp(slong * exp, const double * poly, const slong * poly_exp, slong n, double d, slong dexp);

void revbin1(fmpz * out, const fmpz * in, slong len, slong bits);
void revbin2(fmpz * out, const fmpz * in, slong len, slong bits);

void _fmpz_vec_add_rev(fmpz * in1, fmpz * in2, slong bits);

#ifdef __cplusplus
}
#endif

#endif /* FMPZ_POLY_IMPL_H */
