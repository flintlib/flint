/*
    Copyright (C) 2025 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef QQBAR_IMPL_H
#define QQBAR_IMPL_H

#include "qqbar.h"
#include "fmpz_poly_factor.h"

void qqbar_get_decimal_root_nearest(char ** re_s, char ** im_s, const qqbar_t x, slong default_digits);
void best_rational_fast(slong * p, ulong * q, double x, slong N);
ulong qqbar_try_as_cyclotomic(qqbar_t zeta, fmpq_poly_t poly, const qqbar_t x);

void _qqbar_binary_op_quadratic(qqbar_t res, const qqbar_t x, const qqbar_t y, int op);
int _qqbar_express_in_field_search(fmpq_poly_t R, const qqbar_t x, const qqbar_t y, slong min_prec, slong max_prec);
int _qqbar_binary_op_structured(qqbar_t res, const qqbar_t x, const qqbar_t y, int op);
int _qqbar_subfield_modular_test(const fmpz_poly_t P, const fmpz_poly_t Q, int same);
int _qqbar_binary_op_subfield(qqbar_t res, const qqbar_t x, const qqbar_t y, int op, int swapped, slong min_prec, slong max_prec, int test);

void _qqbar_binary_op_select_factor(qqbar_t res, const fmpz_poly_factor_t fac, const qqbar_t x, const qqbar_t y, int op);
void _qqbar_factor_inflate2_irreducible(fmpz_poly_factor_t fac, const fmpz_poly_t T);
void _qqbar_conjugate_pair_op(qqbar_t res, const qqbar_t x, const qqbar_t y, int op);

#endif
