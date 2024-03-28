/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ARB_FMPZ_POLY_H
#define ARB_FMPZ_POLY_H

#include "fmpz_poly.h"
#include "acb_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define ARB_FMPZ_POLY_ROOTS_VERBOSE 1

void _arb_fmpz_poly_evaluate_acb_horner(acb_t res, const fmpz * f, slong len, const acb_t x, slong prec);
void arb_fmpz_poly_evaluate_acb_horner(acb_t res, const fmpz_poly_t poly, const acb_t x, slong prec);
void _arb_fmpz_poly_evaluate_acb_rectangular(acb_t res, const fmpz * f, slong len, const acb_t x, slong prec);
void arb_fmpz_poly_evaluate_acb_rectangular(acb_t res, const fmpz_poly_t poly, const acb_t x, slong prec);
void _arb_fmpz_poly_evaluate_acb(acb_t res, const fmpz * f, slong len, const acb_t x, slong prec);
void arb_fmpz_poly_evaluate_acb(acb_t res, const fmpz_poly_t poly, const acb_t x, slong prec);

void _arb_fmpz_poly_evaluate_arb_horner(arb_t res, const fmpz * f, slong len, const arb_t x, slong prec);
void arb_fmpz_poly_evaluate_arb_horner(arb_t res, const fmpz_poly_t poly, const arb_t x, slong prec);
void _arb_fmpz_poly_evaluate_arb_rectangular(arb_t res, const fmpz * f, slong len, const arb_t x, slong prec);
void arb_fmpz_poly_evaluate_arb_rectangular(arb_t res, const fmpz_poly_t poly, const arb_t x, slong prec);
void _arb_fmpz_poly_evaluate_arb(arb_t res, const fmpz * f, slong len, const arb_t x, slong prec);
void arb_fmpz_poly_evaluate_arb(arb_t res, const fmpz_poly_t poly, const arb_t x, slong prec);

void arb_fmpz_poly_deflate(fmpz_poly_t res, const fmpz_poly_t poly, ulong deflation);
ulong arb_fmpz_poly_deflation(const fmpz_poly_t poly);

void arb_fmpz_poly_complex_roots(acb_ptr roots, const fmpz_poly_t poly, int flags, slong prec);

void arb_fmpz_poly_gauss_period_minpoly(fmpz_poly_t res, ulong q, ulong n);

#ifdef __cplusplus
}
#endif

#endif
