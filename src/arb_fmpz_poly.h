/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ARB_FMPZ_POLY_H
#define ARB_FMPZ_POLY_H

#ifdef ARB_FMPZ_POLY_INLINES_C
#define ARB_FMPZ_POLY_INLINE
#else
#define ARB_FMPZ_POLY_INLINE static inline
#endif

#include "fmpz_poly.h"
#include "acb_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define ARB_FMPZ_POLY_ROOTS_VERBOSE 1

void _arb_fmpz_poly_evaluate_acb_horner(acb_t res, const fmpz * f, slong len, const acb_t x, slong prec);
void arb_fmpz_poly_evaluate_acb_horner(acb_t res, const fmpz_poly_t f, const acb_t a, slong prec);
void _arb_fmpz_poly_evaluate_acb_rectangular(acb_t res, const fmpz * f, slong len, const acb_t x, slong prec);
void arb_fmpz_poly_evaluate_acb_rectangular(acb_t res, const fmpz_poly_t f, const acb_t a, slong prec);
void _arb_fmpz_poly_evaluate_acb(acb_t res, const fmpz * f, slong len, const acb_t x, slong prec);
void arb_fmpz_poly_evaluate_acb(acb_t res, const fmpz_poly_t f, const acb_t a, slong prec);

void _arb_fmpz_poly_evaluate_arb_horner(arb_t res, const fmpz * f, slong len, const arb_t x, slong prec);
void arb_fmpz_poly_evaluate_arb_horner(arb_t res, const fmpz_poly_t f, const arb_t a, slong prec);
void _arb_fmpz_poly_evaluate_arb_rectangular(arb_t res, const fmpz * f, slong len, const arb_t x, slong prec);
void arb_fmpz_poly_evaluate_arb_rectangular(arb_t res, const fmpz_poly_t f, const arb_t a, slong prec);
void _arb_fmpz_poly_evaluate_arb(arb_t res, const fmpz * f, slong len, const arb_t x, slong prec);
void arb_fmpz_poly_evaluate_arb(arb_t res, const fmpz_poly_t f, const arb_t a, slong prec);

void arb_fmpz_poly_deflate(fmpz_poly_t result, const fmpz_poly_t input, ulong deflation);
ulong arb_fmpz_poly_deflation(const fmpz_poly_t input);

void arb_fmpz_poly_complex_roots(acb_ptr roots, const fmpz_poly_t poly, int flags, slong target_prec);

ARB_FMPZ_POLY_INLINE
void arb_fmpz_poly_cos_minpoly(fmpz_poly_t res, ulong n)
{
    fmpz_poly_cos_minpoly(res, n);
}

void arb_fmpz_poly_gauss_period_minpoly(fmpz_poly_t res, ulong q, ulong n);

#ifdef __cplusplus
}
#endif

#endif
