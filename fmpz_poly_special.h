/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#ifndef FMPZ_POLY_SPECIAL_H
#define FMPZ_POLY_SPECIAL_H

#ifdef __cplusplus
 extern "C" {
#endif

void _fmpz_poly_cyclotomic(fmpz * a, ulong n, mp_ptr factors,
                                        slong num_factors, ulong phi);
void fmpz_poly_cyclotomic(fmpz_poly_t poly, ulong n);

void _fmpz_poly_cos_minpoly(fmpz * f, ulong n);
void fmpz_poly_cos_minpoly(fmpz_poly_t f, ulong n);

void _fmpz_poly_swinnerton_dyer(fmpz * T, ulong n);
void fmpz_poly_swinnerton_dyer(fmpz_poly_t poly, ulong n);

void _fmpz_poly_chebyshev_t(fmpz * coeffs, ulong n);
void fmpz_poly_chebyshev_t(fmpz_poly_t poly, ulong n);

void _fmpz_poly_chebyshev_u(fmpz * coeffs, ulong n);
void fmpz_poly_chebyshev_u(fmpz_poly_t poly, ulong n);

void _fmpz_poly_eta_qexp(fmpz * f, slong e, slong n);
void fmpz_poly_eta_qexp(fmpz_poly_t f, slong e, slong n);

void _fmpz_poly_theta_qexp(fmpz * f, slong e, slong n);
void fmpz_poly_theta_qexp(fmpz_poly_t f, slong e, slong n);

#ifdef __cplusplus
}
#endif

#endif

