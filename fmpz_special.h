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

#ifndef FMPZ_SPECIAL_H
#define FMPZ_SPECIAL_H

#ifdef __cplusplus
 extern "C" {
#endif

/* Primorials */

void fmpz_primorial(fmpz_t res, ulong n);

/* Multiplicative functions */

void _fmpz_euler_phi(fmpz_t res, const fmpz_factor_t fac);
void fmpz_euler_phi(fmpz_t res, const fmpz_t n);

int _fmpz_moebius_mu(const fmpz_factor_t fac);
int fmpz_moebius_mu(const fmpz_t n);

void _fmpz_divisor_sigma(fmpz_t res, const fmpz_factor_t fac, ulong k);
void fmpz_divisor_sigma(fmpz_t res, const fmpz_t n, ulong k);

/* Orthogonal polynomials */

void fmpz_chebyshev_t(fmpz_t y, ulong n, const fmpz_t x);
void fmpz_chebyshev_u(fmpz_t y, ulong n, const fmpz_t x);

#ifdef __cplusplus
}
#endif

#endif

