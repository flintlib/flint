/*============================================================================

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

===============================================================================*/
/******************************************************************************

 Copyright (C) 2010-2011 Fredrik Johansson

******************************************************************************/

#ifndef ARITH_H
#define ARITH_H

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

typedef struct
{
    int sign;
    fmpz * p;
    fmpz * exp;
    long alloc;
    long length;
} fmpz_factor_struct;

typedef fmpz_factor_struct fmpz_factor_t[1];


void fmpz_factor_init(fmpz_factor_t factor);
void fmpz_factor_clear(fmpz_factor_t factor);
void fmpz_factor_print(const fmpz_factor_t factor);

void _fmpz_factor_si(fmpz_factor_t factor, long n);
void _fmpz_factor_fit_length(fmpz_factor_t factor, long len);
void _fmpz_factor_append_ui(fmpz_factor_t factor, ulong p, ulong exp);
void _fmpz_factor_set_length(fmpz_factor_t factor, long newlen);
void _fmpz_factor_extend_factor_n(fmpz_factor_t factor, ulong n);

void fmpz_factor(fmpz_factor_t factor, const fmpz_t n);
void fmpz_unfactor(fmpz_t n, const fmpz_factor_t factor);

void fmpz_primorial(fmpz_t res, long n);
void fmpz_poly_ramanujan_tau(fmpz_poly_t res, long n);
void fmpz_ramanujan_tau(fmpz_t res, const fmpz_t n);
void fmpz_divisors(fmpz_poly_t res, const fmpz_t n);
void fmpz_divisor_sigma(fmpz_t res, const fmpz_t n, ulong k);

int fmpz_moebius_mu(const fmpz_t n);
void fmpz_euler_phi(fmpz_t res, const fmpz_t n);

void _mpq_harmonic_tiny(mpq_t res, long n);
void _mpq_harmonic_balanced(mpq_t res, long a, long b);
void _mpq_harmonic_odd_balanced(mpq_t res, long n);
void mpq_harmonic(mpq_t res, long n);

void _fmpz_stirling2_powsum(fmpz_t s, long n, long k);
void _fmpz_stirling2_powsum_odd(fmpz_t , long n, long k);

void fmpz_stirling1u(fmpz_t s, long n, long k);
void fmpz_stirling1(fmpz_t s, long n, long k);
void fmpz_stirling2(fmpz_t s, long n, long k);

void fmpz_stirling1u_vec(fmpz * row, long n, long klen);
void fmpz_stirling1_vec(fmpz * row, long n, long klen);
void fmpz_stirling2_vec(fmpz * row, long n, long klen);

void fmpz_stirling1u_vec_next(fmpz * row, fmpz * prev, long n, long klen);
void fmpz_stirling1_vec_next(fmpz * row, fmpz * prev, long n, long klen);
void fmpz_stirling2_vec_next(fmpz * row, fmpz * prev, long n, long klen);

void fmpz_stirling1u_mat(fmpz ** rows, long n);
void fmpz_stirling1_mat(fmpz ** rows, long n);
void fmpz_stirling2_mat(fmpz ** rows, long n);

void _fmpz_bernoulli_vec_series(fmpz_t den, fmpz * b, long n);
void _fmpz_bernoulli_vec_recursive(fmpz_t den, fmpz * b, long n);
void fmpz_bernoulli_vec(fmpz_t den, fmpz * num, long n);

void fmpq_poly_bernoulli(fmpq_poly_t poly, long n);

#endif
