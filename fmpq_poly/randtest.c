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

    Copyright (C) 2010 William Hart
    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq_poly.h"
#include "ulong_extras.h"

void fmpq_poly_randtest(fmpq_poly_t poly, flint_rand_t state, 
                        len_t len, mp_bitcnt_t bits)
{
    ulong m;

    m = n_randlimb(state);

    fmpq_poly_fit_length(poly, len);
    _fmpq_poly_set_length(poly, len);

    if (m & 1UL)
    {
        _fmpz_vec_randtest(poly->coeffs, state, len, bits);
    }
    else
    {
        fmpz_t x;

        fmpz_init(x);
        fmpz_randtest(x, state, bits / 2);
        _fmpz_vec_randtest(poly->coeffs, state, len, (bits + 1) / 2);
        _fmpz_vec_scalar_mul_fmpz(poly->coeffs, poly->coeffs, len, x);
        fmpz_clear(x);
    }

    if (m & 2UL)
    {
        fmpz_randtest_not_zero(poly->den, state, FLINT_MAX(bits, 1));
        fmpz_abs(poly->den, poly->den);
        fmpq_poly_canonicalise(poly);
    }
    else
    {
        fmpz_one(poly->den);
        _fmpq_poly_normalise(poly);
    }
}

void fmpq_poly_randtest_unsigned(fmpq_poly_t poly, flint_rand_t state,
                                 len_t len, mp_bitcnt_t bits)
{
    ulong m;

    m = n_randlimb(state);

    fmpq_poly_fit_length(poly, len);
    _fmpq_poly_set_length(poly, len);

    if (m & 1UL)
    {
        _fmpz_vec_randtest_unsigned(poly->coeffs, state, len, bits);
    }
    else
    {
        fmpz_t x;

        fmpz_init(x);
        fmpz_randtest_unsigned(x, state, bits / 2);
        _fmpz_vec_randtest_unsigned(poly->coeffs, state, len, (bits + 1) / 2);
        _fmpz_vec_scalar_mul_fmpz(poly->coeffs, poly->coeffs, len, x);
        fmpz_clear(x);
    }

    if (m & 2UL)
    {
        fmpz_randtest_not_zero(poly->den, state, FLINT_MAX(bits, 1));
        fmpz_abs(poly->den, poly->den);
        fmpq_poly_canonicalise(poly);
    }
    else
    {
        fmpz_one(poly->den);
        _fmpq_poly_normalise(poly);
    }
}

void fmpq_poly_randtest_not_zero(fmpq_poly_t f, flint_rand_t state, 
                                 len_t len, mp_bitcnt_t bits)
{
    if ((bits == 0) | (len == 0))
    {
        printf("Exception (fmpq_poly_randtest_not_zeo). bits == 0.\n");
        abort();
    }

    fmpq_poly_randtest(f, state, len, bits);
    if (f->length == 0) 
        fmpq_poly_set_ui(f, 1);
}

