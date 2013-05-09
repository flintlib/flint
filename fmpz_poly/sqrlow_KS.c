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

    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010, 2011 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void _fmpz_poly_sqrlow_KS(fmpz * res, const fmpz * poly, len_t len, len_t n)
{
    int neg;
    len_t bits, limbs, loglen, sign = 0;
    mp_limb_t *arr_in, *arr_out;

    FMPZ_VEC_NORM(poly, len);

    if (len == 0)
    {
        _fmpz_vec_zero(res, n);
        return;
    }

    neg = (fmpz_sgn(poly + len - 1) > 0) ? 0 : -1;

    if (n > 2 * len - 1)
    {
       _fmpz_vec_zero(res + 2 * len - 1, n - (2 * len - 1));
       n = 2 * len - 1;
    }

    bits = _fmpz_vec_max_bits(poly, len);
    if (bits < 0)
    {
        sign = 1;
        bits = - bits;
    }

    loglen = FLINT_BIT_COUNT(len);
    bits   = 2 * bits + loglen + sign;
    limbs  = (bits * len - 1) / FLINT_BITS + 1;

    arr_in  = flint_calloc(limbs, sizeof(mp_limb_t));
    arr_out = flint_malloc((2 * limbs) * sizeof(mp_limb_t));

    _fmpz_poly_bit_pack(arr_in, poly, len, bits, neg);

    mpn_sqr(arr_out, arr_in, limbs);

    if (sign)
        _fmpz_poly_bit_unpack(res, n, arr_out, bits, 0);
    else
        _fmpz_poly_bit_unpack_unsigned(res, n, arr_out, bits);

    flint_free(arr_in);
    flint_free(arr_out);
}

void fmpz_poly_sqrlow_KS(fmpz_poly_t res, const fmpz_poly_t poly, len_t n)
{
    const len_t len = poly->length;

    if (len == 0 || n == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    if (res == poly)
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, n);
        fmpz_poly_sqrlow_KS(t, poly, n);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
        return;
    }

    fmpz_poly_fit_length(res, n);
    _fmpz_poly_sqrlow_KS(res->coeffs, poly->coeffs, len, n);
    _fmpz_poly_set_length(res, n);
    _fmpz_poly_normalise(res);
}

