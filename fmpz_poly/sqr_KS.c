/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010, 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void
_fmpz_poly_sqr_KS(fmpz *rop, const fmpz *op, slong len)
{
    const slong in_len = len;
    int neg;
    slong bits, limbs, loglen;
    mp_limb_t *arr, *arr3;
    slong sign = 0;

    FMPZ_VEC_NORM(op, len);

    if (!len)
    {
        if (2 * in_len - 1 > 0)
            _fmpz_vec_zero(rop, 2 * in_len - 1);
        return;
    }

    neg = (fmpz_sgn(op + len - 1) > 0) ? 0 : -1;

    bits = _fmpz_vec_max_bits(op, len);
    if (bits < 0)
    {
        sign = 1;
        bits = - bits;
    }

    loglen = FLINT_BIT_COUNT(len);
    bits   = 2 * bits + loglen + sign;
    limbs  = (bits * len - 1) / FLINT_BITS + 1;

    arr = (mp_limb_t *) flint_calloc(limbs, sizeof(mp_limb_t));

    _fmpz_poly_bit_pack(arr, op, len, bits, neg);

    arr3 = (mp_limb_t *) flint_malloc((2 * limbs) * sizeof(mp_limb_t));

    mpn_sqr(arr3, arr, limbs);

    if (sign)
        _fmpz_poly_bit_unpack(rop, 2 * len - 1, arr3, bits, 0);
    else
        _fmpz_poly_bit_unpack_unsigned(rop, 2 * len - 1, arr3, bits);

    if (len < in_len)
        _fmpz_vec_zero(rop + (2 * len - 1), 2 * (in_len - len));

    flint_free(arr);
    flint_free(arr3);
}

void fmpz_poly_sqr_KS(fmpz_poly_t rop, const fmpz_poly_t op)
{
    slong len;

    if (op->length == 0)
    {
        fmpz_poly_zero(rop);
        return;
    }

    len = 2 * op->length - 1;

    if (rop == op)
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, len);
        fmpz_poly_sqr_KS(t, op);
        fmpz_poly_swap(rop, t);
        fmpz_poly_clear(t);
        return;
    }

    fmpz_poly_fit_length(rop, len);
    _fmpz_poly_sqr_KS(rop->coeffs, op->coeffs, op->length);
    _fmpz_poly_set_length(rop, len);
}

