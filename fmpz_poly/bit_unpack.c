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
   
******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

int
_fmpz_poly_bit_unpack(fmpz * poly, long len,
                      mp_srcptr arr, mp_bitcnt_t bit_size, int negate)
{
    mp_bitcnt_t bits = 0;
    mp_size_t limbs = 0;
    mp_bitcnt_t b = bit_size % FLINT_BITS;
    mp_size_t l = bit_size / FLINT_BITS;
    int borrow = 0;
    long i;

    for (i = 0; i < len; i++)
    {
        borrow =
            fmpz_bit_unpack(poly + i, arr + limbs, bits, bit_size, negate,
                            borrow);
        limbs += l;
        bits += b;
        if (bits >= FLINT_BITS)
        {
            bits -= FLINT_BITS;
            limbs++;
        }
    }

    return borrow;
}

void
_fmpz_poly_bit_unpack_unsigned(fmpz * poly, long len,
                               mp_srcptr arr, mp_bitcnt_t bit_size)
{
    mp_bitcnt_t bits = 0;
    mp_size_t limbs = 0;
    mp_bitcnt_t b = bit_size % FLINT_BITS;
    mp_size_t l = bit_size / FLINT_BITS;
    long i;

    for (i = 0; i < len; i++)
    {
        fmpz_bit_unpack_unsigned(poly + i, arr + limbs, bits, bit_size);
        limbs += l;
        bits += b;
        if (bits >= FLINT_BITS)
        {
            bits -= FLINT_BITS;
            limbs++;
        }
    }
}

void
fmpz_poly_bit_unpack_unsigned(fmpz_poly_t poly, const fmpz_t f,
                                        mp_bitcnt_t bit_size)
{
    long len;
    mpz_t tmp;

    if (fmpz_sgn(f) < 0)
    {
        printf("Exception (fmpz_poly_bit_unpack_unsigned). Expected an unsigned value.\n");
        abort();
    }

    if (bit_size == 0 || fmpz_is_zero(f))
    {
        fmpz_poly_zero(poly);
        return;
    }

    len = (fmpz_bits(f) + bit_size - 1) / bit_size;

    mpz_init2(tmp, bit_size*len);
    flint_mpn_zero(tmp->_mp_d, tmp->_mp_alloc);
    fmpz_get_mpz(tmp, f);

    fmpz_poly_fit_length(poly, len);

    _fmpz_poly_bit_unpack_unsigned(poly->coeffs, len, tmp->_mp_d, bit_size);
    _fmpz_poly_set_length(poly, len);
    _fmpz_poly_normalise(poly);

    mpz_clear(tmp);
}


void
fmpz_poly_bit_unpack(fmpz_poly_t poly, const fmpz_t f, mp_bitcnt_t bit_size)
{
    long len;
    mpz_t tmp;
    int negate, borrow;

    if (bit_size == 0 || fmpz_is_zero(f))
    {
        fmpz_poly_zero(poly);
        return;
    }

    /* Round up */
    len = (fmpz_bits(f) + bit_size - 1) / bit_size;
    negate = (fmpz_sgn(f) < 0) ? -1 : 0;

    mpz_init2(tmp, bit_size*len);

    /* TODO: avoid all this wastefulness */
    flint_mpn_zero(tmp->_mp_d, tmp->_mp_alloc);
    fmpz_get_mpz(tmp, f);

    fmpz_poly_fit_length(poly, len + 1);

    borrow = _fmpz_poly_bit_unpack(poly->coeffs, len,
                    tmp->_mp_d, bit_size, negate);

    if (borrow)
    {
        fmpz_set_si(poly->coeffs + len, negate ? -1L : 1L);
        _fmpz_poly_set_length(poly, len + 1);
    }
    else
    {
        _fmpz_poly_set_length(poly, len);
        _fmpz_poly_normalise(poly);
    }

    mpz_clear(tmp);
}
