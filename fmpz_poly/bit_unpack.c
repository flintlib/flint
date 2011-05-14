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

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

void
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
fmpz_poly_bit_unpack(fmpz_poly_t poly, const fmpz_t f, mp_bitcnt_t bit_size)
{
    long len;
    mpz_t tmp;
    __mpz_struct * mpz;
    int negate;

    if (COEFF_IS_MPZ(*f))
    {
        mpz = COEFF_TO_PTR(*f);
    }
    else
    {
        mpz_init(tmp);
        fmpz_get_mpz(tmp, f);
        mpz = tmp;
    }

    negate = (mpz->_mp_size < 0);

    len = fmpz_bits(f) / bit_size + 1;

    fmpz_poly_fit_length(poly, len);

    _fmpz_poly_bit_unpack(poly->coeffs, len, mpz->_mp_d, bit_size, negate);
    _fmpz_poly_set_length(poly, len);
    _fmpz_poly_normalise(poly);

    if (!COEFF_IS_MPZ(*f))
        mpz_clear(tmp);
}

void
fmpz_poly_bit_unpack_unsigned(fmpz_poly_t poly, const fmpz_t f, mp_bitcnt_t bit_size)
{
    long len;
    mpz_t tmp;
    __mpz_struct * mpz;
    int negate;

    if (COEFF_IS_MPZ(*f))
    {
        mpz = COEFF_TO_PTR(*f);
    }
    else
    {
        mpz_init(tmp);
        fmpz_get_mpz(tmp, f);
        mpz = tmp;
    }

    if (mpz->_mp_size < 0)
        abort();

    len = fmpz_bits(f) / bit_size + 1;

    fmpz_poly_fit_length(poly, len);

    _fmpz_poly_bit_unpack_unsigned(poly->coeffs, len, mpz->_mp_d, bit_size);
    _fmpz_poly_set_length(poly, len);
    _fmpz_poly_normalise(poly);

    if (!COEFF_IS_MPZ(*f))
        mpz_clear(tmp);
}