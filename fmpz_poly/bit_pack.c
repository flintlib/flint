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

void
_fmpz_poly_bit_pack(mp_ptr arr, const fmpz * poly, len_t len,
                    mp_bitcnt_t bit_size, int negate)
{
    mp_bitcnt_t bits = 0;
    mp_size_t limbs = 0;
    mp_bitcnt_t b = bit_size % FLINT_BITS;
    mp_size_t l = bit_size / FLINT_BITS;
    int borrow = 0;
    len_t i;

    for (i = 0; i < len; i++)
    {
        borrow =
            fmpz_bit_pack(arr + limbs, bits, bit_size, poly + i, negate,
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
fmpz_poly_bit_pack(fmpz_t f, const fmpz_poly_t poly,
                   mp_bitcnt_t bit_size)
{
    len_t len;
    __mpz_struct * mpz;
    len_t i, d;
    int negate;

    len = fmpz_poly_length(poly);

    if (len == 0 || bit_size == 0)
    {
        fmpz_zero(f);
        return;
    }

    mpz = _fmpz_promote(f);
    mpz_realloc2(mpz, len * bit_size);
    d = mpz->_mp_alloc;

    flint_mpn_zero(mpz->_mp_d, d);

    if (fmpz_sgn(fmpz_poly_lead(poly)) < 0)
        negate = -1;
    else
        negate = 0;

    _fmpz_poly_bit_pack(mpz->_mp_d, poly->coeffs, len, bit_size, negate);

    for (i = d - 1; i >= 0; i--)
    {
        if (mpz->_mp_d[i] != 0)
            break;
    }
    d = i + 1;

    mpz->_mp_size = d;
    _fmpz_demote_val(f);

    if (negate)
        fmpz_neg(f, f);
}
