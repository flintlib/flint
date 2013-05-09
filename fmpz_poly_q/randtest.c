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

    Copyright (C) 2010, 2011 Sebastian Pancratz
   
******************************************************************************/

#include "fmpz_poly_q.h"

void fmpz_poly_q_randtest(fmpz_poly_q_t poly, flint_rand_t state,
                          len_t len1, mp_bitcnt_t bits1, 
                          len_t len2, mp_bitcnt_t bits2)
{
    len2  = FLINT_MAX(len2, 1);
    bits2 = FLINT_MAX(bits2, 1);

    fmpz_poly_randtest(poly->num, state, len1, bits1);
    fmpz_poly_randtest_not_zero(poly->den, state, len2, bits2);
    fmpz_poly_q_canonicalise(poly);
}

void fmpz_poly_q_randtest_not_zero(fmpz_poly_q_t poly, flint_rand_t state, 
                                   len_t len1, mp_bitcnt_t bits1, 
                                   len_t len2, mp_bitcnt_t bits2)
{
    len1  = FLINT_MAX(len1, 1);
    len2  = FLINT_MAX(len2, 1);
    bits1 = FLINT_MAX(bits1, 1);
    bits2 = FLINT_MAX(bits2, 1);

    fmpz_poly_randtest_not_zero(poly->num, state, len1, bits1);
    fmpz_poly_randtest_not_zero(poly->den, state, len2, bits2);
    fmpz_poly_q_canonicalise(poly);
}
