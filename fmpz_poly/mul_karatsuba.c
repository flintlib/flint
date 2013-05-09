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

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

/*
   Implements karatsuba multiplication. There is no basecase crossover, so 
   this is only efficient when the coefficients are large (the main usage 
   case).

   The algorithm is the "odd/even" Karatsuba algorithm. Let 
   f(x) = f1(x^2) + x*f2(x^2), g(x) = g1(x^2) + x*g2(x^2), then 
   
   f(x)*g(x) = f1(x^2)*g1(x^2) + x^2*f2(x^2)*g2(x^2) 
               + x*((f1(x^2) + f2(x^2))*(g1(x^2) + g2(x^2)) 
                    - f1(x^2)*g1(x^2) - f2(x^2)*g2(x^2)).
   
   Thus only three multiplications are performed (and numerous additions 
   and subtractions).
   
   Instead of working with polynomials with the usual ordering, reverse 
   binary ordering is used, i.e. for length 2^3 (zero padded) terms of 
   degree 110 and 011 in binary are swapped, etc.

   The advantage of working in this format is that the first half of the 
   coefficients of f will be the coefficients of f1, and the second half, 
   those of f2, etc. This applies right down the recursion. The only tricky 
   bit is when multiplying by x. One must undo the revbin to shift by one 
   term to the left.
*/

void _fmpz_poly_mul_kara_recursive(fmpz * out, fmpz * rev1, fmpz * rev2,
                                   fmpz * temp, len_t bits);

/*
   Switches the coefficients of poly in of length len into a 
   poly out of length 2^bits.
 */
void
revbin1(fmpz * out, const fmpz * in, len_t len, len_t bits)
{
    len_t i;
    for (i = 0; i < len; i++)
        out[n_revbin(i, bits)] = in[i];
}

/*
   Switches the coefficients of poly in of length 2^bits into a
   poly out of length len.
 */
void
revbin2(fmpz * out, const fmpz * in, len_t len, len_t bits)
{
    len_t i;
    for (i = 0; i < len; i++)
        out[i] = in[n_revbin(i, bits)];
}

/* in1 += x*in2 assuming both in1 and in2 are revbin'd. */
void
_fmpz_vec_add_rev(fmpz * in1, fmpz * in2, len_t bits)
{
    len_t i;
    for (i = 0; i < (1L << bits) - 1; i++)
    {
        len_t j = n_revbin(n_revbin(i, bits) + 1, bits);
        fmpz_add(in1 + j, in1 + j, in2 + i);
    }
}

/*
   Recursive Karatsuba assuming polynomials are in revbin format.
   
   Assumes rev1 and rev2 are both of length 2^bits and that temp has 
   space for 2^bits coefficients.
 */
void
_fmpz_poly_mul_kara_recursive(fmpz * out, fmpz * rev1, fmpz * rev2,
                              fmpz * temp, len_t bits)
{
    len_t length = (1L << bits);
    len_t m = length / 2;

    if (length == 1)
    {
        fmpz_mul(out, rev1, rev2);
        fmpz_zero(out + 1);
        return;
    }

    _fmpz_vec_add(temp, rev1, rev1 + m, m);
    _fmpz_vec_add(temp + m, rev2, rev2 + m, m);

    _fmpz_poly_mul_kara_recursive(out, rev1, rev2, temp + 2 * m, bits - 1);

    _fmpz_poly_mul_kara_recursive(out + length, temp, temp + m, temp + 2 * m,
                                  bits - 1);

    _fmpz_poly_mul_kara_recursive(temp, rev1 + m, rev2 + m, temp + 2 * m,
                                  bits - 1);

    _fmpz_vec_sub(out + length, out + length, out, length);
    _fmpz_vec_sub(out + length, out + length, temp, length);

    _fmpz_vec_add_rev(out, temp, bits);
}

/* Assumes poly1 and poly2 are not length 0 and len1 >= len2. */
void
_fmpz_poly_mul_karatsuba(fmpz * res, const fmpz * poly1,
                         len_t len1, const fmpz * poly2, len_t len2)
{
    fmpz *rev1, *rev2, *out, *temp;
    len_t length, loglen = 0;

    if (len1 == 1)
    {
        fmpz_mul(res, poly1, poly2);
        return;
    }

    while ((1L << loglen) < len1)
        loglen++;
    length = (1L << loglen);

    rev1 = (fmpz *) flint_calloc(4 * length, sizeof(fmpz *));
    rev2 = rev1 + length;
    out  = rev1 + 2 * length;
    temp = _fmpz_vec_init(2 * length);

    revbin1(rev1, poly1, len1, loglen);
    revbin1(rev2, poly2, len2, loglen);

    _fmpz_poly_mul_kara_recursive(out, rev1, rev2, temp, loglen);

    _fmpz_vec_zero(res, len1 + len2 - 1);
    revbin2(res, out, len1 + len2 - 1, loglen + 1);

    _fmpz_vec_clear(temp, 2 * length);
    flint_free(rev1);
}

void
fmpz_poly_mul_karatsuba(fmpz_poly_t res,
                        const fmpz_poly_t poly1, const fmpz_poly_t poly2)
{
    len_t len_out;

    if ((poly1->length == 0) || (poly2->length == 0))
    {
        fmpz_poly_zero(res);
        return;
    }

    len_out = poly1->length + poly2->length - 1;

    fmpz_poly_fit_length(res, len_out);

    if (poly1->length >= poly2->length)
        _fmpz_poly_mul_karatsuba(res->coeffs, poly1->coeffs, poly1->length,
                                 poly2->coeffs, poly2->length);
    else
        _fmpz_poly_mul_karatsuba(res->coeffs, poly2->coeffs, poly2->length,
                                 poly1->coeffs, poly1->length);

    _fmpz_poly_set_length(res, len_out);
}
