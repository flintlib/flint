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
/****************************************************************************

   Copyright (C) 2010 William Hart
   
*****************************************************************************/

#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

/**
   Implements karatsuba multiplication. There is no basecase crossover, so this is
   only efficient when the coefficients are large (the main usage case).

   The algorithm is the "odd/even" Karatsuba algorithm. Let f(x) = f1(x^2) + x*f2(x^2),
   g(x) = g1(x^2) + x*g2(x^2), then f(x)*g(x) = f1(x^2)*g1(x^2) + x^2*f2(x^2)*g2(x^2)
   + x*((f1(x^2) + f2(x^2))*(g1(x^2) + g2(x^2)) - f1(x^2)*g1(x^2) - f2(x^2)*g2(x^2)).

   Thus only three multiplications are performed (and numerous additions and subtractions).
   
   Instead of working with polynomials with the usual ordering, reverse binary ordering 
   is used, i.e. for length 2^3 (zero padded) terms of degree 110 and 011 in binary are
   swapped, etc.

   The advantage of working in this format is that the first half of the coefficients of f
   will be the coefficients of f1, and the second half, those of f2, etc. This applies 
   right down the recursion. The only tricky bit is when multiplying by x. One must undo 
   the revbin to shift by one term to the left.
*/

const ulong revtab1[2] = { 0, 1 };
const ulong revtab2[4] = { 0, 2, 1, 3 };
const ulong revtab3[8] = { 0, 4, 2, 6, 1, 5, 3, 7 };
const ulong revtab4[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };

const ulong * revtab[5] = { NULL, revtab1, revtab2, revtab3, revtab4 };

// computes the reverse binary of a binary number of the given number of bits
ulong revbin(ulong in, ulong bits)
{
   ulong out = 0, i;
   
   if (bits < 3)
      return revtab[bits][in];

   for (i = 0; i < bits; i++)
   {   
      out <<= 1;
      out += (in & 1);
      in >>= 1;
   }

   return out;
}

// switches the coefficients of poly in of length len into a
// poly out of length 2^bits
void revbin1(fmpz * out, const fmpz * in, ulong len, ulong bits)
{
   ulong i;
 
   for (i = 0; i < len; i++)
      out[revbin(i, bits)] = in[i];
}

// switches the coefficients of poly in of length 2^bits into a
// poly out of length len
void revbin2(fmpz * out, const fmpz * in, ulong len, ulong bits)
{
   ulong i;
 
   for (i = 0; i < len; i++)
      out[i] = in[revbin(i, bits)];
}

// in1 += x*in2 assuming both in1 and in2 are revbin'd
void _fmpz_vec_add_rev(fmpz * in1, fmpz * in2, ulong bits)
{
   ulong i, j;
   
   for (i = 0; i < (1L<<bits) - 1; i++)
   {
      j = revbin(revbin(i, bits) + 1, bits);
      fmpz_add(in1 + j, in1 + j, in2 + i);
   }
}

// recursive Karatsuba assuming polynomials are in revbin format
// Assumes rev1 and rev2 are both of length 2^bits and that temp has 
// space for 2^bits coefficients
void _fmpz_poly_mul_kara_recursive(fmpz * out, fmpz * rev1, fmpz * rev2, fmpz * temp, ulong bits)
{
   ulong length = (1L<<bits);
   ulong m = length/2, i;
   
   if (length == 1)
   {
      fmpz_mul(out, rev1, rev2);
      fmpz_zero(out + 1);
      return;
   }

   _fmpz_vec_add(temp, rev1, rev1 + m, m);
   _fmpz_vec_add(temp + m, rev2, rev2 + m, m);

   _fmpz_poly_mul_kara_recursive(out, rev1, rev2, temp + 2*m, bits - 1);
   
   _fmpz_poly_mul_kara_recursive(out + length, temp, temp + m, temp + 2*m, bits - 1);

   _fmpz_poly_mul_kara_recursive(temp, rev1 + m, rev2 + m, temp + 2*m, bits - 1);

   _fmpz_vec_sub(out + length, out + length, out, length); 
   _fmpz_vec_sub(out + length, out + length, temp, length); 

   _fmpz_vec_add_rev(out, temp, bits);
}

// Assumes poly1 and poly2 are not length 0 and len1 >= len2
void _fmpz_poly_mul_karatsuba(fmpz * res, const fmpz * poly1, 
							    ulong len1, const fmpz * poly2, ulong len2)
{
   fmpz * rev1, * rev2, * out, * temp;
   ulong loglen = 0, length;

   while ((1L<<loglen) < len1) loglen++;
   length = (1L<<loglen);

   rev1 = (fmpz *) calloc(2*length, sizeof(fmpz *)); 
   out = (fmpz *) calloc(2*length, sizeof(fmpz *)); 
   temp = (fmpz *) calloc(2*length, sizeof(fmpz *)); 
   
   rev2 = rev1 + length;
   
   revbin1(rev1, poly1, len1, loglen);
   revbin1(rev2, poly2, len2, loglen);

   _fmpz_poly_mul_kara_recursive(out, rev1, rev2, temp, loglen);  

   revbin2(res, out, len1 + len2 - 1, loglen + 1);
   
   _fmpz_vec_clear(temp, 2*length);
   free(out);
   free(rev1);
}

void fmpz_poly_mul_karatsuba(fmpz_poly_t res, 
                         const fmpz_poly_t poly1, const fmpz_poly_t poly2)
{
   ulong len_out;
   
   if ((poly1->length == 0) || (poly2->length == 0))  
   {
      fmpz_poly_zero(res);
      return;
   }

   len_out = poly1->length + poly2->length - 1;

   fmpz_poly_fit_length(res, len_out);

   if (poly1->length >= poly2->length)
      _fmpz_poly_mul_karatsuba(res->coeffs, poly1->coeffs, poly1->length, poly2->coeffs, poly2->length);
   else
      _fmpz_poly_mul_karatsuba(res->coeffs, poly2->coeffs, poly2->length, poly1->coeffs, poly1->length); 
   
   _fmpz_poly_set_length(res, len_out);
   _fmpz_poly_normalise(res);
}
