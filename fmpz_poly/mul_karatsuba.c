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

ulong revbin(ulong in, ulong bits)
{
   ulong out = 0, i;
   
   for (i = 0; i < bits; i++)
   {   
      out <<= 1;
      out += (in & 1);
      in >>= 1;
   }

   return out;
}

void revbin1(fmpz * out, const fmpz * in, ulong len, ulong bits)
{
   ulong i;
 
   for (i = 0; i < len; i++)
      out[revbin(i, bits)] = in[i];
}

void revbin2(fmpz * out, const fmpz * in, ulong len, ulong bits)
{
   ulong i;
 
   for (i = 0; i < len; i++)
      out[i] = in[revbin(i, bits)];
}

void _fmpz_vec_add_rev(fmpz * in1, fmpz * in2, ulong bits)
{
   ulong i, j;
   
   for (i = 0; i < (1L<<bits) - 1; i++)
   {
      j = revbin(revbin(i, bits) + 1, bits);
      fmpz_add(in1 + j, in1 + j, in2 + i);
   }
}

// Assumes rev2 immediately follows rev1
void _fmpz_poly_mul_kara_recursive(fmpz * out, fmpz * rev1, fmpz * rev2, ulong length, fmpz * temp, ulong bits)
{
   ulong m = length/2, i;
   
   if (length == 1)
   {
      fmpz_mul(out, rev1, rev2);
      fmpz_zero(out + 1);
      return;
   }

   _fmpz_vec_add(temp, rev1, rev1 + m, m);
   _fmpz_vec_add(temp + m, rev2, rev2 + m, m);

   _fmpz_poly_mul_kara_recursive(out, rev1, rev2, m, temp + 2*m, bits - 1);
   
   _fmpz_poly_mul_kara_recursive(out + length, temp, temp + m, m, temp + 2*m, bits - 1);

   _fmpz_poly_mul_kara_recursive(temp, rev1 + m, rev2 + m, m, temp + 2*m, bits - 1);

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

   _fmpz_poly_mul_kara_recursive(out, rev1, rev2, length, temp, loglen);  

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
