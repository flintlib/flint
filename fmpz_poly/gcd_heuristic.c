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

    Copyright (C) 2011 William Hart
   
******************************************************************************/

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

int _fmpz_poly_divides(fmpz * q, const fmpz * a, 
                       long len1, const fmpz * b, long len2)
{
   fmpz * r = _fmpz_vec_init(len1);
   
   _fmpz_poly_divrem(q, r, a, len1, b, len2);
   while ((len2) && r[len2 - 1] == 0) len2--;
   
   _fmpz_vec_clear(r, len1);

   return (len2 == 0);
}

mp_size_t mpn_gcd_adjusted(mp_ptr arrayg, 
    mp_ptr array1, mp_size_t limbs1, mp_ptr array2, mp_size_t limbs2)
{
   mp_size_t s1 = 0, s2 = 0, m, b1, b2, mb, len1, len2, leng;
   mp_ptr in1, in2;
   mp_limb_t cy;

   while (array1[s1] == 0) s1++; len1 = limbs1 - s1;
   while (array2[s2] == 0) s2++; len2 = limbs2 - s2;
   m = FLINT_MIN(s1, s2);
   mpn_zero(arrayg, m);

   b1 = mpn_scan1(array1 + s1, 0);
   b2 = mpn_scan1(array2 + s2, 0);
   mb = FLINT_MIN(b1, b2);
   if (s1 > s2) mb = b2;
   if (s2 > s1) mb = b1;

   if (b1 == 0)
      in1 = array1 + s1;
   else
   {
      in1 = malloc(len1*sizeof(mp_limb_t));
      mpn_rshift(in1, array1 + s1, len1, b1);
      len1 -= (in1[len1 - 1] == 0); 
   }

   if (b2 == 0)
      in2 = array2 + s2;
   else
   {
      in2 = malloc(len2*sizeof(mp_limb_t));
      mpn_rshift(in2, array2 + s2, len2, b2);
      len2 -= (in2[len2 - 1] == 0); 
   }
   
   if (len1 >= len2)
      leng = mpn_gcd(arrayg + m, in1, len1, in2, len2);
   else 
      leng = mpn_gcd(arrayg + m, in2, len2, in1, len1);

   if (mb)
   {
      cy = mpn_lshift(arrayg + m, arrayg + m, leng, mb);
      if (cy)
         arrayg[m + leng++] = cy;
   }

   if (b1) free(in1);
   if (b2) free(in2);

   return m + leng;
}

/* Assumes len1 != 0 != len2 */

int
_fmpz_poly_gcd_heuristic(fmpz * res, const fmpz * poly1, long len1, 
                                        const fmpz * poly2, long len2)
{
	ulong bits1, bits2, max_bits, pack_bits, bound_bits, bits_R, bits_Q;
   ulong limbs1, limbs2, limbsg, pack_limbs;
   ulong bound_limbs, qlimbs, tlimbs, hlimbs, rlimbs;
   ulong log1, log_length;
   long sign1, sign2, length, length2;
	fmpz_t ac, bc, d, rc;
   fmpz * A, * B, * R, * Q, * prod;
   mp_ptr array1, array2, arrayg, temp, q, r;
   int divides, ok;

   fmpz_init(ac);
   fmpz_init(bc);
   fmpz_init(d);
   
	_fmpz_poly_content(ac, poly1, len1);
   _fmpz_poly_content(bc, poly2, len2);

   fmpz_gcd(d, ac, bc);

   if (len1 == 1 || len2 == 1)
   {
      fmpz_set(res, d);

      fmpz_clear(ac);
      fmpz_clear(bc);
	   fmpz_clear(d);

      return 1;
   }
   
   A = _fmpz_vec_init(len1);
   B = _fmpz_vec_init(len2);
   _fmpz_vec_scalar_divexact_fmpz(A, poly1, len1, ac);
   _fmpz_vec_scalar_divexact_fmpz(B, poly2, len2, bc);
   fmpz_clear(ac);
   fmpz_clear(bc);
	   
	if (len2 == 2) /* if len1 == 2 then so does len2 */
	{
		Q = _fmpz_vec_init(len1 - len2 + 1);
		if (_fmpz_poly_divides(Q, A, len1, B, len2))
        {
		   _fmpz_vec_scalar_mul_fmpz(res, B, len2, d);
            if (fmpz_sgn(res + 1) < 0)
            {
                fmpz_neg(res, res);
                fmpz_neg(res + 1, res + 1);
            }
        }
		else  
        {
			fmpz_set(res, d);
            fmpz_zero(res + 1);
        }

		fmpz_clear(d);
		_fmpz_vec_clear(A, len1);
      _fmpz_vec_clear(B, len2);
      _fmpz_vec_clear(Q, len1 - len2 + 1);
      
      return 1;
	}
	
   R = _fmpz_vec_init(FLINT_MIN(len1, len2));
   
	bits1 = FLINT_ABS(_fmpz_vec_max_bits(A, len1));
	bits2 = FLINT_ABS(_fmpz_vec_max_bits(B, len2));
	max_bits = FLINT_MAX(bits1, bits2);

	sign1 = (long) fmpz_sgn(A + len1 - 1);
	sign2 = (long) fmpz_sgn(B + len2 - 1);
   
	/*
	   This bound ensures that if R | A and R | B with R primitive then R is the 
		gcd of A and B. The bound is taken from 
		http://arxiv.org/abs/cs/0206032v1
	*/
			
	bound_bits = FLINT_MIN(bits1, bits2) + 6; 
	bound_limbs = (bound_bits - 1)/FLINT_BITS + 1; 
	
	pack_limbs = (max_bits - 1)/FLINT_BITS + 1;
   if (bound_limbs > pack_limbs) pack_limbs = bound_limbs;

	array1 = calloc(len1*pack_limbs, sizeof(mp_limb_t));
   array2 = calloc(len2*pack_limbs, sizeof(mp_limb_t));
   arrayg = calloc(len2*pack_limbs, sizeof(mp_limb_t));
	
#if FLINT64
	if ((bits1 < 32) && (bits2 < 32) && (bound_bits < 32))
	{
	   pack_bits = FLINT_MAX(bits1, bits2) + 6;
		if (pack_bits > 32) pack_bits = 32;
   } else
#endif
   if ((bits1 < FLINT_BITS) && (bits2 < FLINT_BITS) && (bound_bits < FLINT_BITS))
		pack_bits = FLINT_BITS;
	else
		pack_bits = FLINT_BITS*pack_limbs;
		
   _fmpz_poly_bit_pack(array1, A, len1, pack_bits, sign1);
	limbs1 = (pack_bits*len1 - 1)/FLINT_BITS + 1;
   while (!array1[limbs1 - 1]) limbs1--;

   _fmpz_poly_bit_pack(array2, B, len2, pack_bits, sign2);
	limbs2 = (pack_bits*len2 - 1)/FLINT_BITS + 1;
	while (!array2[limbs2 - 1]) limbs2--;
	
	limbsg = mpn_gcd_adjusted(arrayg, array1, limbs1, array2, limbs2);
	
   /* may have one extra coeff due to 1 0 -x being packed as 0 -1 -x */
   length = (limbsg*FLINT_BITS)/pack_bits + 1; 
   length = FLINT_MIN(length, len2);
   
   R = _fmpz_vec_init(length);
   
   mpn_zero(res, len2);
   _fmpz_poly_bit_unpack(R, length, arrayg, pack_bits, 0);
   while (R[length - 1] == 0) length--;

   Q = _fmpz_vec_init(len1); 
   
	fmpz_init(rc);
	_fmpz_poly_content(rc, R, length);

   if (!fmpz_is_one(rc)) 
   {
      temp = malloc(limbsg*sizeof(mp_limb_t));

      if (fmpz_size(rc) == 1 || -fmpz_size(rc) == 1) 
      {
         mpn_divmod_1(temp, arrayg, limbsg, fmpz_get_ui(rc));
         tlimbs = limbsg - (temp[limbsg - 1] == 0);
      }
	   else 
      {
         __mpz_struct * mpz_ptr = COEFF_TO_PTR(*rc);
         hlimbs = FLINT_ABS(mpz_ptr->_mp_size);
         mpn_tdiv_q(temp, arrayg, limbsg, mpz_ptr->_mp_d, hlimbs);
         tlimbs = limbsg - hlimbs + 1;
         tlimbs -= (temp[tlimbs - 1] == 0);
      }
   } else 
   {
      temp = arrayg;
      tlimbs = limbsg;
   }

   qlimbs = FLINT_MAX(limbs1, limbs2) + pack_limbs + 1;
	q = malloc(qlimbs*sizeof(mp_limb_t));
   r = malloc(tlimbs*sizeof(mp_limb_t));

	divides = 0;

	mpn_zero(q, qlimbs);
   mpn_tdiv_qr(q, r, 0, array1, limbs1, temp, tlimbs);
   rlimbs = tlimbs;
   while ((rlimbs) && r[rlimbs - 1] == 0) rlimbs--;
   
   if (rlimbs == 0) /* division was exact */
	{
      length2 = (qlimbs*FLINT_BITS)/pack_bits + 1;
      length2 = FLINT_MIN(length2, len1);
      _fmpz_poly_bit_unpack(Q, length2, q, pack_bits, 0);
      while (Q[length2 - 1] == 0) length2--;
      
      bits_R = FLINT_ABS(_fmpz_vec_max_bits(R, length));
		bits_Q = FLINT_ABS(_fmpz_vec_max_bits(Q, length2));
		log1 = FLINT_BIT_COUNT(length);
		log_length = FLINT_MIN(log1, FLINT_BIT_COUNT(length2));
       
	   ok = 0;
		if (bits_R + bits_Q + log_length < pack_bits)
			ok = 1;
		else
		{
			prod = _fmpz_vec_init(length + length2 - 1);
         _fmpz_vec_scalar_divexact_fmpz(R, R, length, rc);
		   if (length2 >= length)
            _fmpz_poly_mul(prod, Q, length2, R, length);
		   else
            _fmpz_poly_mul(prod, R, length, Q, length2);
         if (sign1 < 0L) _fmpz_vec_neg(prod, prod, length + length2 - 1);
		}
      
		if (ok || (length + length2 - 1 == len1 && _fmpz_vec_equal(Q, A, len1)))
		{
         if (!ok) _fmpz_vec_clear(prod, length + length2 - 1);

         mpn_zero(q, qlimbs);
         mpn_tdiv_qr(q, r, 0, array2, limbs2, temp, tlimbs);
         rlimbs = tlimbs;
         while ((rlimbs) && r[rlimbs - 1] == 0) rlimbs--;
         
         if (rlimbs == 0) /* division was exact */
	      {
            if (ok) _fmpz_vec_scalar_divexact_fmpz(R, R, length, rc);
		      
            _fmpz_poly_bit_unpack(Q, length2, q, pack_bits, 0);
            while (Q[length2 - 1] == 0) length2--;
            
            bits_Q = FLINT_ABS(_fmpz_vec_max_bits(Q, length2));
				log_length = FLINT_MIN(log1, FLINT_BIT_COUNT(length2));

				ok = 0;
				if (bits_R + bits_Q + log_length < pack_bits)
		      {
			      ok = 1;
		      } else
				{
		         prod = _fmpz_vec_init(length + length2 - 1);
               if (length2 >= length)
                  _fmpz_poly_mul(prod, Q, length2, R, length);
		         else
                  _fmpz_poly_mul(prod, R, length, Q, length2);
               if (sign2 < 0L) _fmpz_vec_neg(prod, prod, length + length2 - 1);
				}
				if (ok || (length + length2 - 1 == len2 && _fmpz_vec_equal(Q, B, len2))) 
               divides = 1;
            if (!ok) _fmpz_vec_clear(prod, length + length2 - 1);
			} 
		} 
	}

   free(q); 
	free(r); 
	free(arrayg); 
	free(array1); 
	free(array2); 
	fmpz_clear(rc); 
	if (temp != arrayg) 
      free(temp); 

	_fmpz_vec_clear(Q, len1);
	_fmpz_vec_clear(A, len1);
	_fmpz_vec_clear(B, len2);
	
   if (divides)
	{
		_fmpz_vec_scalar_mul_fmpz(res, R, length, d);
		fmpz_clear(d);
      _fmpz_vec_clear(R, length);
		
      return 1;
	} else
	{
	   fmpz_clear(d);
	   _fmpz_vec_clear(R, length);
		
      return 0;
	}
}

int
fmpz_poly_gcd_heuristic(fmpz_poly_t res,
                           const fmpz_poly_t poly1, const fmpz_poly_t poly2)
{
    const long len1 = poly1->length;
    const long len2 = poly2->length;
    long rlen;
    int done = 0;

    if (len1 == 0)
    {
        if (len2 == 0)
            fmpz_poly_zero(res);
        else
        {
            if (fmpz_sgn(poly2->coeffs + (len2 - 1)) > 0)
                fmpz_poly_set(res, poly2);
            else
                fmpz_poly_neg(res, poly2);
        }
        return 1;
    }
    else
    {
        if (len2 == 0)
        {
            if (fmpz_sgn(poly1->coeffs + (len1 - 1)) > 0)
                fmpz_poly_set(res, poly1);
            else
                fmpz_poly_neg(res, poly1);
            return 1;
        }
    }

    rlen = FLINT_MIN(len1, len2);

    if (res == poly1 || res == poly2)
    {
       fmpz_poly_t temp;
       fmpz_poly_init2(temp, rlen);
       if (len1 >= len2)
          done = _fmpz_poly_gcd_heuristic(temp->coeffs, poly1->coeffs, len1,
                                    poly2->coeffs, len2);
       else
          done = _fmpz_poly_gcd_heuristic(temp->coeffs, poly2->coeffs, len2,
                                    poly1->coeffs, len1);
       fmpz_poly_swap(temp, res);
       fmpz_poly_clear(temp);
    }
    else
    {
       fmpz_poly_fit_length(res, rlen);
       if (len1 >= len2)
          done = _fmpz_poly_gcd_heuristic(res->coeffs, poly1->coeffs, len1,
                                    poly2->coeffs, len2);
       else
          done = _fmpz_poly_gcd_heuristic(res->coeffs, poly2->coeffs, len2,
                                    poly1->coeffs, len1);
    }
    
    if (done)
    {
       _fmpz_poly_set_length(res, rlen);
       _fmpz_poly_normalise(res);
    }

    return done;
}
