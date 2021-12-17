/*
    Copyright (C) 2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "mpn_extras.h"

/* 
   Divide (arrayg, limbsg) by the positive value gc in-place and
   return the number of limbs written
*/
mp_size_t flint_mpn_tdiv_q_fmpz_inplace(mp_ptr arrayg, mp_size_t limbsg, fmpz_t gc)
{
   if (fmpz_size(gc) == 1) 
   {
      mpn_divmod_1(arrayg, arrayg, limbsg, fmpz_get_ui(gc));
      return limbsg - (arrayg[limbsg - 1] == 0);
   }
	else 
   {
      mp_size_t tlimbs;
      __mpz_struct * mgc = COEFF_TO_PTR(*gc);
      
      mp_ptr temp = flint_malloc(limbsg*sizeof(mp_limb_t));
      flint_mpn_copyi(temp, arrayg, limbsg);
      
      mpn_tdiv_q(arrayg, temp, limbsg, mgc->_mp_d, mgc->_mp_size);
      tlimbs = limbsg - mgc->_mp_size + 1;
      tlimbs -= (arrayg[tlimbs - 1] == 0);
      
      flint_free(temp);
      return tlimbs;
   } 
}

/* 
   Returns 1 if sign * (G, glen) * (Q, qlen) = (P, len), else returns 0.
   Temp requires space for glen + qlen - 1 coefficients
*/
int multiplies_out(fmpz * P, slong len, const fmpz * Q, slong qlen, 
                   const fmpz * G, slong glen, slong sign, fmpz * temp)
{
   int divides = 0;

   /* multiply out */
   if (qlen >= glen)
      _fmpz_poly_mul(temp, Q, qlen, G, glen);
   else
      _fmpz_poly_mul(temp, G, glen, Q, qlen);
   if (sign < WORD(0)) _fmpz_vec_neg(temp, temp, glen + qlen - 1);

   /* check if quotient really was exact */
   divides = (glen + qlen - 1 == len && _fmpz_vec_equal(temp, P, len));

   return divides;
}
   
/* Assumes len1 != 0 != len2 */
int
_fmpz_poly_gcd_heuristic(fmpz * res, const fmpz * poly1, slong len1, 
                                        const fmpz * poly2, slong len2)
{
   slong bits1, bits2, max_bits, pack_bits, bound_bits, bits_G, bits_Q;
   ulong limbs1, limbs2, limbsg, pack_limbs, qlimbs, qlimbs2;
   ulong log_glen, log_length;
   slong sign1, sign2, glen, qlen, qlen2;
	fmpz_t ac, bc, d, gc;
   fmpz * A, * B, * G, * Q, * t;
   mp_ptr array1, array2, arrayg, q, temp;
   int divides;

   fmpz_init(ac);
   fmpz_init(bc);
   fmpz_init(d);
   
	/* compute gcd of content of poly1 and poly2 */
   _fmpz_poly_content(ac, poly1, len1);
   _fmpz_poly_content(bc, poly2, len2);
   fmpz_gcd(d, ac, bc);

   /* special case, one of the polys is a constant */
   if (len2 == 1) /* if len1 == 1 then so does len2 */
   {
      fmpz_set(res, d);

      fmpz_clear(ac);
      fmpz_clear(bc);
	   fmpz_clear(d);

      return 1;
   }
   
   /* divide poly1 and poly2 by their content */
   A = _fmpz_vec_init(len1);
   B = _fmpz_vec_init(len2);
   _fmpz_vec_scalar_divexact_fmpz(A, poly1, len1, ac);
   _fmpz_vec_scalar_divexact_fmpz(B, poly2, len2, bc);
   fmpz_clear(ac);
   fmpz_clear(bc);
	   
	/* special case, one of the polys is length 2 */
   if (len2 == 2) /* if len1 == 2 then so does len2 */
	{
		Q = _fmpz_vec_init(len1 - len2 + 1);
		if (_fmpz_poly_divides(Q, A, len1, B, 2))
      {
		   _fmpz_vec_scalar_mul_fmpz(res, B, 2, d);
         if (fmpz_sgn(res + 1) < 0)
            _fmpz_vec_neg(res, res, 2);
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
	
   /* 
      Determine how many bits (pack_bits) to pack into. The bound 
      bound_bits ensures that if G | A and G | B with G primitive 
      then G is the gcd of A and B. The bound is taken from 
      https://arxiv.org/abs/cs/0206032v1
   */
   bits1 = FLINT_ABS(_fmpz_vec_max_bits(A, len1));
   bits2 = FLINT_ABS(_fmpz_vec_max_bits(B, len2));
   /* 
      always extra bit for signs whether polys are signed or not, since we don't
      know if any purported gcds/quotients will be signed
   */
   max_bits = FLINT_MAX(bits1, bits2) + 1; 
   			
	/* 
	   the +6 is chosen heuristically for performance; the theorem
	   is satisfied with +3 (including a bit for signs)
	*/
	bound_bits = FLINT_MIN(bits1, bits2) + 6; 
	pack_bits = FLINT_MAX(bound_bits, max_bits); /* need to pack original polys */
   pack_limbs = (pack_bits - 1)/FLINT_BITS + 1;
   
	if (pack_bits >= 32) /* pack into multiples of limbs if >= 32 bits */
      pack_bits = FLINT_BITS*pack_limbs;
		
   /* allocate space to pack into */
   limbs1 = (pack_bits*len1 - 1)/FLINT_BITS + 1;
   limbs2 = (pack_bits*len2 - 1)/FLINT_BITS + 1;
	array1 = flint_calloc(limbs1, sizeof(mp_limb_t));
   array2 = flint_calloc(limbs2, sizeof(mp_limb_t));
   arrayg = flint_calloc(limbs2, sizeof(mp_limb_t));
   
   /* pack first poly and normalise */
   sign1 = (slong) fmpz_sgn(A + len1 - 1);
	_fmpz_poly_bit_pack(array1, A, len1, pack_bits, sign1);
	while (array1[limbs1 - 1] == 0) limbs1--;

   /* pack second poly and normalise */
   sign2 = (slong) fmpz_sgn(B + len2 - 1);
   _fmpz_poly_bit_pack(array2, B, len2, pack_bits, sign2);
	while (array2[limbs2 - 1] == 0) limbs2--;
	
	/* compute integer GCD */
   limbsg = flint_mpn_gcd_full(arrayg, array1, limbs1, array2, limbs2);
	
   /* 
      Make space for unpacked gcd. May have one extra coeff due to 
      1 0 -x being packed as 0 -1 -x. 
   */
   glen = FLINT_MIN((limbsg*FLINT_BITS)/pack_bits + 1, len2); 
   G = _fmpz_vec_init(glen);

   /* 
      clear bits after g in arrayg so they are not inadvertently
      pulled into G after bit unpacking
   */
   flint_mpn_zero(arrayg + limbsg, limbs2-limbsg);

   /* unpack gcd */
   _fmpz_poly_bit_unpack(G, glen, arrayg, pack_bits, 0);
   while (G[glen - 1] == 0) glen--;
   
	/* divide by any content */
   fmpz_init(gc);
	_fmpz_poly_content(gc, G, glen);

   if (!fmpz_is_one(gc)) 
      limbsg = flint_mpn_tdiv_q_fmpz_inplace(arrayg, limbsg, gc);

   /* make space for quotient and remainder of both polys by gcd */
   qlimbs = limbs1 - limbsg + 1;
   qlen = FLINT_MIN(len1, (qlimbs*FLINT_BITS)/pack_bits + 1);
   qlimbs2 = limbs2 - limbsg + 1;
   qlen2 = FLINT_MIN(len2, (qlimbs2*FLINT_BITS)/pack_bits + 1);
   qlimbs = (FLINT_MAX(qlen, qlen2)*pack_bits - 1)/FLINT_BITS + 1;
   q = flint_calloc(qlimbs, sizeof(mp_limb_t));
   temp = flint_malloc(limbsg*sizeof(mp_limb_t));
	divides = 0;

   if (flint_mpn_divides(q, array1, limbs1, arrayg, limbsg, temp)) 
	{
      /* unpack quotient of first poly by gcd */
      Q = _fmpz_vec_init(len1); 
      t = _fmpz_vec_init(len1 + glen);
      _fmpz_poly_bit_unpack(Q, qlen, q, pack_bits, 0);
      while (Q[qlen - 1] == 0) qlen--;
            
      /* divide by content */
      _fmpz_vec_scalar_divexact_fmpz(G, G, glen, gc);
		
      /* check if we really need to multiply out to check for exact quotient */
      bits_G = FLINT_ABS(_fmpz_vec_max_bits(G, glen));
	  bits_Q = FLINT_ABS(_fmpz_vec_max_bits(Q, qlen));
	  log_glen = FLINT_BIT_COUNT(glen);
	  log_length = FLINT_MIN(log_glen, FLINT_BIT_COUNT(qlen));
       
	  /* allow one bit for signs */
	  divides = (bits_G + bits_Q + log_length < pack_bits);
     
      if (!divides) /* need to multiply out to check exact quotient */
         divides = multiplies_out(A, len1, Q, qlen, G, glen, sign1, t);

		if (divides) /* quotient really was exact */
		{
           divides = 0;
           flint_mpn_zero(q, qlimbs);
          
         if (flint_mpn_divides(q, array2, limbs2, arrayg, limbsg, temp)) 
	      {
            /* unpack quotient of second poly by gcd */
            _fmpz_poly_bit_unpack(Q, qlen2, q, pack_bits, 0);
            while (Q[qlen2 - 1] == 0) qlen2--;
			
			/* check if we really need to multiply out to check for exact quotient */
            bits_Q = FLINT_ABS(_fmpz_vec_max_bits(Q, qlen2));
		    log_length = FLINT_MIN(log_glen, FLINT_BIT_COUNT(qlen2));

			/* allow one bit for signs */
			divides = (bits_G + bits_Q + log_length < pack_bits);
       
            if (!divides) /* we need to multiply out */
               divides = multiplies_out(B, len2, Q, qlen2, G, glen, sign1, t);
			} 
		} 

      _fmpz_vec_clear(t, len1 + glen);
      _fmpz_vec_clear(Q, len1);
	}

   flint_free(q); 
	flint_free(temp); 
	flint_free(arrayg); 
	flint_free(array1); 
	flint_free(array2); 
	fmpz_clear(gc); 
	
	_fmpz_vec_clear(A, len1);
	_fmpz_vec_clear(B, len2);
	
   /* we found the gcd, so multiply by content */
   if (divides)
   {
	   _fmpz_vec_zero(res + glen, len2 - glen);
      _fmpz_vec_scalar_mul_fmpz(res, G, glen, d);
   }
		
   fmpz_clear(d);
   _fmpz_vec_clear(G, glen);
		
   return divides;
}

int
fmpz_poly_gcd_heuristic(fmpz_poly_t res, const fmpz_poly_t poly1,
              const fmpz_poly_t poly2)
{
    if (poly1->length < poly2->length)
    {
        return fmpz_poly_gcd_heuristic(res, poly2, poly1);
    }
    else /* len1 >= len2 >= 0 */
    {
        const slong len1 = poly1->length;
        const slong len2 = poly2->length;
        int done = 1; /* len1 = 0 or len2 = 0 need done = 1 */

        if (len1 == 0) /* len1 = len2 = 0 */
        {
            fmpz_poly_zero(res);
        } 
        else if (len2 == 0) /* len1 > len2 = 0 */
        {
            if (fmpz_sgn(poly1->coeffs + (len1 - 1)) > 0)
                fmpz_poly_set(res, poly1);
            else
                fmpz_poly_neg(res, poly1);
        }
        else /* len1 >= len2 >= 1 */
        {
            /* underscore function automatically takes care of aliasing */
           
            fmpz_poly_fit_length(res, len2);
                
            done = _fmpz_poly_gcd_heuristic(res->coeffs, poly1->coeffs, len1,
                                    poly2->coeffs, len2);
            
            if (done)
            {
                _fmpz_poly_set_length(res, len2);
                _fmpz_poly_normalise(res);
            }
        }

        return done;
    }
}
