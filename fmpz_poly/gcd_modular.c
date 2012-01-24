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
#include "mpn_extras.h"


void _fmpz_poly_gcd_modular(fmpz * res, const fmpz * poly1, long len1, 
                                        const fmpz * poly2, long len2)
{
   mp_bitcnt_t bits1, bits2, nb1, nb2, gbits, bits_small, pbits, curr_bits = 0, new_bits;   
   fmpz_t ac, bc, hc, d, g, eval_A, eval_B, eval_GCD, modulus;
   fmpz * A, * B, * Q, * lead_A, * lead_B;
   mp_ptr a, b, h;
   mp_limb_t p, h_inv, g_mod;
   nmod_t mod;
   long i, n, hlen, bound;
   int g_pm1;
   
   fmpz_init(ac);
   fmpz_init(bc);
   fmpz_init(d);

   /* compute gcd of content of poly1 and poly2 */
   _fmpz_vec_content(ac, poly1, len1);
   _fmpz_vec_content(bc, poly2, len2);
	fmpz_gcd(d, ac, bc);

   /* special case, one of the polys is a constant */
   if (len2 == 1) /* if len1 == 1 then so does len2 */
   {
      fmpz_set(res, d);

      fmpz_clear(ac);
      fmpz_clear(bc);
	   fmpz_clear(d);

      return;
   }
   
   /* divide poly1 and poly2 by their content */
   A = _fmpz_vec_init(len1);
   B = _fmpz_vec_init(len2);
   _fmpz_vec_scalar_divexact_fmpz(A, poly1, len1, ac);
   _fmpz_vec_scalar_divexact_fmpz(B, poly2, len2, bc);
   fmpz_clear(ac);
   fmpz_clear(bc);
	
	/* get bound on size of gcd coefficients */
   lead_A = A + len1 - 1;
   lead_B = B + len2 - 1;

	bits1 = FLINT_ABS(_fmpz_vec_max_bits(A, len1));
	bits2 = FLINT_ABS(_fmpz_vec_max_bits(B, len2));

	if (len1 < 64 && len2 < 64)
	{
		nb1 = _fmpz_poly_2norm_normalised_bits(A, len1);
		nb2 = _fmpz_poly_2norm_normalised_bits(B, len2);
	} else /* approximate to save time */
	{
		nb1 = (2*bits1 + FLINT_BIT_COUNT(len1) + 1)/2 - fmpz_bits(lead_A) + 1;
		nb2 = (2*bits2 + FLINT_BIT_COUNT(len2) + 1)/2 - fmpz_bits(lead_B) + 1;
	}
   
   /* get gcd of leading coefficients */
   fmpz_init(g);
   fmpz_gcd(g, lead_A, lead_B);
   gbits = fmpz_bits(g);
   g_pm1 = fmpz_is_pm1(g);
   
	/* evaluate -A at -1 */
   fmpz_init(eval_A);
   for (i = 0; i < len1; i++)
	{
      if (i & 1) fmpz_add(eval_A, eval_A, A + i);
		else fmpz_sub(eval_A, eval_A, A + i);
	}

   /* evaluate -B at -1 */
   fmpz_init(eval_B);
   for (i = 0; i < len2; i++)
	{
      if (i & 1) fmpz_add(eval_B, eval_B, B + i);
		else fmpz_sub(eval_B, eval_B, B + i);
	}

   /* compute the gcd of eval(-A, -1) and eval(-B, -1) */
   fmpz_init(eval_GCD);
	fmpz_gcd(eval_GCD, eval_A, eval_B);

   /* compute a heuristic bound after which we should begin checking if we're done */
   bits_small = FLINT_MAX(fmpz_bits(eval_GCD), fmpz_bits(g));
	if (bits_small < 2L) bits_small = 2;

	fmpz_clear(eval_GCD);
	fmpz_clear(eval_A);
	fmpz_clear(eval_B);
   
   /* set size of first prime */
   pbits = FLINT_BITS - 1;
   p = (1UL<<(FLINT_BITS - 1));

   fmpz_init(modulus);
   
   Q = _fmpz_vec_init(len1);

   /* make space for polynomials mod p */
   a = _nmod_vec_init(len1);
   b = _nmod_vec_init(len2);
   h = _nmod_vec_init(len2);

   /* zero entire output */
   _fmpz_vec_zero(res, len2);

   n = len2; /* current bound on length of result */
   bound = len2 + FLINT_MIN(nb1, nb2) + gbits + 1; /* initialise bound */

   for (;;)
   {
      /* get new prime and initialise modulus */
      do { p = n_nextprime(p, 0); }
      while (!fmpz_fdiv_ui(g, p));
      nmod_init(&mod, p);
		
      /* reduce polynomials modulo p */
      _fmpz_vec_get_nmod_vec(a, A, len1, mod);
      _fmpz_vec_get_nmod_vec(b, B, len2, mod);
      
      /* compute gcd over Z/pZ */
      hlen = _nmod_poly_gcd(h, a, len1, b, len2, mod);
      
      if (hlen == 1) /* gcd is 1 */
      {
         fmpz_one(res);
         break; 
      }
      
      if (hlen > n + 1) /* discard this prime */
         continue;      
      
      /* scale new polynomial mod p appropriately */
      if (g_pm1) _nmod_poly_make_monic(h, h, hlen, mod);
      else
      {
         h_inv = n_invmod(h[hlen - 1], mod.n);
         g_mod = fmpz_fdiv_ui(g, mod.n);
         h_inv = n_mulmod2_preinv(h_inv, g_mod, mod.n, mod.ninv);
         _nmod_vec_scalar_mul_nmod(h, h, hlen, h_inv, mod);
      }
      
      if (hlen <= n) /* we have a new bound on size of result */
      {
         _fmpz_vec_set_nmod_vec(res, h, hlen, mod);
         
         if (g_pm1)
         {
            /* are we done? */
            if (_fmpz_poly_divides(Q, A, len1, res, hlen) && _fmpz_poly_divides(Q, B, len2, res, hlen)) 
				   break;
         } else
         {
            /*
				   Bound is easily derived from Thm 5.3 and Cor 5.4 of 
				   http://compalg.inf.elte.hu/~tony/Informatikai-Konyvtar/03-Algorithms%20of%20Informatics%201,%202,%203/CompAlg29May.pdf
					The + 1 is to allow for signed coefficients after Chinese Remaindering
				*/
				bound = hlen + FLINT_MIN(nb1, nb2) + gbits + 1;

            if (pbits >= bound) /* if we reach the bound with one prime */
            { 
               fmpz_t hc;
               fmpz_init(hc);
               _fmpz_vec_content(hc, res, hlen);
               
               /* divide by content */
               _fmpz_vec_scalar_divexact_fmpz(res, res, hlen, hc);
               fmpz_clear(hc); 

               break;
            }

				if (pbits >= bits_small) /* if one prime is already big enough to check */
				{
               fmpz_t hc;
               fmpz_init(hc);
               
               /* divide by content */
               _fmpz_vec_content(hc, res, hlen);
               _fmpz_vec_scalar_divexact_fmpz(res, res, hlen, hc);
               
               /* are we done? */
               if (_fmpz_poly_divides(Q, A, len1, res, hlen) && _fmpz_poly_divides(Q, B, len2, res, hlen)) 
					{
						fmpz_clear(hc); 
						break;
					}

					/* no, so multiply by content again */
               _fmpz_vec_scalar_mul_fmpz(res, res, hlen, hc);
					fmpz_clear(hc); 
				}
         }

         curr_bits = FLINT_ABS(_fmpz_vec_max_bits(res, hlen));
         fmpz_set_ui(modulus, p);
         n = hlen - 1; /* if we reach this we have a new bound on length of result */
         continue;
      }
      
      _fmpz_poly_CRT_ui(res, res, hlen, modulus, h, hlen, mod.n, mod.ninv, 1);
      fmpz_mul_ui(modulus, modulus, mod.n);

      new_bits = FLINT_ABS(_fmpz_vec_max_bits(res, hlen));
      if (new_bits == curr_bits || fmpz_bits(modulus) >= bits_small)
      {
			if (!g_pm1)
         {
            fmpz_init(hc);
            _fmpz_vec_content(hc, res, hlen);
               
            /* divide by content */
            _fmpz_vec_scalar_divexact_fmpz(res, res, hlen, hc);      
         }

         if (fmpz_bits(modulus) > bound)
         {
            if (!g_pm1) fmpz_clear(hc);
            break;
         }
         
         /* are we done? */
         if (_fmpz_poly_divides(Q, A, len1, res, hlen) && _fmpz_poly_divides(Q, B, len2, res, hlen)) 
         {
            if (!g_pm1) fmpz_clear(hc);
            break;
         }

         if (!g_pm1) 
         {            
            /* no, so multiply by content again */
            _fmpz_vec_scalar_mul_fmpz(res, res, hlen, hc);
			   fmpz_clear(hc);
         }
      }
      
      curr_bits = new_bits;
   }
   
   fmpz_clear(modulus);
   fmpz_clear(g); 

   _nmod_vec_clear(a);
   _nmod_vec_clear(b);
   _nmod_vec_clear(h); 
        
   /* finally multiply by content */
   _fmpz_vec_scalar_mul_fmpz(res, res, hlen, d);
   fmpz_clear(d);

   _fmpz_vec_clear(A, len1);
   _fmpz_vec_clear(B, len2);
   _fmpz_vec_clear(Q, len1);
}

void
fmpz_poly_gcd_modular(fmpz_poly_t res,
                           const fmpz_poly_t poly1, const fmpz_poly_t poly2)
{
    const long len1 = poly1->length;
    const long len2 = poly2->length;
    long rlen;
    
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

        return;
    }
    else if (len2 == 0)
    {
        if (fmpz_sgn(poly1->coeffs + (len1 - 1)) > 0)
            fmpz_poly_set(res, poly1);
        else
            fmpz_poly_neg(res, poly1);

        return;
    }

    rlen = FLINT_MIN(len1, len2);

    if (res == poly1 || res == poly2)
    {
       fmpz_poly_t temp;
       fmpz_poly_init2(temp, rlen);
       if (len1 >= len2)
          _fmpz_poly_gcd_modular(temp->coeffs, poly1->coeffs, len1,
                                    poly2->coeffs, len2);
       else
          _fmpz_poly_gcd_modular(temp->coeffs, poly2->coeffs, len2,
                                    poly1->coeffs, len1);
       fmpz_poly_swap(temp, res);
       fmpz_poly_clear(temp);
    }
    else
    {
       fmpz_poly_fit_length(res, rlen);
       if (len1 >= len2)
          _fmpz_poly_gcd_modular(res->coeffs, poly1->coeffs, len1,
                                    poly2->coeffs, len2);
       else
          _fmpz_poly_gcd_modular(res->coeffs, poly2->coeffs, len2,
                                    poly1->coeffs, len1);
    }
    
    _fmpz_poly_set_length(res, rlen);
    _fmpz_poly_normalise(res);
}
