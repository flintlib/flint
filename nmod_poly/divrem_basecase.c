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

#include <mpir.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

void _nmod_poly_divrem_basecase_1(mp_ptr Q, mp_ptr R, 
		mp_srcptr A, ulong A_len, mp_srcptr B, ulong B_len, nmod_t mod)
{
   mp_limb_t lead_inv = n_invmod(B[B_len - 1], mod.n);
   mp_limb_t c;

   mp_ptr R_sub;

   ulong coeff = A_len - 1;
   mp_ptr coeff_Q = Q - B_len + 1;

   mpn_copyi(R, A, A_len);
   
   while (coeff + 1 >= B_len)
   {
      R[coeff] = n_mod2_preinv(R[coeff], mod.n, mod.ninv);
  
      while ((coeff + 1 >= B_len) && (R[coeff] == 0L))
      {
         coeff_Q[coeff--] = 0L;
         if (coeff + 1 >= B_len) 
            R[coeff] = n_mod2_preinv(R[coeff], mod.n, mod.ninv);
      }
      
      if (coeff + 1 >= B_len)
      {
         coeff_Q[coeff] = n_mulmod2_preinv(R[coeff], lead_inv, mod.n, mod.ninv); 
         
		 c = n_negmod(coeff_Q[coeff], mod.n);

         R_sub = R + coeff - B_len + 1;
         if (B_len > 1) mpn_addmul_1(R_sub, B, B_len - 1, c);
		 
         coeff--;
      }
   }
   
   while (coeff + 1 > 0)
   {
      R[coeff] = n_mod2_preinv(R[coeff], mod.n, mod.ninv);
      coeff--;
   }
}

void _nmod_poly_divrem_basecase_2(mp_ptr Q, mp_ptr R, 
		mp_srcptr A, ulong A_len, mp_srcptr B, ulong B_len, nmod_t mod)
{
   ulong i;
   mp_limb_t lead_inv = n_invmod(B[B_len - 1], mod.n);
   mp_limb_t c;
   mp_limb_t r_coeff;

   mp_ptr R_sub;
 
   mp_ptr B2 = nmod_vec_init(2*B_len - 2);
   for (i = 0; i < B_len - 1; i++)
   {
      B2[2*i] = B[i];
	  B2[2*i + 1] = 0;
   }

   mp_ptr R2 = nmod_vec_init(2*A_len);
   for (i = 0; i < A_len; i++)
   {
      R2[2*i] = A[i];
	  R2[2*i + 1] = 0;
   }

   ulong coeff = A_len - 1;
   mp_ptr coeff_Q = Q - B_len + 1;
  
   while (coeff + 1 >= B_len)
   {
      r_coeff = n_ll_mod_preinv(R2[2*coeff+1], R2[2*coeff], mod.n, mod.ninv);
  
      while ((coeff + 1 >= B_len) && (r_coeff == 0L))
      {
         coeff_Q[coeff--] = 0L;
         if (coeff + 1 >= B_len) 
            r_coeff = n_ll_mod_preinv(R2[2*coeff+1], R2[2*coeff], mod.n, mod.ninv);
      }
      
      if (coeff + 1 >= B_len)
      {
         coeff_Q[coeff] = n_mulmod2_preinv(r_coeff, lead_inv, mod.n, mod.ninv); 
         
		 c = n_negmod(coeff_Q[coeff], mod.n);

         R_sub = R2 + 2*(coeff - B_len + 1);
         if (B_len > 1) mpn_addmul_1(R_sub, B2, 2*B_len - 2, c);
		 
         coeff--;
      }
   }

   while (coeff + 1 > 0)
   {
      R[coeff] = n_ll_mod_preinv(R2[2*coeff+1], R2[2*coeff], mod.n, mod.ninv);
	  coeff--;
   }

   nmod_vec_free(B2);
   nmod_vec_free(R2);
}

void _nmod_poly_divrem_basecase_3(mp_ptr Q, mp_ptr R, 
		mp_srcptr A, ulong A_len, mp_srcptr B, ulong B_len, nmod_t mod)
{
   ulong i;
   mp_limb_t lead_inv = n_invmod(B[B_len - 1], mod.n);
   mp_limb_t c;
   mp_limb_t r_coeff;

   mp_ptr R_sub;
 
   mp_ptr B3 = nmod_vec_init(3*B_len - 3);
   for (i = 0; i < B_len - 1; i++)
   {
      B3[3*i] = B[i];
	  B3[3*i + 1] = 0;
	  B3[3*i + 2] = 0;
   }

   mp_ptr R3 = nmod_vec_init(3*A_len);
   for (i = 0; i < A_len; i++)
   {
      R3[3*i] = A[i];
	  R3[3*i + 1] = 0;
	  R3[3*i + 2] = 0;
   }

   ulong coeff = A_len - 1;
   mp_ptr coeff_Q = Q - B_len + 1;
  
   while (coeff + 1 >= B_len)
   {
	  r_coeff = n_lll_mod_preinv(R3[3*coeff+2], R3[3*coeff+1], R3[3*coeff], mod.n, mod.ninv);
  
      while ((coeff + 1 >= B_len) && (r_coeff == 0L))
      {
         coeff_Q[coeff--] = 0L;
         if (coeff + 1 >= B_len) 
            r_coeff = n_lll_mod_preinv(R3[3*coeff+2], R3[3*coeff+1], R3[3*coeff], mod.n, mod.ninv);
      }
      
      if (coeff + 1 >= B_len)
      {
         coeff_Q[coeff] = n_mulmod2_preinv(r_coeff, lead_inv, mod.n, mod.ninv); 
         
		 c = n_negmod(coeff_Q[coeff], mod.n);

         R_sub = R3 + 3*(coeff - B_len + 1);
         if (B_len > 1) mpn_addmul_1(R_sub, B3, 3*B_len - 3, c);
		 
         coeff--;
      }
   }

   while (coeff + 1 > 0)
   {
      R[coeff] = n_lll_mod_preinv(R3[3*coeff+2], R3[3*coeff+1], R3[3*coeff], mod.n, mod.ninv);
      coeff--;
   }

   nmod_vec_free(B3);
   nmod_vec_free(R3);
}

void _nmod_poly_divrem_basecase(mp_ptr Q, mp_ptr R, 
		mp_srcptr A, ulong A_len, mp_srcptr B, ulong B_len, nmod_t mod)
{
   ulong bits = 2*(FLINT_BITS - mod.norm) + FLINT_BIT_COUNT(A_len - B_len + 1);

   if (bits <= FLINT_BITS)
      _nmod_poly_divrem_basecase_1(Q, R, A, A_len, B, B_len, mod);
   else if (bits <= 2*FLINT_BITS)
      _nmod_poly_divrem_basecase_2(Q, R, A, A_len, B, B_len, mod);
   else
      _nmod_poly_divrem_basecase_3(Q, R, A, A_len, B, B_len, mod);
}

void nmod_poly_divrem_basecase(nmod_poly_t Q, nmod_poly_t R, nmod_poly_t A, nmod_poly_t B)
{
   mp_ptr Q_coeffs, R_coeffs;
   nmod_poly_t t1, t2;

   if (B->length == 0)
   {
      printf("Exception: Divide by zero in nmod_poly_divrem_basecase\n");
      abort();      
   }
   
   if (A->length < B->length)
   {
      nmod_poly_set(R, A);
      nmod_poly_zero(Q);
      
      return;
   }

   if ((Q == A) || (Q == B))
   {
	   nmod_poly_init2_preinv(t1, B->mod.n, B->mod.ninv, A->length - B->length + 1);
	   Q_coeffs = t1->coeffs;
   } else
   {
	  nmod_poly_fit_length(Q, A->length - B->length + 1);
      Q_coeffs = Q->coeffs;	  
   }
   
   if ((R == A) || (R == B))
   {
	  nmod_poly_init2_preinv(t2, B->mod.n, B->mod.ninv, A->length);
	  R_coeffs = t2->coeffs;
   } else
   {
	  nmod_poly_fit_length(R, A->length);
      R_coeffs = R->coeffs;	  
   }
   
   _nmod_poly_divrem_basecase(Q_coeffs, R_coeffs, A->coeffs, A->length, 
	                                                      B->coeffs, B->length, B->mod);
   
   if ((Q == A) || (Q == B))
   {
      t1->length = A->length - B->length + 1;
      nmod_poly_swap(Q, t1);
	  nmod_poly_clear(t1);
   } else
      Q->length = A->length - B->length + 1;

   if ((R == A) || (R == B))
   {
      t2->length = B->length - 1;
      nmod_poly_swap(R, t2);
	  nmod_poly_clear(t2);
   } else
      R->length = B->length - 1;
   
   _nmod_poly_normalise(Q);
   _nmod_poly_normalise(R);
}
