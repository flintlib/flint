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

#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "mpn_extras.h"

long _nmod_poly_xgcd_euclidean(mp_ptr res, mp_ptr s, mp_ptr t, 
              mp_srcptr poly1, long len1, mp_srcptr poly2, long len2, nmod_t mod)
{
   mp_ptr Q, R, u1, u2, v1, v2, prod, A, B;
   mp_limb_t p;
   long u1_len, u2_len, v1_len, v2_len, Q_len, R_len, A_len, B_len, len;
   int steps;
   
   /* clear s, t and res */
   mpn_zero(res, len2);
   mpn_zero(s, len2);
   mpn_zero(t, len1);
   
   p = mod.n;
  
   /* initialise arrays to store intermediate info */
   Q = _nmod_vec_init(len1);
   R = _nmod_vec_init(len2);
   u1 = _nmod_vec_init(len1);
   u2 = _nmod_vec_init(len1);
   v1 = _nmod_vec_init(len1);
   v2 = _nmod_vec_init(len1);
   prod = _nmod_vec_init(len1);
   mpn_zero(u1, len1);
   mpn_zero(u2, len1);
   mpn_zero(v1, len1);
   mpn_zero(v2, len1);
   mpn_zero(prod, len1);
   
   u1[0] = 1;
   u1_len = 1;
   u2_len = 0;
   v2[0] = 1;
   v2_len = 1;
   v1_len = 0;    

   steps = 0;

   A = (mp_ptr) poly1;
   A_len = len1;
   B = (mp_ptr) poly2;
   B_len = len2;
  
   while (B_len > 1)
   {
      _nmod_poly_divrem(Q, R, A, A_len, B, B_len, mod);
      Q_len = A_len - B_len + 1;
      R_len = B_len - 1;
      MPN_NORM(Q, Q_len);
      MPN_NORM(R, R_len);

      if (Q_len && u2_len) /* u1 = u1 - Q*u2 */
      {
          if (Q_len >= u2_len)
              _nmod_poly_mul(prod, Q, Q_len, u2, u2_len, mod);
          else
              _nmod_poly_mul(prod, u2, u2_len, Q, Q_len, mod);
          _nmod_poly_sub(u1, u1, u1_len, prod, u2_len + Q_len - 1, mod);
          u1_len = FLINT_MAX(u1_len, u2_len + Q_len - 1);
          MPN_NORM(u1, u1_len);
      }

      MPN_SWAP(u1, u1_len, u2, u2_len);
      
      if (Q_len && v2_len) /* v1 = v1 - Q*v2 */
      {
          if (Q_len >= v2_len)
              _nmod_poly_mul(prod, Q, Q_len, v2, v2_len, mod);
          else
              _nmod_poly_mul(prod, v2, v2_len, Q, Q_len, mod);
          _nmod_poly_sub(v1, v1, v1_len, prod, v2_len + Q_len - 1, mod);
          v1_len = FLINT_MAX(v1_len, v2_len + Q_len - 1);
          MPN_NORM(v1, v1_len);
      }

      MPN_SWAP(v1, v1_len, v2, v2_len);
      
      MPN_SWAP(A, A_len, B, B_len);
      MPN_SWAP(B, B_len, R, R_len);

      if (steps < 2)
      {
          R = _nmod_vec_init(B_len); /* initialise to original R_len */ 
          steps++;
      }
   }
   
   if (B_len == 1) 
   {      
      MPN_SWAP(u1, u1_len, u2, u2_len);
      MPN_SWAP(v1, v1_len, v2, v2_len);

      mpn_copyi(res, B, B_len);
      len = B_len;
   } else 
   {
      mpn_copyi(res, A, A_len);
      len = A_len;
   }

   mpn_copyi(s, u1, u1_len);
   mpn_copyi(t, v1, v1_len);
      
   if (steps == 2) 
      _nmod_vec_free(A);

   _nmod_vec_free(u1);
   _nmod_vec_free(u2);
   _nmod_vec_free(v1);
   _nmod_vec_free(v2);
   _nmod_vec_free(prod);
   _nmod_vec_free(B);
   _nmod_vec_free(R);
   _nmod_vec_free(Q);

   return len;
}

void
nmod_poly_xgcd_euclidean(nmod_poly_t G, nmod_poly_t S, nmod_poly_t T,
                 const nmod_poly_t A, const nmod_poly_t B)
{
    nmod_poly_t tG, tS, tT;
    mp_ptr g, s, t;
    long A_len, B_len, len;
    mp_limb_t inv;

    B_len = B->length;
    A_len = A->length;
    
    if (A_len == 0)
    {
        if (B_len == 0) 
        {
            nmod_poly_zero(G);
            nmod_poly_zero(S);
            nmod_poly_zero(T);
        }
        else 
        {
            inv = n_invmod(B->coeffs[B_len - 1], B->mod.n);
            nmod_poly_scalar_mul_nmod(G, B, inv);
            nmod_poly_zero(S);
            nmod_poly_set_coeff_ui(T, 0, inv);
            T->length = 1;
        }
        return;
    } 
    else if (B_len == 0)
    {
        inv = n_invmod(A->coeffs[A_len - 1], A->mod.n);
        nmod_poly_scalar_mul_nmod(G, A, inv);
        nmod_poly_zero(T);
        nmod_poly_set_coeff_ui(S, 0, inv);
        S->length = 1;
        return;
    }

    if (A_len == 1)
    {
        nmod_poly_set_coeff_ui(G, 0, 1);
        G->length = 1;
        nmod_poly_zero(T);
        inv = n_invmod(A->coeffs[0], A->mod.n);
        nmod_poly_set_coeff_ui(S, 0, inv);
        S->length = 1;
        return;
    }

    if (B_len == 1)
    {
        nmod_poly_set_coeff_ui(G, 0, 1);
        G->length = 1;
        nmod_poly_zero(S);
        inv = n_invmod(B->coeffs[0], B->mod.n);
        nmod_poly_set_coeff_ui(T, 0, inv);
        T->length = 1;
        return;
    }

    if (G == A || G == B)
    {
        nmod_poly_init2(tG, A->mod.n, FLINT_MIN(A_len, B_len));
        g = tG->coeffs;
    }
    else
    {
        nmod_poly_fit_length(G, FLINT_MIN(A_len, B_len));
        g = G->coeffs;
    }

    if (S == A || S == B)
    {
        nmod_poly_init2(tS, A->mod.n, B_len);
        s = tS->coeffs;
    }
    else
    {
        nmod_poly_fit_length(S, B_len);
        s = S->coeffs;
    }

    if (T == A || T == B)
    {
        nmod_poly_init2(tT, A->mod.n, A_len);
        t = tT->coeffs;
    }
    else
    {
        nmod_poly_fit_length(T, A_len);
        t = T->coeffs;
    }

    if (A_len >= B_len)
        len = _nmod_poly_xgcd_euclidean(g, s, t, A->coeffs, A_len,
                            B->coeffs, B_len, A->mod);
    else
        len = _nmod_poly_xgcd_euclidean(g, t, s, B->coeffs, B_len,
                            A->coeffs, A_len, A->mod);

    if (G == A || G == B)
    {
        nmod_poly_swap(tG, G);
        nmod_poly_clear(tG);
    }
    
    if (S == A || S == B)
    {
        nmod_poly_swap(tS, S);
        nmod_poly_clear(tS);
    }
    
    if (T == A || T == B)
    {
        nmod_poly_swap(tT, T);
        nmod_poly_clear(tT);
    }
    
    MPN_NORM(S->coeffs, B_len);
    S->length = B_len;

    MPN_NORM(T->coeffs, A_len);
    T->length = A_len;

    G->length = len;

    inv = n_invmod(G->coeffs[len - 1], A->mod.n);
    nmod_poly_scalar_mul_nmod(G, G, inv);
    nmod_poly_scalar_mul_nmod(S, S, inv);
    nmod_poly_scalar_mul_nmod(T, T, inv);
}
