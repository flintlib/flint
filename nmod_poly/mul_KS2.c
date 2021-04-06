/*
    Copyright (C) 2007, 2008 David Harvey (zn_poly)
    Copyright (C) 2013 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

/*
   Multiplication/squaring using Kronecker substitution at 2^b and -2^b.
*/
void
_nmod_poly_mul_KS2(mp_ptr res, mp_srcptr op1, slong n1,
                  mp_srcptr op2, slong n2, nmod_t mod)
{
   int sqr, v3m_neg;
   ulong bits, b, w;
   slong n1o, n1e, n2o, n2e, n3o, n3e, n3, k1, k2, k3;
   mp_ptr v1_buf0, v2_buf0, v1_buf1, v2_buf1, v1_buf2, v2_buf2;
   mp_ptr v1o, v1e, v1p, v1m, v2o, v2e, v2p, v2m, v3o, v3e, v3p, v3m;
   mp_ptr z;
   TMP_INIT;

   if (n2 == 1)
   {
      /* code below needs n2 > 1, so fall back on scalar multiplication */
      _nmod_vec_scalar_mul_nmod(res, op1, n1, op2[0], mod);
      return;
   }

   TMP_START;

   sqr = (op1 == op2 && n1 == n2);

   /* bits in each output coefficient */
   bits = 2 * (FLINT_BITS - mod.norm) + FLINT_CLOG2(n2);
   
   /* we're evaluating at x = B and -B, where B = 2^b, and b = ceil(bits / 2) */
   b = (bits + 1) / 2;

   /* number of ulongs required to store each output coefficient */
   w = (2*b - 1)/FLINT_BITS + 1;
   
   /* 
      Write f1(x) = f1e(x^2) + x * f1o(x^2)
            f2(x) = f2e(x^2) + x * f2o(x^2)
             h(x) =  he(x^2) + x *  ho(x^2)
      "e" = even, "o" = odd
   */

   n1o = n1 / 2;
   n1e = n1 - n1o;

   n2o = n2 / 2;
   n2e = n2 - n2o;

   n3 = n1 + n2 - 1;    /* length of h */
   n3o = n3 / 2;
   n3e = n3 - n3o;

   /*
      f1(B) and |f1(-B)| are at most ((n1 - 1) * b + mod->bits) bits long.
      However, when evaluating f1e(B^2) and B * f1o(B^2) the bitpacking
      routine needs room for the last chunk of 2b bits. Therefore we need to
      allow room for (n1 + 1) * b bits. Ditto for f2.
   */
   k1 = ((n1 + 1)*b - 1)/FLINT_BITS + 1;
   k2 = ((n2 + 1)*b - 1)/FLINT_BITS + 1;
   k3 = k1 + k2;

   /* allocate space */
   v1_buf0 = TMP_ALLOC(sizeof(mp_limb_t) * 3 * k3); /* k1 limbs */
   v2_buf0 = v1_buf0 + k1;         /* k2 limbs */
   v1_buf1 = v2_buf0 + k2;         /* k1 limbs */
   v2_buf1 = v1_buf1 + k1;         /* k2 limbs */
   v1_buf2 = v2_buf1 + k2;         /* k1 limbs */
   v2_buf2 = v1_buf2 + k1;         /* k2 limbs */

   /*
      arrange overlapping buffers to minimise memory use
      "p" = plus, "m" = minus
   */
   v1e = v1_buf0;
   v2e = v2_buf0;
   v1o = v1_buf1;
   v2o = v2_buf1;
   v1p = v1_buf2;
   v2p = v2_buf2;
   v1m = v1_buf0;
   v2m = v2_buf0;
   v3m = v1_buf1;
   v3p = v1_buf0;
   v3e = v1_buf2;
   v3o = v1_buf0;
   
   z = TMP_ALLOC(sizeof(mp_limb_t) * w * n3e);
   
   if (!sqr)
   {
      /* multiplication version */

      /* evaluate f1e(B^2) and B * f1o(B^2) */
      _nmod_poly_KS2_pack(v1e, op1, n1e, 2, 2 * b, 0, k1);
      _nmod_poly_KS2_pack(v1o, op1 + 1, n1o, 2, 2 * b, b, k1);

      /* evaluate f2e(B^2) and B * f2o(B^2) */
      _nmod_poly_KS2_pack(v2e, op2, n2e, 2, 2 * b, 0, k2);
      _nmod_poly_KS2_pack(v2o, op2 + 1, n2o, 2, 2 * b, b, k2);

      /*
         compute f1(B) = f1e(B^2) + B * f1o(B^2)
             and f2(B) = f2e(B^2) + B * f2o(B^2)
      */
      mpn_add_n(v1p, v1e, v1o, k1);
      mpn_add_n(v2p, v2e, v2o, k2);

      /*
         compute |f1(-B)| = |f1e(B^2) - B * f1o(B^2)|
             and |f2(-B)| = |f2e(B^2) - B * f2o(B^2)|
      */
      v3m_neg  = signed_mpn_sub_n(v1m, v1e, v1o, k1);
      v3m_neg ^= signed_mpn_sub_n(v2m, v2e, v2o, k2);

      /*
         compute  h(B)   =  f1(B)   *  f2(B)
         compute |h(-B)| = |f1(-B)| * |f2(-B)|
         v3m_neg is set if h(-B) is negative
      */
      mpn_mul(v3m, v1m, k1, v2m, k2);
      mpn_mul(v3p, v1p, k1, v2p, k2);
   }
   else
   {
      /* squaring version */

      /* evaluate f1e(B^2) and B * f1o(B^2) */
      _nmod_poly_KS2_pack(v1e, op1, n1e, 2, 2 * b, 0, k1);
      _nmod_poly_KS2_pack(v1o, op1 + 1, n1o, 2, 2 * b, b, k1);

      /* compute f1(B) = f1e(B^2) + B * f1o(B^2) */
      mpn_add_n(v1p, v1e, v1o, k1);

      /* compute |f1(-B)| = |f1e(B^2) - B * f1o(B^2)| */
      signed_mpn_sub_n(v1m, v1e, v1o, k1);

      /*
         compute h(B)  = f1(B)^2
         compute h(-B) = f1(-B)^2
         v3m_neg is cleared (since f1(-B)^2 is never negative)
      */
      mpn_sqr(v3m, v1m, k1);
      mpn_sqr(v3p, v1p, k1);
      v3m_neg = 0;
   }
   
   /* 
      he(B^2) and B * ho(B^2) are both at most b * (n3 + 1) bits long (since
      the coefficients don't overlap). The buffers used below are at least
      b * (n1 + n2 + 2) = b * (n3 + 3) bits long. So we definitely have
      enough room for 2 * he(B^2) and 2 * B * ho(B^2).
   */

   /* compute 2 * he(B^2) = h(B) + h(-B) */
   if (v3m_neg) 
      mpn_sub_n(v3e, v3p, v3m, k3);
   else
      mpn_add_n(v3e, v3p, v3m, k3);

   /* unpack coefficients of he, and reduce mod m */
   _nmod_poly_KS2_unpack(z, v3e, n3e, 2 * b, 1);
   _nmod_poly_KS2_reduce(res, 2, z, n3e, w, mod);
   
   /* compute 2 * b * ho(B^2) = h(B) - h(-B) */
   if (v3m_neg)
      mpn_add_n(v3o, v3p, v3m, k3);
   else 
      mpn_sub_n(v3o, v3p, v3m, k3);
   
   /* unpack coefficients of ho, and reduce mod m */
   _nmod_poly_KS2_unpack(z, v3o, n3o, 2 * b, b + 1);
   _nmod_poly_KS2_reduce(res + 1, 2, z, n3o, w, mod);

   TMP_END;
}

void
nmod_poly_mul_KS2(nmod_poly_t res,
                 const nmod_poly_t poly1, const nmod_poly_t poly2)
{
    slong len_out;

    if ((poly1->length == 0) || (poly2->length == 0))
    {
        nmod_poly_zero(res);
        return;
    }

    len_out = poly1->length + poly2->length - 1;

    if (res == poly1 || res == poly2)
    {
        nmod_poly_t temp;
        nmod_poly_init2_preinv(temp, poly1->mod.n, poly1->mod.ninv, len_out);
        if (poly1->length >= poly2->length)
            _nmod_poly_mul_KS2(temp->coeffs, poly1->coeffs, poly1->length,
                              poly2->coeffs, poly2->length,
                              poly1->mod);
        else
            _nmod_poly_mul_KS2(temp->coeffs, poly2->coeffs, poly2->length,
                              poly1->coeffs, poly1->length,
                              poly1->mod);
        nmod_poly_swap(res, temp);
        nmod_poly_clear(temp);
    }
    else
    {
        nmod_poly_fit_length(res, len_out);
        if (poly1->length >= poly2->length)
            _nmod_poly_mul_KS2(res->coeffs, poly1->coeffs, poly1->length,
                              poly2->coeffs, poly2->length,
                              poly1->mod);
        else
            _nmod_poly_mul_KS2(res->coeffs, poly2->coeffs, poly2->length,
                              poly1->coeffs, poly1->length,
                              poly1->mod);
    }

    res->length = len_out;
    _nmod_poly_normalise(res);
}
