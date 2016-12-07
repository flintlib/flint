/*
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"
#include "hashmap.h"

void _fmpz_mpoly_addmul_array1_slong1(ulong * poly1, 
                 const slong * poly2, const ulong * exp2, slong len2,
                           const slong * poly3, const ulong * exp3, slong len3)
{
   slong i, j;
   ulong * c2;

   for (i = 0; i < len2; i++)
   {
      c2 = poly1 + (slong) exp2[i];

      for (j = 0; j < len3; j++)
         c2[(slong) exp3[j]] += poly2[i]*poly3[j];
   }
}

void _fmpz_mpoly_addmul_array1_slong(ulong * poly1, 
                 const slong * poly2, const ulong * exp2, slong len2,
                           const slong * poly3, const ulong * exp3, slong len3)
{
   slong i, j;
   ulong cy;
   ulong p[2]; /* for products of coefficients */
   ulong * c2, * c;

   for (i = 0; i < len2; i++)
   {
      c2 = poly1 + 3*((slong) exp2[i]);

      for (j = 0; j < len3; j++)
      {
         c = c2 + 3*((slong) exp3[j]);

         smul_ppmm(p[1], p[0], poly2[i], poly3[j]);
         add_sssaaaaaa(cy, c[1], c[0], 0, c[1], c[0], 0, p[1], p[0]);
         c[2] += (0 <= (slong) p[1]) ? cy : cy - 1;
      }
   }
}

void _fmpz_mpoly_addmul_array1_slong2(ulong * poly1, 
                 const slong * poly2, const ulong * exp2, slong len2,
                           const slong * poly3, const ulong * exp3, slong len3)
{
   slong i, j;
   ulong p[2]; /* for products of coefficients */
   ulong * c2, * c;

   for (i = 0; i < len2; i++)
   {
      c2 = poly1 + 2*((slong) exp2[i]);

      for (j = 0; j < len3; j++)
      {
         c = c2 + 2*((slong) exp3[j]);

         smul_ppmm(p[1], p[0], poly2[i], poly3[j]);
         add_ssaaaa(c[1], c[0], c[1], c[0], p[1], p[0]);
      }
   }
}

void _fmpz_mpoly_addmul_array1_fmpz(fmpz * poly1, 
                 const fmpz * poly2, const ulong * exp2, slong len2,
                           const fmpz * poly3, const ulong * exp3, slong len3)
{
   slong i, j;
   fmpz * c2, * c;

   for (i = 0; i < len2; i++)
   {
      c2 = poly1 + (slong) exp2[i];

      for (j = 0; j < len3; j++)
      {
         c = c2 + (slong) exp3[j];

         fmpz_addmul(c, poly2 + i, poly3 + i);
      }
   }
}

/* this function destroys poly2 and starts writing at index k */
slong _fmpz_mpoly_from_ulong_array(fmpz ** poly1, ulong ** exp1, slong * alloc, 
              ulong * poly2, const slong * mults, slong N, slong bits, slong k)
{
   slong i, j;
   ulong exp;
   ulong * c;
   slong * prods;
   int negate;
   fmpz * p1 = *poly1;
   ulong * e1 = *exp1;
   slong shift = FLINT_BITS - N*bits;
   TMP_INIT;

   TMP_START;

   prods = (slong *) TMP_ALLOC((N + 1)*sizeof(slong));

   prods[0] = 1;
   for (i = 1; i <= N; i++)
     prods[i] = mults[i - 1]*prods[i - 1];
   
   for (i = 0; i < prods[N]; i++)
   {
      c = poly2 + i*3;

      if (c[0] != 0 || c[1] != 0 || c[2] != 0)
      {
         if (k >= *alloc)
         {
            p1 = (fmpz *) flint_realloc(p1, 2*(*alloc)*sizeof(fmpz));
            e1 = (ulong *) flint_realloc(e1, 2*(*alloc)*sizeof(ulong));
            flint_mpn_zero(p1 + *alloc, *alloc);
            (*alloc) *= 2;
         }

         exp = 0;
         
         for (j = 0; j < N; j++)
            exp += (i % prods[j + 1])/prods[j] << bits*j;

         e1[k] = exp << shift;
         
         negate = 0;

         if (0 > (slong) c[2])
         {
            c[0] = ~c[0];
            c[1] = ~c[1];
            c[2] = ~c[2];
            add_sssaaaaaa(c[2], c[1], c[0], c[2], c[1], c[0], 0, 0, 1);
            negate = 1;
         } 

         fmpz_set_ui(p1 + k, c[2]);
         fmpz_mul_2exp(p1 + k, p1 + k, FLINT_BITS);
         fmpz_add_ui(p1 + k, p1 + k, c[1]);
         fmpz_mul_2exp(p1 + k, p1 + k, FLINT_BITS);
         fmpz_add_ui(p1 + k, p1 + k, c[0]);
      
         if (negate)
            fmpz_neg(p1 + k, p1 + k);
         
         k++;
      }
   }

   *poly1 = p1;
   *exp1 = e1;

   TMP_END;

   return k;
}

/* this function destroys poly2 and starts writing at index k */
slong _fmpz_mpoly_from_ulong_array2(fmpz ** poly1, ulong ** exp1, slong * alloc, 
              ulong * poly2, const slong * mults, slong N, slong bits, slong k)
{
   slong i, j;
   ulong exp;
   ulong * c;
   slong * prods;
   int negate;
   fmpz * p1 = *poly1;
   ulong * e1 = *exp1;
   slong shift = FLINT_BITS - N*bits;
   TMP_INIT;

   TMP_START;

   prods = (slong *) TMP_ALLOC((N + 1)*sizeof(slong));

   prods[0] = 1;
   for (i = 1; i <= N; i++)
     prods[i] = mults[i - 1]*prods[i - 1];
   
   for (i = 0; i < prods[N]; i++)
   {
      c = poly2 + i*2;

      if (c[0] != 0 || c[1] != 0)
      {
         if (k >= *alloc)
         {
            p1 = (fmpz *) flint_realloc(p1, 2*(*alloc)*sizeof(fmpz));
            e1 = (ulong *) flint_realloc(e1, 2*(*alloc)*sizeof(ulong));
            flint_mpn_zero(p1 + *alloc, *alloc);
            (*alloc) *= 2;
         }

         exp = 0;
         
         for (j = 0; j < N; j++)
            exp += (i % prods[j + 1])/prods[j] << bits*j;

         e1[k] = exp << shift;
         
         negate = 0;

         if (0 > (slong) c[1])
         {
            c[0] = ~c[0];
            c[1] = ~c[1];
            add_ssaaaa(c[1], c[0], c[1], c[0], 0, 1);
            negate = 1;
         } 

         fmpz_set_ui(p1 + k, c[1]);
         fmpz_mul_2exp(p1 + k, p1 + k, FLINT_BITS);
         fmpz_add_ui(p1 + k, p1 + k, c[0]);
         
         if (negate)
            fmpz_neg(p1 + k, p1 + k);
         
         k++;
      }
   }

   *poly1 = p1;
   *exp1 = e1;

   TMP_END;

   return k;
}

/* this function destroys poly2 and starts writing at index k */
slong _fmpz_mpoly_from_ulong_array1(fmpz ** poly1, ulong ** exp1, slong * alloc, 
              ulong * poly2, const slong * mults, slong N, slong bits, slong k)
{
   slong i, j;
   ulong exp;
   ulong * c;
   slong * prods;
   fmpz * p1 = *poly1;
   ulong * e1 = *exp1;
   slong shift = FLINT_BITS - N*bits;
   TMP_INIT;

   TMP_START;

   prods = (slong *) TMP_ALLOC((N + 1)*sizeof(slong));

   prods[0] = 1;
   for (i = 1; i <= N; i++)
     prods[i] = mults[i - 1]*prods[i - 1];
   
   for (i = 0; i < prods[N]; i++)
   {
      c = poly2 + i;

      if (c[0] != 0)
      {
         if (k >= *alloc)
         {
            p1 = (fmpz *) flint_realloc(p1, 2*(*alloc)*sizeof(fmpz));
            e1 = (ulong *) flint_realloc(e1, 2*(*alloc)*sizeof(ulong));
            flint_mpn_zero(p1 + *alloc, *alloc);
            (*alloc) *= 2;
         }

         exp = 0;
         
         for (j = 0; j < N; j++)
            exp += (i % prods[j + 1])/prods[j] << bits*j;

         e1[k] = exp << shift;
         
         fmpz_set_si(p1 + k, c[0]);
         
         k++;
      }
   }

   *poly1 = p1;
   *exp1 = e1;

   TMP_END;

   return k;
}

slong _fmpz_mpoly_from_fmpz_array(fmpz ** poly1, ulong ** exp1, slong * alloc, 
               fmpz * poly2, const slong * mults, slong N, slong bits, slong k)
{
   slong i, j;
   ulong exp;
   fmpz * c;
   slong * prods;
   fmpz * p1 = *poly1;
   ulong * e1 = *exp1;
   slong shift = FLINT_BITS - N*bits;
   TMP_INIT;

   TMP_START;

   prods = (slong *) TMP_ALLOC((N + 1)*sizeof(slong));

   prods[0] = 1;
   for (i = 1; i <= N; i++)
      prods[i] = mults[i - 1]*prods[i - 1];

   for (i = 0; i < prods[N]; i++)
   {
      c = poly2 + i;

      if (!fmpz_is_zero(c))
      {
         if (k >= *alloc)
         {
            p1 = (fmpz *) flint_realloc(p1, 2*(*alloc)*sizeof(fmpz));
            e1 = (ulong *) flint_realloc(e1, 2*(*alloc)*sizeof(ulong));
            flint_mpn_zero(p1 + *alloc, *alloc);
            (*alloc) *= 2;
         }

         exp = 0;
         
         for (j = 0; j < N; j++)
            exp += (i % prods[j + 1])/prods[j] << bits*j;

         e1[k] = exp << shift;
         
         fmpz_set(p1 + k, poly2 + i);
         
         k++;
      }
   }

   *poly1 = p1;
   *exp1 = e1;

   TMP_END;

   return k;
}

void _fmpz_mpoly_pack_exponents_tight(ulong * exp1, const ulong * exp2,
                           slong len, const slong * mults, slong N, slong bits)
{
   slong i, j;
   ulong e1, e2;
   ulong mask = (UWORD(1) << bits) - 1;
   slong shift = FLINT_BITS - N*bits;

   for (i = 0; i < len; i++)
   {
      e2 = exp2[i] >> shift;
      e1 = ((e2 >> (N - 1)*bits) & mask);
      
      for (j = N - 2; j >= 0; j--)
      {
         e1 *= mults[j];
         e1 += ((e2 >> j*bits) & mask);
      }

      exp1[i] = e1;
   }
}

/*
   use array multiplication to set poly1 to poly2*poly3 in N variables, given
   a list of multipliers to tightly pack exponents and a number of bits for the
   fields of the exponents of the result, assuming no aliasing
*/
slong _fmpz_mpoly_mul_array(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                  const fmpz_mpoly_t poly3, slong * mults, slong bits, slong N)
{
   slong i, bits1, bits2, bits3;
   ulong * exp2, * exp3;
   slong prod, len;
   int small, sign;
   TMP_INIT;

   TMP_START;

   exp2 = (ulong *) TMP_ALLOC(poly2->length*sizeof(ulong));
   exp3 = (ulong *) TMP_ALLOC(poly3->length*sizeof(ulong));

   _fmpz_mpoly_pack_exponents_tight(exp2, poly2->exps, poly2->length,
                                                        mults, N, poly2->bits);

   _fmpz_mpoly_pack_exponents_tight(exp3, poly3->exps, poly3->length,
                                                        mults, N, poly3->bits);

   bits2 = fmpz_mpoly_max_bits(poly2);
   bits3 = fmpz_mpoly_max_bits(poly2);

   sign = (bits2 < 0) || (bits3 < 0);

   bits1 = FLINT_ABS(bits2) + FLINT_ABS(bits3) +
          FLINT_BIT_COUNT(FLINT_MIN(poly2->length, poly3->length)) + sign;

   small = FLINT_ABS(bits2) <= 62 && FLINT_ABS(bits3) <= 62;

   prod = 1;
   for (i = 0; i < N; i++)
      prod *= mults[i];

   if (small && bits1 <= FLINT_BITS)
   {
      ulong * p1 = (ulong *) TMP_ALLOC(prod*sizeof(ulong));

      for (i = 0; i < prod; i++)
         p1[i] = 0;

      _fmpz_mpoly_addmul_array1_slong1(p1, 
                   (slong *) poly2->coeffs, exp2, poly2->length,
                                 (slong *) poly3->coeffs, exp3, poly3->length);

      len = _fmpz_mpoly_from_ulong_array1(&poly1->coeffs, &poly1->exps, &poly1->alloc, 
                                                        p1, mults, N, bits, 0);
   } else if (small && bits1 <= 2*FLINT_BITS)
   {
      ulong * p1 = (ulong *) TMP_ALLOC(2*prod*sizeof(ulong));

      for (i = 0; i < 2*prod; i++)
         p1[i] = 0;

      _fmpz_mpoly_addmul_array1_slong2(p1, 
                   (slong *) poly2->coeffs, exp2, poly2->length,
                                 (slong *) poly3->coeffs, exp3, poly3->length);

      len = _fmpz_mpoly_from_ulong_array2(&poly1->coeffs, &poly1->exps, &poly1->alloc, 
                                                        p1, mults, N, bits, 0);
   } else if (small)
   {
      ulong * p1 = (ulong *) TMP_ALLOC(3*prod*sizeof(ulong));

      for (i = 0; i < 3*prod; i++)
         p1[i] = 0;

      _fmpz_mpoly_addmul_array1_slong(p1, 
                   (slong *) poly2->coeffs, exp2, poly2->length,
                                 (slong *) poly3->coeffs, exp3, poly3->length);

      len = _fmpz_mpoly_from_ulong_array(&poly1->coeffs, &poly1->exps, &poly1->alloc, 
                                                        p1, mults, N, bits, 0);      
   } else
   {
      fmpz * p1 = (fmpz *) TMP_ALLOC(prod*sizeof(fmpz));

      for (i = 0; i < prod; i++)
         p1[i] = 0;

      _fmpz_mpoly_addmul_array1_fmpz(p1, 
                            poly2->coeffs, exp2, poly2->length,
                                           poly3->coeffs, exp3, poly3->length);

      len = _fmpz_mpoly_from_fmpz_array(&poly1->coeffs, &poly1->exps, &poly1->alloc, 
                                                        p1, mults, N, bits, 0);
   }
  
   TMP_END;

   return len;
}

int fmpz_mpoly_mul_array(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                          const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx)
{
   slong i, bits, exp_bits, N, len = 0, array_size;
   ulong * max_degs2;
   ulong * max_degs3;
   ulong max2 = 0, max3 = 0, max;
   int res = 1;

   TMP_INIT;

   if (poly2->length == 0 || poly3->length == 0)
   {
      fmpz_mpoly_zero(poly1, ctx);

      return 1;
   }

   TMP_START;

   max_degs2 = (ulong *) TMP_ALLOC(ctx->n*sizeof(ulong));
   max_degs3 = (ulong *) TMP_ALLOC(ctx->n*sizeof(ulong));

   fmpz_mpoly_max_degrees(max_degs2, poly2, ctx);
   fmpz_mpoly_max_degrees(max_degs3, poly3, ctx);

   for (i = 0; i < ctx->n; i++)
   {
      if (max_degs2[i] > max2)
         max2 = max_degs2[i];

      if (max_degs3[i] > max3)
         max3 = max_degs3[i];
   }

   max = max2 + max3;
   if (max < max2 || 0 > (slong) max)
      flint_throw(FLINT_EXPOF, "Exponent overflow in fmpz_mpoly_mul_array");

   bits = FLINT_BIT_COUNT(max);

   exp_bits = 8;
   while (bits >= exp_bits)
      exp_bits *= 2;

   N = (exp_bits*ctx->n - 1)/FLINT_BITS + 1;

   if (N != 1)
   {
      res = 0;
      goto cleanup;
   }

   array_size = 1;
   for (i = 0; i < ctx->n; i++)
   {
      max_degs2[i] += max_degs3[i] + 1;
      array_size *= max_degs2[i];  
   }

   if (array_size > WORD(100000000))
   {
      res = 0;
      goto cleanup;
   }

   if (poly1 == poly2 || poly1 == poly3)
   {
      fmpz_mpoly_t temp;

      fmpz_mpoly_init2(temp, poly2->length + poly3->length - 1, ctx);
      fmpz_mpoly_fit_bits(temp, exp_bits, ctx);

      if (poly2->length >= poly3->length)
         len = _fmpz_mpoly_mul_array(temp, poly3, poly2,
                                        (slong *) max_degs2, exp_bits, ctx->n);
      else
         len = _fmpz_mpoly_mul_array(temp, poly2, poly3,
                                        (slong *) max_degs2, exp_bits, ctx->n);

      fmpz_mpoly_swap(temp, poly1, ctx);

      fmpz_mpoly_clear(temp, ctx);
   } else
   {
      fmpz_mpoly_fit_length(poly1, poly2->length + poly3->length - 1, ctx);
      fmpz_mpoly_fit_bits(poly1, exp_bits, ctx);

      if (poly2->length >= poly3->length)
         len = _fmpz_mpoly_mul_array(poly1, poly3, poly2,
                                        (slong *) max_degs2, exp_bits, ctx->n);
      else
         len = _fmpz_mpoly_mul_array(poly1, poly2, poly3,
                                        (slong *) max_degs2, exp_bits, ctx->n);
   }

   _fmpz_mpoly_set_length(poly1, len, ctx);

cleanup:

   TMP_END;

   return res;
}
