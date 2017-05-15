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

/* improve locality */
#define BLOCK 128
#define MAX_ARRAY_SIZE (WORD(300000))

void _fmpz_mpoly_addmul_array1_slong1(ulong * poly1, 
                 const slong * poly2, const ulong * exp2, slong len2,
                           const slong * poly3, const ulong * exp3, slong len3)
{
   slong ii, i, jj, j;
   ulong * c2;

   for (ii = 0; ii < len2 + BLOCK; ii += BLOCK)
   {
      for (jj = 0; jj < len3 + BLOCK; jj += BLOCK)
      {
         for (i = ii; i < FLINT_MIN(ii + BLOCK, len2); i++)
         {
            c2 = poly1 + (slong) exp2[i];

            if (poly2[i] != 0)
            {
               for (j = jj; j < FLINT_MIN(jj + BLOCK, len3); j++)
               {
                  c2[(slong) exp3[j]] += poly2[i]*poly3[j];
               }
            }
         }
      }
   }
}

void _fmpz_mpoly_addmul_array1_slong(ulong * poly1, 
                 const slong * poly2, const ulong * exp2, slong len2,
                           const slong * poly3, const ulong * exp3, slong len3)
{
   slong ii, i, jj, j;
   ulong cy;
   ulong p[2]; /* for products of coefficients */
   ulong * c2, * c;

   for (ii = 0; ii < len2 + BLOCK; ii += BLOCK)
   {
      for (jj = 0; jj < len3 + BLOCK; jj += BLOCK)
      {
         for (i = ii; i < FLINT_MIN(ii + BLOCK, len2); i++)
         {
            c2 = poly1 + 3*((slong) exp2[i]);

            if (poly2[i] != 0)
            {
               for (j = jj; j < FLINT_MIN(jj + BLOCK, len3); j++)
               {
                  c = c2 + 3*((slong) exp3[j]);

                  smul_ppmm(p[1], p[0], poly2[i], poly3[j]);
                  add_sssaaaaaa(cy, c[1], c[0], 0, c[1], c[0], 0, p[1], p[0]);
                  c[2] += (0 <= (slong) p[1]) ? cy : cy - 1;
               }
            }
         }
      }
   }
}

void _fmpz_mpoly_addmul_array1_slong2(ulong * poly1, 
                 const slong * poly2, const ulong * exp2, slong len2,
                           const slong * poly3, const ulong * exp3, slong len3)
{
   slong ii, i, jj, j;
   ulong p[2]; /* for products of coefficients */
   ulong * c2, * c;

   for (ii = 0; ii < len2 + BLOCK; ii += BLOCK)
   {
      for (jj = 0; jj < len3 + BLOCK; jj += BLOCK)
      {
         for (i = ii; i < FLINT_MIN(ii + BLOCK, len2); i++)
         {
            c2 = poly1 + 2*((slong) exp2[i]);

            if (poly2[i] != 0)
            {
               for (j = jj; j < FLINT_MIN(jj + BLOCK, len3); j++)
               {
                  c = c2 + 2*((slong) exp3[j]);

                  smul_ppmm(p[1], p[0], poly2[i], poly3[j]);
                  add_ssaaaa(c[1], c[0], c[1], c[0], p[1], p[0]);
               }
            }
         }
      }
   }
}

void _fmpz_mpoly_addmul_array1_fmpz(fmpz * poly1, 
                 const fmpz * poly2, const ulong * exp2, slong len2,
                           const fmpz * poly3, const ulong * exp3, slong len3)
{
   slong ii, i, jj, j;
   fmpz * c2, * c;

   for (ii = 0; ii < len2 + BLOCK; ii += BLOCK)
   {
      for (jj = 0; jj < len3 + BLOCK; jj += BLOCK)
      {
         for (i = ii; i < FLINT_MIN(ii + BLOCK, len2); i++)
         {
            c2 = poly1 + (slong) exp2[i];

            if (poly2[i] != 0)
            {
               for (j = jj; j < FLINT_MIN(jj + BLOCK, len3); j++)
               {
                  c = c2 + (slong) exp3[j];
                  fmpz_addmul(c, poly2 + i, poly3 + i);
               }
            }
         }
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
              ulong * poly2, const slong * mults, slong num, slong bits, slong k)
{
   slong i, j;
   ulong exp;
   ulong * c;
   slong * prods;
   int negate;
   fmpz * p1 = *poly1;
   ulong * e1 = *exp1;
   slong shift = FLINT_BITS - num*bits;
   TMP_INIT;

   TMP_START;

   prods = (slong *) TMP_ALLOC((num + 1)*sizeof(slong));

   prods[0] = 1;
   for (i = 1; i <= num; i++)
     prods[i] = mults[i - 1]*prods[i - 1];
   
   for (i = 0; i < prods[num]; i++)
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
         
         for (j = 0; j < num; j++)
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
              ulong * poly2, const slong * mults, slong num, slong bits, slong k)
{
   slong i, j;
   ulong exp;
   ulong * c;
   slong * prods;
   fmpz * p1 = *poly1;
   ulong * e1 = *exp1;
   slong shift = FLINT_BITS - num*bits;
   TMP_INIT;

   TMP_START;

   prods = (slong *) TMP_ALLOC((num + 1)*sizeof(slong));

   prods[0] = 1;
   for (i = 1; i <= num; i++)
     prods[i] = mults[i - 1]*prods[i - 1];
   
   for (i = 0; i < prods[num]; i++)
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
         
         for (j = 0; j < num; j++)
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
               fmpz * poly2, const slong * mults, slong num, slong bits, slong k)
{
   slong i, j;
   ulong exp;
   fmpz * c;
   slong * prods;
   fmpz * p1 = *poly1;
   ulong * e1 = *exp1;
   slong shift = FLINT_BITS - num*bits;
   TMP_INIT;

   TMP_START;

   prods = (slong *) TMP_ALLOC((num + 1)*sizeof(slong));

   prods[0] = 1;
   for (i = 1; i <= num; i++)
      prods[i] = mults[i - 1]*prods[i - 1];

   for (i = 0; i < prods[num]; i++)
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
         
         for (j = 0; j < num; j++)
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

/*
   use array multiplication to set poly1 to poly2*poly3 in num + 1 variables,
   given a list of multipliers to tightly pack exponents and a number of bits
   for the fields of the exponents of the result, assuming no aliasing
   classical multiplication in main variable, array multiplication for
   multivariate coefficients
*/
slong _fmpz_mpoly_mul_array_chunked(fmpz ** poly1, ulong ** exp1,
        slong * alloc, const fmpz * poly2, const ulong * exp2, slong len2,
                       const fmpz * poly3, const ulong * exp3, slong len3, 
                                          slong * mults, slong num, slong bits)
{
   slong i, j, k = 0, len, l1, l2, l3, prod, bits1, bits2 = 0, bits3 = 0;
   slong shift = FLINT_BITS - bits;
   slong * i2, * i3, * n2, * n3, * b2, * maxb2, * b3, * maxb3;
   ulong * e2, * e3, * p1;
   ulong hi, lo;
   int small;
   TMP_INIT;

   prod = 1;
   for (i = 0; i < num; i++)
      prod *= mults[i];

   l2 = 1 + (slong) (exp2[len2 - 1] >> shift);
   l3 = 1 + (slong) (exp3[len3 - 1] >> shift);

   TMP_START;

   /* compute indices and lengths of coefficients of polys in main variable */

   i2 = (slong *) TMP_ALLOC(l2*sizeof(slong));
   n2 = (slong *) TMP_ALLOC(l2*sizeof(slong));
   b2 = (slong *) TMP_ALLOC(l2*sizeof(slong));
   i3 = (slong *) TMP_ALLOC(l3*sizeof(slong));
   n3 = (slong *) TMP_ALLOC(l3*sizeof(slong));
   b3 = (slong *) TMP_ALLOC(l3*sizeof(slong));
   maxb2 = (slong *) TMP_ALLOC(l2*sizeof(slong));
   maxb3 = (slong *) TMP_ALLOC(l3*sizeof(slong));
   
   i2[0] = 0;
   j = 0;
   for (i = 0; i < l2 - 1; i++)
   {
      while (j < len2 && i == (slong) (exp2[j] >> shift))
         j++;

      i2[i + 1] = j;
      n2[i] = j - i2[i];
   }
   n2[l2 - 1] = len2 - j;

   i3[0] = 0;
   j = 0;
   for (i = 0; i < l3 - 1; i++)
   {
      while (j < len3 && i == (slong) (exp3[j] >> shift))
         j++;

      i3[i + 1] = j;
      n3[i] = j - i3[i];
   }
   n3[l3 - 1] = len3 - j;

   /* pack input exponents tightly */

   e2 = (ulong *) TMP_ALLOC(len2*sizeof(ulong));
   e3 = (ulong *) TMP_ALLOC(len3*sizeof(ulong));

   mpoly_pack_monomials_tight(e2, exp2, len2, mults, num, 1, bits);

   mpoly_pack_monomials_tight(e3, exp3, len3, mults, num, 1, bits);

   /* work out max bits for each coeff and optimal bits */

   for (i = 0; i < l2; i++)
   {
      hi = lo = 0;

      maxb2[i] = 0;

      for (j = 0; j < n2[i]; j++)
      {
         slong bit2 = fmpz_get_si(poly2 + i2[i] + j);
         ulong u2 = (ulong) FLINT_ABS(bit2);

         if (FLINT_BIT_COUNT(u2) > maxb2[i])
            maxb2[i] = FLINT_BIT_COUNT(u2);

         add_ssaaaa(hi, lo, hi, lo, UWORD(0), u2);
      }

      if (bits2 < maxb2[i])
         bits2 = maxb2[i];

      if (hi != 0)
         b2[i] = FLINT_BIT_COUNT(hi) + FLINT_BITS;
      else
         b2[i] = FLINT_BIT_COUNT(lo);
   }

   for (i = 0; i < l3; i++)
   {
      hi = lo = 0;

      maxb3[i] = 0;

      for (j = 0; j < n3[i]; j++)
      {
         slong bit3 = fmpz_get_si(poly3 + i3[i] + j);
         ulong u3 = (ulong) FLINT_ABS(bit3);

         if (FLINT_BIT_COUNT(u3) > maxb3[i])
            maxb3[i] = FLINT_BIT_COUNT(u3);

         add_ssaaaa(hi, lo, hi, lo, UWORD(0), u3);
      }

      if (bits3 < maxb3[i])
         bits3 = maxb3[i];

      if (hi != 0)
         b3[i] = FLINT_BIT_COUNT(hi) + FLINT_BITS;
      else
         b3[i] = FLINT_BIT_COUNT(lo);
   }

   small = bits2 <= (FLINT_BITS - 2) &&
           bits3 <= (FLINT_BITS - 2);

   /* classical multiplication one output coeff at a time */

   l1 = l2 + l3 - 1;

   if (small) /* not used in large case */
   {
      p1 = (ulong *) TMP_ALLOC(3*prod*sizeof(ulong));

      for (i = 0; i < l1; i++)
      {
         slong num1 = 0;
         bits1 = 0;

         for (j = 0; j < l2 && j <= i; j++)
         {
            if (i - j < l3)
            {
               bits1 = FLINT_MAX(bits1, FLINT_MIN(b2[j] + maxb3[i - j], b3[j] + maxb2[i - j]));
               num1++;
            }
         }

         bits1 += FLINT_BIT_COUNT(num1) + 1; /* includes one bit for sign */

         if (bits1 <= FLINT_BITS)
         {
            for (j = 0; j < prod; j++)
               p1[j] = 0;

            for (j = 0; j < l2 && j <= i; j++)
            {
               if (i - j < l3)
               {
                  _fmpz_mpoly_addmul_array1_slong1(p1, 
                     (slong *) poly2 + i2[j], e2 + i2[j], n2[j],
                     (slong *) poly3 + i3[i - j], e3 + i3[i - j], n3[i - j]);

               }
            }

            len = _fmpz_mpoly_from_ulong_array1(poly1, exp1, alloc, 
                                                      p1, mults, num, bits, k) - k;

            for (j = 0; j < len; j++)
               (*exp1)[k + j] = ((*exp1)[k + j] >> bits) + (i << shift);

            k += len;
         } else if (bits1 <= 2*FLINT_BITS)
         {
            for (j = 0; j < 2*prod; j++)
               p1[j] = 0;

            for (j = 0; j < l2 && j <= i; j++)
            {
               if (i - j < l3)
               {
                  _fmpz_mpoly_addmul_array1_slong2(p1, 
                     (slong *) poly2 + i2[j], e2 + i2[j], n2[j],
                     (slong *) poly3 + i3[i - j], e3 + i3[i - j], n3[i - j]);

               }
            }

            len = _fmpz_mpoly_from_ulong_array2(poly1, exp1, alloc, 
                                                      p1, mults, num, bits, k) - k;

            for (j = 0; j < len; j++)
               (*exp1)[k + j] = ((*exp1)[k + j] >> bits) + (i << shift);

            k += len;
         } else
         {
            for (j = 0; j < 3*prod; j++)
               p1[j] = 0;

            for (j = 0; j < l2 && j <= i; j++)
            {
               if (i - j < l3)
               {
                  _fmpz_mpoly_addmul_array1_slong(p1, 
                     (slong *) poly2 + i2[j], e2 + i2[j], n2[j],
                     (slong *) poly3 + i3[i - j], e3 + i3[i - j], n3[i - j]);

               }
            }

            len = _fmpz_mpoly_from_ulong_array(poly1, exp1, alloc, 
                                                      p1, mults, num, bits, k) - k;

            for (j = 0; j < len; j++)
               (*exp1)[k + j] = ((*exp1)[k + j] >> bits) + (i << shift);

            k += len;
         }
      }
   } else
   {
      fmpz * p1 = (fmpz *) TMP_ALLOC(prod*sizeof(fmpz));
     
      for (i = 0; i < l1; i++)
      {
         for (j = 0; j < prod; j++)
            p1[j] = 0;

         for (j = 0; j < l2 && j <= i; j++)
         {
            if (i - j < l3)
            {
               _fmpz_mpoly_addmul_array1_fmpz(p1, 
                   poly2 + i2[j], e2 + i2[j], n2[j],
                   poly3 + i3[i - j], e3 + i3[i - j], n3[i - j]);

             }
          }

          len = _fmpz_mpoly_from_fmpz_array(poly1, exp1, alloc, 
                                                      p1, mults, num, bits, k) - k;

          for (j = 0; j < len; j++)
             (*exp1)[k + j] = ((*exp1)[k + j] >> bits) + (i << shift);

          k += len;
       }
   }

   TMP_END;

   return k;
}

/*
   use array multiplication to set poly1 to poly2*poly3 in num variables, given
   a list of multipliers to tightly pack exponents and a number of bits for the
   fields of the exponents of the result, assuming no aliasing
*/
slong _fmpz_mpoly_mul_array(fmpz ** poly1, ulong ** exp1, slong * alloc,
                         const fmpz * poly2, const ulong * exp2, slong len2, 
                         const fmpz * poly3, const ulong * exp3, slong len3, 
                                          slong * mults, slong num, slong bits)
{
   slong i, bits1, bits2, bits3;
   ulong * e2, * e3;
   slong prod, len;
   int small, sign;
   ulong hi = 0, lo = 0; /* two words to count bits */
   TMP_INIT;

   prod = 1;
   for (i = 0; i < num; i++)
      prod *= mults[i];

   if (prod > MAX_ARRAY_SIZE)
      return _fmpz_mpoly_mul_array_chunked(poly1, exp1, alloc,
                   poly2, exp2, len2, poly3, exp3, len3, mults, num - 1, bits);

   TMP_START;

   e2 = (ulong *) TMP_ALLOC(len2*sizeof(ulong));
   e3 = (ulong *) TMP_ALLOC(len3*sizeof(ulong));

   mpoly_pack_monomials_tight(e2, exp2, len2, mults, num, 0, bits);

   mpoly_pack_monomials_tight(e3, exp3, len3, mults, num, 0, bits);

   bits2 = _fmpz_vec_max_bits(poly2, len2);
   bits3 = _fmpz_vec_max_bits(poly3, len3);

   sign = (bits2 < 0) || (bits3 < 0);

   small = FLINT_ABS(bits2) <= (FLINT_BITS - 2) &&
           FLINT_ABS(bits3) <= (FLINT_BITS - 2);

   bits1 = -1; /* not used in large case */

   if (small)
   {
      if (len2 < len3)
      {
         for (i = 0; i < len2; i++)
         {
            slong b2 = fmpz_get_si(poly2 + i);

            add_ssaaaa(hi, lo, hi, lo, UWORD(0), (ulong) FLINT_ABS(b2));
         }

         bits1 = FLINT_ABS(bits3) + sign;
      } else
      {
         for (i = 0; i < len3; i++)
         {
            slong b3 = fmpz_get_si(poly3 + i);

            add_ssaaaa(hi, lo, hi, lo, UWORD(0), (ulong) FLINT_ABS(b3));
         }

         bits1 = FLINT_ABS(bits2) + sign;
      }

      if (hi != 0)
         bits1 += FLINT_BIT_COUNT(hi) + FLINT_BITS;
      else
         bits1 += FLINT_BIT_COUNT(lo);
   } 

   if (small && bits1 <= FLINT_BITS)
   {
      ulong * p1 = (ulong *) TMP_ALLOC(prod*sizeof(ulong));

      for (i = 0; i < prod; i++)
         p1[i] = 0;

      _fmpz_mpoly_addmul_array1_slong1(p1, 
                         (slong *) poly2, e2, len2, (slong *) poly3, e3, len3);

      len = _fmpz_mpoly_from_ulong_array1(poly1, exp1, alloc, 
                                                      p1, mults, num, bits, 0);
   } else if (small && bits1 <= 2*FLINT_BITS)
   {
      ulong * p1 = (ulong *) TMP_ALLOC(2*prod*sizeof(ulong));

      for (i = 0; i < 2*prod; i++)
         p1[i] = 0;

      _fmpz_mpoly_addmul_array1_slong2(p1, 
                         (slong *) poly2, e2, len2, (slong *) poly3, e3, len3);

      len = _fmpz_mpoly_from_ulong_array2(poly1, exp1, alloc, 
                                                      p1, mults, num, bits, 0);
   } else if (small)
   {
      ulong * p1 = (ulong *) TMP_ALLOC(3*prod*sizeof(ulong));

      for (i = 0; i < 3*prod; i++)
         p1[i] = 0;

      _fmpz_mpoly_addmul_array1_slong(p1, 
                         (slong *) poly2, e2, len2, (slong *) poly3, e3, len3);

      len = _fmpz_mpoly_from_ulong_array(poly1, exp1, alloc, 
                                                      p1, mults, num, bits, 0);      
   } else
   {
      fmpz * p1 = (fmpz *) TMP_ALLOC(prod*sizeof(fmpz));

      for (i = 0; i < prod; i++)
         p1[i] = 0;

      _fmpz_mpoly_addmul_array1_fmpz(p1, poly2, e2, len2, poly3, e3, len3);

      len = _fmpz_mpoly_from_fmpz_array(poly1, exp1, alloc, 
                                                      p1, mults, num, bits, 0);
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
   ulong * exp2, * exp3;
   int free2, free3;
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
   for (i = 0; i < ctx->n - 1; i++)
   {
      max_degs2[i] += max_degs3[i] + 1;
      array_size *= max_degs2[i];  
   }
   max_degs2[ctx->n - 1] += max_degs3[ctx->n - 1] + 1;

   if (array_size > MAX_ARRAY_SIZE)
   {
      res = 0;
      goto cleanup;
   }

   exp2 = mpoly_unpack_monomials(exp_bits, poly2->exps, 
                                           poly2->length, ctx->n, poly2->bits);

   free2 = exp2 != poly2->exps;

   exp3 = mpoly_unpack_monomials(exp_bits, poly3->exps, 
                                           poly3->length, ctx->n, poly3->bits);
   
   free3 = exp3 != poly3->exps;

   if (poly1 == poly2 || poly1 == poly3)
   {
      fmpz_mpoly_t temp;

      fmpz_mpoly_init2(temp, poly2->length + poly3->length - 1, ctx);
      fmpz_mpoly_fit_bits(temp, exp_bits, ctx);

      if (poly2->length >= poly3->length)
         len = _fmpz_mpoly_mul_array(&temp->coeffs, &temp->exps, &temp->alloc, 
                                           poly3->coeffs, exp3, poly3->length,
                                           poly2->coeffs, exp2, poly2->length,
                                        (slong *) max_degs2, ctx->n, exp_bits);
      else
         len = _fmpz_mpoly_mul_array(&temp->coeffs, &temp->exps, &temp->alloc, 
                                           poly2->coeffs, exp2, poly2->length,
                                           poly3->coeffs, exp3, poly3->length,
                                        (slong *) max_degs2, ctx->n, exp_bits);

      fmpz_mpoly_swap(temp, poly1, ctx);

      fmpz_mpoly_clear(temp, ctx);
   } else
   {
      fmpz_mpoly_fit_length(poly1, poly2->length + poly3->length - 1, ctx);
      fmpz_mpoly_fit_bits(poly1, exp_bits, ctx);

      if (poly2->length >= poly3->length)
         len = _fmpz_mpoly_mul_array(&poly1->coeffs, &poly1->exps, &poly1->alloc,
                                            poly3->coeffs, exp3, poly3->length,
                                            poly2->coeffs, exp2, poly2->length,
                                        (slong *) max_degs2, ctx->n, exp_bits);
      else
         len = _fmpz_mpoly_mul_array(&poly1->coeffs, &poly1->exps, &poly1->alloc,
                                            poly2->coeffs, exp2, poly2->length, 
                                            poly3->coeffs, exp3, poly3->length,
                                        (slong *) max_degs2, ctx->n, exp_bits);
   }

   _fmpz_mpoly_set_length(poly1, len, ctx);

   if (free2)
      flint_free(exp2);

   if (free3)
      flint_free(exp3);

cleanup:

   TMP_END;

   return res;
}
