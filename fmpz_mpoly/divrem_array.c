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

slong _fmpz_mpoly_divrem_array_tight(slong * lenr,
 fmpz ** polyq, ulong ** expq, slong * allocq, slong len0,
       fmpz ** polyr, ulong ** expr, slong * allocr, slong len1,
                  const fmpz * poly2, const ulong * exp2, slong len2,
                        const fmpz * poly3, const ulong * exp3, slong len3,
                                          slong * mults, slong num, slong bits)
{
   slong i, j, q, r, prod, bits1, bits2, bits3, k = len0, l = len1;
   slong max3 = (slong) exp3[len3 - 1]; /* largest exponent in poly3 */
   slong min3 = (slong) exp3[0]; /* smallest exponent in poly3 */
   slong * prods;
   fmpz c3 = poly3[0];
   ulong u3 = ((ulong) FLINT_ABS(c3)) >> 1;
   fmpz * p1 = *polyq, * p2 = *polyr;
   ulong * e1 = *expq, * e2 = *expr;
   int small;
   TMP_INIT;

   TMP_START;

   prods = (slong *) TMP_ALLOC((num + 1)*sizeof(slong));

   prods[0] = 1;
   for (i = 1; i <= num; i++)
      prods[i] = mults[i - 1]*prods[i - 1];

   prod = prods[num];

   bits2 = _fmpz_vec_max_bits(poly2, len2);
   bits3 = _fmpz_vec_max_bits(poly3, len3);
   /* we assume a bound of FLINT_BITS - 2 for coefficients of the quotient */
   bits1 = FLINT_ABS(bits3) + FLINT_BITS + FLINT_BIT_COUNT(len3) - 2;

   small = FLINT_ABS(bits2) <= bits1 && FLINT_ABS(bits3) <= FLINT_BITS - 2;
   bits1 += 2; /* incr. so poly2 - q*poly3 doesn't overflow and for sign */

   if (small && bits1 <= 2*FLINT_BITS)
   {
      ulong * t2 = (ulong *) TMP_ALLOC(2*prod*sizeof(ulong));

      for (i = 0; i < 2*prod; i++)
         t2[i] = 0;

      _fmpz_mpoly_to_ulong_array2(t2, poly2, exp2, len2);

      for (i = prod - 1; i >= max3; i--)
      {
         ulong * ptr = t2 + 2*i;
         ulong p[2];

         if (ptr[0] != 0 || ptr[1] != 0)
         {
            if (0 > (slong) ptr[1])
               mpn_neg(p, ptr, 2);
            else
               flint_mpn_copyi(p, ptr, 2);

            /* remainder monomial */
            if (!mpoly_monomial_divides_tight(i, max3, prods, num))
            {
               if (l >= *allocr)
               {
                  p2 = (fmpz *) flint_realloc(p2, 2*(*allocr)*sizeof(fmpz));
                  e2 = (ulong *) flint_realloc(e2, 2*(*allocr)*sizeof(ulong));
                  flint_mpn_zero(p2 + *allocr, *allocr);
                  (*allocr) *= 2;
               }

               fmpz_set_ui(p2 + l, p[1]);
               fmpz_mul_2exp(p2 + l, p2 + l, FLINT_BITS);
               fmpz_add_ui(p2 + l, p2 + l, p[0]);
               if (0 > (slong) ptr[1])
                  fmpz_neg(p2 + l, p2 + l);

               e2[l++] = i;
            } else /* monomials can be divided exactly */
            {
               if (u3 < p[1] || (u3 == 0 && 0 > (slong) p[0])) /* quotient too large */
               {
                  for (j = len0; j < k; j++)
                     _fmpz_demote(p1 + j);
                  for (j = len1; j < l; j++)
                     _fmpz_demote(p2 + j);
                  k = len0;
                  l = len1;

                  goto big;
               }

               sdiv_qrnnd(q, r, ptr[1], ptr[0], c3);

               if (COEFF_IS_MPZ(q)) /* quotient too large */
               {
                  for (j = len0; j < k; j++)
                     _fmpz_demote(p1 + j);
                  for (j = len1; j < l; j++)
                     _fmpz_demote(p2 + j);
                  k = len0;
                  l = len1;

                  goto big;
               }

               if (r != 0) /* not an exact division */
               {
                  if (l >= *allocr)
                  {
                     p2 = (fmpz *) flint_realloc(p2, 2*(*allocr)*sizeof(fmpz));
                     e2 = (ulong *) flint_realloc(e2, 2*(*allocr)*sizeof(ulong));
                     flint_mpn_zero(p2 + *allocr, *allocr);
                     (*allocr) *= 2;
                  }

                  fmpz_set_si(p2 + l, (slong) r);  

                  e2[l++] = i;
               }

               _fmpz_mpoly_submul_array1_slong2_1(t2, q, i - max3 + min3,
                                                            poly3, exp3, len3);

               if (k >= *allocq)
               {
                  p1 = (fmpz *) flint_realloc(p1, 2*(*allocq)*sizeof(fmpz));
                  e1 = (ulong *) flint_realloc(e1, 2*(*allocq)*sizeof(ulong));
                  flint_mpn_zero(p1 + *allocq, *allocq);
                  (*allocq) *= 2;
               }
            
               p1[k] = q;
               e1[k++] = i - max3;
            }         
         }
      }

      for ( ; i >= 0; i--)
      {
         ulong * ptr = t2 + 2*i;

         if (ptr[0] != 0 || ptr[1] != 0)  /* not an exact division */
         {
            if (l >= *allocr)
            {
               p2 = (fmpz *) flint_realloc(p2, 2*(*allocr)*sizeof(fmpz));
               e2 = (ulong *) flint_realloc(e2, 2*(*allocr)*sizeof(ulong));
               flint_mpn_zero(p2 + *allocr, *allocr);
               (*allocr) *= 2;
            }

            fmpz_set_ui(p2 + l, ptr[1]);
            fmpz_mul_2exp(p2 + l, p2 + l, FLINT_BITS);
            fmpz_add_ui(p2 + l, p2 + l, ptr[0]);
            if (0 > (slong) ptr[1])
               fmpz_neg(p2 + l, p2 + l);

            e2[l++] = i;
         }
      }
   }

   if (k == len0 && l == len1 && small)
   {
      ulong * t2 = (ulong *) TMP_ALLOC(3*prod*sizeof(ulong));

      for (i = 0; i < 3*prod; i++)
         t2[i] = 0;

      _fmpz_mpoly_to_ulong_array(t2, poly2, exp2, len2);

      for (i = prod - 1; i >= max3; i--)
      {
         ulong * ptr = t2 + 3*i;
         ulong p[3];

         if (ptr[0] != 0 || ptr[1] != 0 || ptr[2] != 0)
         {
            if (0 > (slong) ptr[2])
               mpn_neg(p, ptr, 3);
            else
               flint_mpn_copyi(p, ptr, 3);

            /* remainder monomial */
            if (!mpoly_monomial_divides_tight(i, max3, prods, num))
            {
               if (l >= *allocr)
               {
                  p2 = (fmpz *) flint_realloc(p2, 2*(*allocr)*sizeof(fmpz));
                  e2 = (ulong *) flint_realloc(e2, 2*(*allocr)*sizeof(ulong));
                  flint_mpn_zero(p2 + *allocr, *allocr);
                  (*allocr) *= 2;
               }

               fmpz_set_ui(p2 + l, p[2]);
               fmpz_mul_2exp(p2 + l, p2 + l, FLINT_BITS);
               fmpz_add_ui(p2 + l, p2 + l, p[1]);
               fmpz_mul_2exp(p2 + l, p2 + l, FLINT_BITS);
               fmpz_add_ui(p2 + l, p2 + l, p[0]);
               if (0 > (slong) ptr[2])
                  fmpz_neg(p2 + l, p2 + l);

               e2[l++] = i;
            } else /* monomials can be divided exact */
            {
               if (p[2] > 0 || u3 < p[1] || (u3 == 0 && 0 > (slong) p[0])) /* quotient too large */
               {
                  for (j = len0; j < k; j++)
                     _fmpz_demote(p1 + j);
                  for (j = len1; j < l; j++)
                     _fmpz_demote(p2 + j);
                  k = len0;
                  l = len1;

                  goto big;
               }

               sdiv_qrnnd(q, r, ptr[1], ptr[0], c3);

               if (COEFF_IS_MPZ(q)) /* quotient too large */
               {
                  for (j = len0; j < k; j++)
                     _fmpz_demote(p1 + j);
                  for (j = len1; j < l; j++)
                     _fmpz_demote(p2 + j);
                  k = len0;
                  l = len1;

                  goto big;
               }

               if (r != 0) /* remainder */ 
               {
                  if (l >= *allocr)
                  {
                     p2 = (fmpz *) flint_realloc(p2, 2*(*allocr)*sizeof(fmpz));
                     e2 = (ulong *) flint_realloc(e2, 2*(*allocr)*sizeof(ulong));
                     flint_mpn_zero(p2 + *allocr, *allocr);
                     (*allocr) *= 2;
                  }

                  fmpz_set_si(p2 + l, (slong) r);  

                  e2[l++] = i;
               }

               _fmpz_mpoly_submul_array1_slong_1(t2, q, i - max3 + min3,
                                                            poly3, exp3, len3);

               if (k >= *allocq)
               {
                  p1 = (fmpz *) flint_realloc(p1, 2*(*allocq)*sizeof(fmpz));
                  e1 = (ulong *) flint_realloc(e1, 2*(*allocq)*sizeof(ulong));
                  flint_mpn_zero(p1 + *allocq, *allocq);
                  (*allocq) *= 2;
               }
            
               p1[k] = q;
               e1[k++] = i - max3;
            }
         }
      }

      for ( ; i >= 0; i--)
      {
         ulong * ptr = t2 + 3*i;

         if (ptr[0] != 0 || ptr[1] != 0 || ptr[2] != 0) 
         {
            /* not an exact division */

            if (l >= *allocr)
            {
               p2 = (fmpz *) flint_realloc(p2, 2*(*allocr)*sizeof(fmpz));
               e2 = (ulong *) flint_realloc(e2, 2*(*allocr)*sizeof(ulong));
               flint_mpn_zero(p2 + *allocr, *allocr);
               (*allocr) *= 2;
            }

            fmpz_set_ui(p2 + l, ptr[2]);
            fmpz_mul_2exp(p2 + l, p2 + l, FLINT_BITS);
            fmpz_add_ui(p2 + l, p2 + l, ptr[1]);
            fmpz_mul_2exp(p2 + l, p2 + l, FLINT_BITS);
            fmpz_add_ui(p2 + l, p2 + l, ptr[0]);
            if (0 > (slong) ptr[2])
               fmpz_neg(p2 + l, p2 + l);

            e2[l++] = i;
         }
      }
   }

big:

   if (k == len0 && l == len1)
   {
      fmpz * t2 = (fmpz *) TMP_ALLOC(prod*sizeof(fmpz));
      fmpz_t fq, fr;

      fmpz_init(fq);
      fmpz_init(fr);

      for (i = 0; i < prod; i++)
         fmpz_init(t2 + i);

      _fmpz_mpoly_to_fmpz_array(t2, poly2, exp2, len2);
      
      for (i = prod - 1; i >= max3; i--)
      {
         if (!fmpz_is_zero(t2 + i))
         {
            /* remainder monomial */
            if (!mpoly_monomial_divides_tight(i, max3, prods, num))
            {
               if (l >= *allocr)
               {
                  p2 = (fmpz *) flint_realloc(p2, 2*(*allocr)*sizeof(fmpz));
                  e2 = (ulong *) flint_realloc(e2, 2*(*allocr)*sizeof(ulong));
                  flint_mpn_zero(p2 + *allocr, *allocr);
                  (*allocr) *= 2;
               }

               fmpz_set(p2 + l, t2 + i);
               
               e2[l++] = i;
            } else /* monomials can be divided exactly */
            {
               fmpz_fdiv_qr(fq, fr, t2 + i, poly2 + len2 - 1);

               if (!fmpz_is_zero(fr)) /* not an exact division */
               {
                  if (l >= *allocr)
                  {
                     p2 = (fmpz *) flint_realloc(p2, 2*(*allocr)*sizeof(fmpz));
                     e2 = (ulong *) flint_realloc(e2, 2*(*allocr)*sizeof(ulong));
                     flint_mpn_zero(p2 + *allocr, *allocr);
                     (*allocr) *= 2;
                  }

                  fmpz_set(p2 + l, fr);  

                  e2[l++] = i;
               }

               _fmpz_mpoly_submul_array1_fmpz_1(t2, fq, i - min3,
                                                            poly3, exp3, len3);

               if (k >= *allocq)
               {
                  p1 = (fmpz *) flint_realloc(p1, 2*(*allocq)*sizeof(fmpz));
                  e1 = (ulong *) flint_realloc(e1, 2*(*allocq)*sizeof(ulong));
                  flint_mpn_zero(p1 + *allocq, *allocq);
                  (*allocq) *= 2;
               }
            
               fmpz_set(p1 + k, fq);
               e1[k++] = i - max3;
            }
         }
      }

      for ( ; i >= 0; i--)
      {
         if (!fmpz_is_zero(t2 + i))
         {
            /* remainder */

            if (l >= *allocr)
            {
               p2 = (fmpz *) flint_realloc(p2, 2*(*allocr)*sizeof(fmpz));
               e2 = (ulong *) flint_realloc(e2, 2*(*allocr)*sizeof(ulong));
               flint_mpn_zero(p2 + *allocr, *allocr);
               (*allocr) *= 2;
            }

            fmpz_set(p2 + l, t2 + i);
            e2[l++] = i;
         }
      }

      fmpz_clear(fq);
      fmpz_clear(fr);

      for (i = 0; i < prod; i++)
         fmpz_clear(t2 + i);
   }

   (*polyq) = p1;
   (*expq) = e1;
   (*polyr) = p2;
   (*expr) = e2;

   (*lenr) = l - len1;

   TMP_END;

   return k - len0;
}

/*
   use array division to set q, r to poly2/poly3 in num + 1 variables,
   given a list of multipliers to tightly pack exponents and a number of bits
   for the fields of the exponents of the result, assuming no aliasing
   classical exact division in main variable, array multiplication (submul)
   for multivariate coefficients
*/
slong _fmpz_mpoly_divrem_array_chunked(slong * lenr,
            fmpz ** polyq, ulong ** expq, slong * allocq,
                 fmpz ** polyr, ulong ** expr, slong * allocr, 
                const fmpz * poly2, const ulong * exp2, slong len2, 
        const fmpz * poly3, const ulong * exp3, slong len3, slong * mults, 
                                                         slong num, slong bits)
{
   slong i, j, k, l = 0, prod, len = 0, l1, l2, l3;
   slong bits1, bits2, bits3 = 0, tlen, talloc;
   slong shift = FLINT_BITS - bits;
   slong * i1, * i2, * i3, * n1, * n2, * n3;
   slong * b1, * b3, * maxb1, * maxb3;
   ulong * e2, * e3, * texp, * p2;
   fmpz * temp;
   int small;
   TMP_INIT;

   prod = 1;
   for (i = 0; i < num; i++)
      prod *= mults[i];

   l2 = 1 + (slong) (exp2[len2 - 1] >> shift);
   l3 = 1 + (slong) (exp3[len3 - 1] >> shift);

   l1 = l2 - l3 + 1;

   TMP_START;

   /* compute indices and lengths of coefficients of polys in main variable */

   i1 = (slong *) TMP_ALLOC(l1*sizeof(slong));
   n1 = (slong *) TMP_ALLOC(l1*sizeof(slong));
   i2 = (slong *) TMP_ALLOC(l2*sizeof(slong));
   n2 = (slong *) TMP_ALLOC(l2*sizeof(slong));
   i3 = (slong *) TMP_ALLOC(l3*sizeof(slong));
   n3 = (slong *) TMP_ALLOC(l3*sizeof(slong));
   b1 = (slong *) TMP_ALLOC(l1*sizeof(slong));
   maxb1 = (slong *) TMP_ALLOC(l1*sizeof(slong));
   b3 = (slong *) TMP_ALLOC(l3*sizeof(slong));
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

   /* work out max bits for each coeff and optimal bits */

   for (i = 0; i < l3; i++)
   {
      _fmpz_mpoly_chunk_max_bits(b3, maxb3, poly3, i3, n3, i);

      if (bits3 < maxb3[i])
         bits3 = maxb3[i];
   }

   /* pack input coefficients tightly */

   e2 = (ulong *) TMP_ALLOC(len2*sizeof(ulong));
   e3 = (ulong *) TMP_ALLOC(len3*sizeof(ulong));

   mpoly_pack_monomials_tight(e2, exp2, len2, mults, num, 1, bits);

   mpoly_pack_monomials_tight(e3, exp3, len3, mults, num, 1, bits);

   bits2 = _fmpz_vec_max_bits(poly2, len2);
   /* we assume a bound of FLINT_BITS - 2 for coefficients of the quotient */
   bits1 = FLINT_ABS(bits3) + FLINT_BITS + FLINT_BIT_COUNT(len3) - 2;

   small = FLINT_ABS(bits2) <= bits1 && FLINT_ABS(bits3) <= FLINT_BITS - 2;

   /* space for copy of leading "coefficient" of poly2 */

   temp = (fmpz *) flint_calloc(n2[l2 - 1] + 1, sizeof(fmpz));
   texp = (ulong *) flint_malloc((n2[l2 - 1] + 1)*sizeof(ulong));
   talloc = n2[l2 - 1] + 1; /* plus one so doubling always increases size */

   p2 = (ulong *) TMP_ALLOC(3*prod*sizeof(ulong));

   if (small)
   {
      for (i = l2 - 1; i >= 0; i--)
      {
         slong num1 = 0;
         bits1 = 0;

         if (i != l2 - 1)
         {
            for (j = 0; j < l2 - i - 1 && j < l1; j++)
            {
               k = i - l1 + j + 1;

               if (k < l3 && k >= 0)
               {
                  bits1 = FLINT_MAX(bits1,
                                FLINT_MIN(b1[j] + maxb3[k], maxb1[j] + b3[k]));
                  num1++;
               }
            }

            bits1 += FLINT_BIT_COUNT(num1);
            bits1 = FLINT_MAX(FLINT_ABS(bits2), bits1);

            bits1 += 2; /* bit for sign and so a - q*b doesn't overflow */
         } else
            bits1 = FLINT_ABS(bits2) + 1; /* extra bit for sign */

         if (bits1 <= FLINT_BITS)
         {
            for (j = 0; j < prod; j++)
               p2[j] = 0;

            _fmpz_mpoly_to_ulong_array1(p2, poly2 + i2[i], e2 + i2[i], n2[i]);

            /* submuls */

            for (j = 0; j < l2 - i - 1 && j < l1; j++)
            {
               k = i - l1 + j + 1;

               if (k < l3 && k >= 0)
               {
                  _fmpz_mpoly_submul_array1_slong1(p2, (*polyq) + i1[j],
                     (*expq) + i1[j], n1[j], poly3 + i3[k], e3 + i3[k], n3[k]);
               }
            }

            tlen = _fmpz_mpoly_from_ulong_array1(&temp, &texp, &talloc, 
                                                      p2, mults, num, bits, 0);
         } else if (bits1 <= 2*FLINT_BITS)
         {
            for (j = 0; j < 2*prod; j++)
               p2[j] = 0;

            _fmpz_mpoly_to_ulong_array2(p2, poly2 + i2[i], e2 + i2[i], n2[i]);

            /* submuls */

            for (j = 0; j < l2 - i - 1 && j < l1; j++)
            {
               k = i - l1 + j + 1;

               if (k < l3 && k >= 0)
               {
                  _fmpz_mpoly_submul_array1_slong2(p2, (*polyq) + i1[j],
                     (*expq) + i1[j], n1[j], poly3 + i3[k], e3 + i3[k], n3[k]);
               }
            }

            tlen = _fmpz_mpoly_from_ulong_array2(&temp, &texp, &talloc, 
                                                      p2, mults, num, bits, 0);
         } else /* <= 3*FLINT_BITS */
         {
            for (j = 0; j < 3*prod; j++)
               p2[j] = 0;

            _fmpz_mpoly_to_ulong_array(p2, poly2 + i2[i], e2 + i2[i], n2[i]);

            /* submuls */

            for (j = 0; j < l2 - i - 1 && j < l1; j++)
            {
               k = i - l1 + j + 1;

               if (k < l3 && k >= 0)
                  _fmpz_mpoly_submul_array1_slong(p2, (*polyq) + i1[j],
                     (*expq) + i1[j], n1[j], poly3 + i3[k], e3 + i3[k], n3[k]);
            }

            tlen = _fmpz_mpoly_from_ulong_array(&temp, &texp, &talloc, 
                                                      p2, mults, num, bits, 0);
         }

         if (tlen != 0) /* nonzero "coefficient" */
         {
            if (l2 - i - 1 < l1) /* potentially a quotient with remainder */
            {
               mpoly_pack_monomials_tight(texp, texp, tlen, mults, num, 0, bits);

               i1[l2 - i - 1] = len;

               n1[l2 - i - 1] = _fmpz_mpoly_divrem_array_tight(lenr, polyq,
                 expq, allocq, len, polyr, expr, allocr, l, temp, texp,
                     tlen, poly3 + i3[l3 - 1], e3 + i3[l3 - 1], n3[l3 - 1],
                                                             mults, num, bits);

               /* check the quotient didn't have large coefficients */
               if (FLINT_ABS(_fmpz_vec_max_bits((*polyq) + len, n1[l2 - i - 1])) >
                                                             FLINT_BITS - 2)
               {
                  for (j = 0; j < len; j++)
                     _fmpz_demote((*polyq) + j);
                  for (j = 0; j < l; j++)
                     _fmpz_demote((*polyr) + j);
                  len = 0;
                  l = 0;

                  goto big;
               }

               _fmpz_mpoly_chunk_max_bits(b1, maxb1, *polyq, i1, n1, l2 - i - 1);

               len += n1[l2 - i - 1];
               l += *lenr;
            } else /* remainder terms only */
            {
               if (l + tlen > *allocr)
               {
                  slong alloc = FLINT_MAX(l + tlen, 2*(*allocr));

                  (*polyr) = (fmpz *) flint_realloc(*polyr, alloc*sizeof(fmpz));
                  (*expr) = (ulong *) flint_realloc(*expr, alloc*sizeof(ulong));
                  flint_mpn_zero(*polyr + *allocr, alloc - *allocr);
                  (*allocr) = alloc;
               }

               for (j = 0; j < tlen; j++)
               {
                  fmpz_set(*polyr + l + j, temp + j);
                  (*expr)[l + j] = texp[j] + (i << shift);
               }

               l += tlen;
            }
         } else if (l2 - i - 1 < l1)
         {
            i1[l2 - i - 1] = len;
            n1[l2 - i - 1] = 0;
         }
      }
   }

big:

   if (len == 0)
   {
      fmpz * p2 = (fmpz *) TMP_ALLOC(prod*sizeof(fmpz));

      for (j = 0; j < prod; j++)
            fmpz_init(p2 + j);
      
      for (i = l2 - 1; i >= 0; i--)
      {
         for (j = 0; j < prod; j++)
            fmpz_zero(p2 + j);

         _fmpz_mpoly_to_fmpz_array(p2, poly2 + i2[i], e2 + i2[i], n2[i]);

         /* submuls */

         for (j = 0; j < l2 - i - 1 && j < l1; j++)
         {
            k = i - l1 + j + 1;

            if (k < l3 && k >= 0)
               _fmpz_mpoly_submul_array1_fmpz(p2, (*polyq) + i1[j],
                    (*expq) + i1[j], n1[j], poly3 + i3[k], e3 + i3[k], n3[k]);
         }

         tlen = _fmpz_mpoly_from_fmpz_array(&temp, &texp, &talloc, 
                                                      p2, mults, num, bits, 0);

         if (tlen != 0) /* nonzero "coefficient" */
         {
            if (l2 - i - 1 < l1) /* potentially a quotient with remainder */
            {
               mpoly_pack_monomials_tight(texp, texp, tlen, mults, num, 0, bits);

               i1[l2 - i - 1] = len;
            
               n1[l2 - i - 1] = _fmpz_mpoly_divrem_array_tight(lenr, polyq,
                      expq, allocq, len, polyr, expr, allocr, l, temp, texp, 
                             tlen, poly3 + i3[l3 - 1], e3 + i3[l3 - 1],
                                                 n3[l3 - 1], mults, num, bits);

               _fmpz_mpoly_chunk_max_bits(b1, maxb1, *polyq, i1, n1, l2 - i - 1);

               len += n1[l2 - i - 1];
               l += *lenr;
            } else /* remainder terms */
            {
               if (l + tlen > *allocr)
               {
                  slong alloc = FLINT_MAX(l + tlen, 2*(*allocr));

                  (*polyr) = (fmpz *) flint_realloc(*polyr, alloc*sizeof(fmpz));
                  (*expr) = (ulong *) flint_realloc(*expr, alloc*sizeof(ulong));
                  flint_mpn_zero(*polyr + *allocr, alloc - *allocr);
                  (*allocr) = alloc;
               }

               for (j = 0; j < tlen; j++)
               {
                  fmpz_set(*polyr + l + j, temp + j);
                  (*expr)[l + j] = texp[j] + (i << shift);
               }

               l += tlen;
            }
         } else if (l2 - i - 1 < l1)
         {
            i1[l2 - i - 1] = len;
            n1[l2 - i - 1] = 0;
         }
      }

      for (j = 0; j < prod; j++)
            fmpz_clear(p2 + j);
   }

   if (len != 0)
   {
      mpoly_unpack_monomials_tight((*expq), (*expq), len, mults, num, 1, bits);

      /* put main variable back in quotient */
      for (i = 0; i < l1; i++)
      {
         for (j = 0; j < n1[l1 - i - 1]; j++)
         {
            (*expq)[i1[l1 - i - 1] + j] += (i << shift);
         }
      }
   }

   flint_free(temp);
   flint_free(texp);

   TMP_END;

   *lenr = l;

   return len;
}

/*
   use array division to set q, r to poly2/poly3 in num variables, given
   a list of multipliers to tightly pack exponents and a number of bits for the
   fields of the exponents of the result, assuming no aliasing
*/
slong _fmpz_mpoly_divrem_array(slong * lenr,
       fmpz ** polyq, ulong ** expq, slong * allocq,
              fmpz ** polyr, ulong ** expr, slong * allocr, 
                const fmpz * poly2, const ulong * exp2, slong len2, 
        const fmpz * poly3, const ulong * exp3, slong len3, slong * mults, 
                                                         slong num, slong bits)
{
   slong i;
   ulong * e2, * e3;
   slong len, prod;
   TMP_INIT;

   prod = 1;
   for (i = 0; i < num; i++)
      prod *= mults[i];

   if (prod > MAX_ARRAY_SIZE)
      return _fmpz_mpoly_divrem_array_chunked(lenr, polyq, expq, allocq,
                                     polyr, expr, allocr, poly2, exp2, len2,
                                      poly3, exp3, len3, mults, num - 1, bits);

   TMP_START;

   e2 = (ulong *) TMP_ALLOC(len2*sizeof(ulong));
   e3 = (ulong *) TMP_ALLOC(len3*sizeof(ulong));

   mpoly_pack_monomials_tight(e2, exp2, len2, mults, num, 0, bits);

   mpoly_pack_monomials_tight(e3, exp3, len3, mults, num, 0, bits);

   len = _fmpz_mpoly_divrem_array_tight(lenr, polyq, expq, allocq, 0,
                                  polyr, expr, allocr, 0, poly2, e2, len2,
                                            poly3, e3, len3, mults, num, bits);

   mpoly_unpack_monomials_tight((*expq), (*expq), len, mults, num, 0, bits);
   mpoly_unpack_monomials_tight((*expr), (*expr), *lenr, mults, num, 0, bits);

   TMP_END;

   return len;
}

int fmpz_mpoly_divrem_array(fmpz_mpoly_t q, fmpz_mpoly_t r,
                    const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3, 
                                                    const fmpz_mpoly_ctx_t ctx)
{
   slong i, exp_bits, N, lenq = 0, lenr = 0, array_size;
   ulong * max_degs2, * max_degs3;
   ulong * exp2, * exp3, * maxexp;
   int free2 = 0, free3 = 0;
   fmpz_mpoly_t temp1, temp2;
   fmpz_mpoly_struct * tq, * tr;
   int res = 0;

   TMP_INIT;

   if (poly3->length == 0)
      flint_throw(FLINT_DIVZERO, "Divide by zero in fmpz_mpoly_divrem_array");

   if (poly2->length == 0)
   {
      fmpz_mpoly_zero(q, ctx);
      fmpz_mpoly_zero(r, ctx);

      return 1;
   }

   TMP_START;

   maxexp = (ulong *) TMP_ALLOC(ctx->n*sizeof(ulong));

   max_degs2 = (ulong *) TMP_ALLOC(ctx->n*sizeof(ulong));
   max_degs3 = (ulong *) TMP_ALLOC(ctx->n*sizeof(ulong));

   fmpz_mpoly_max_degrees(max_degs2, poly2, ctx);
   fmpz_mpoly_max_degrees(max_degs3, poly3, ctx);

   for (i = 0; i < ctx->n; i++)
      maxexp[i] = FLINT_MAX(max_degs2[i], max_degs3[i]);

   exp_bits = FLINT_MAX(poly2->bits, poly3->bits);

   N = (exp_bits*ctx->n - 1)/FLINT_BITS + 1;

   if (N != 1)
      goto cleanup;

   array_size = 1;
   for (i = 0; i < ctx->n - 1; i++)
   {
      max_degs2[i] = maxexp[i] + 1;
      array_size *= max_degs2[i];
   }  
   max_degs2[ctx->n - 1] = maxexp[ctx->n - 1] + 1;
   
   if (array_size > MAX_ARRAY_SIZE)
      goto cleanup;

   exp2 = mpoly_unpack_monomials(exp_bits, poly2->exps, 
                                           poly2->length, ctx->n, poly2->bits);

   free2 = exp2 != poly2->exps;

   exp3 = mpoly_unpack_monomials(exp_bits, poly3->exps, 
                                           poly3->length, ctx->n, poly3->bits);
   
   free3 = exp3 != poly3->exps;

   if (q == poly2 || q == poly3)
   {
      fmpz_mpoly_init2(temp1, FLINT_MAX(poly2->length/poly3->length + 1, 1),
                                                                          ctx);
      fmpz_mpoly_fit_bits(temp1, exp_bits, ctx);

      tq = temp1;
   } else
   {
      fmpz_mpoly_fit_length(q, FLINT_MAX(poly2->length/poly3->length + 1, 1),
                                                                          ctx);
      fmpz_mpoly_fit_bits(q, exp_bits, ctx);

      tq = q;
   }

   if (r == poly2 || r == poly3)
   {
      fmpz_mpoly_init2(temp2, poly3->length, ctx);
      fmpz_mpoly_fit_bits(temp2, exp_bits, ctx);

      tr = temp2;
   } else
   {
      fmpz_mpoly_fit_length(r, poly3->length, ctx);
      fmpz_mpoly_fit_bits(r, exp_bits, ctx);

      tr = r;
   }

   lenq = _fmpz_mpoly_divrem_array(&lenr, &tq->coeffs, &tq->exps,
        &tq->alloc, &tr->coeffs, &tr->exps, &tr->alloc, poly2->coeffs,
                  exp2, poly2->length, poly3->coeffs, exp3, poly3->length,
                                        (slong *) max_degs2, ctx->n, exp_bits);

   res = 1;

   if (q == poly2 || q == poly3)
   {
      fmpz_mpoly_swap(temp1, q, ctx);

      fmpz_mpoly_clear(temp1, ctx);
   } 

   if (r == poly2 || r == poly3)
   {
      fmpz_mpoly_swap(temp2, r, ctx);

      fmpz_mpoly_clear(temp2, ctx);
   } 

   _fmpz_mpoly_set_length(q, lenq, ctx);
   _fmpz_mpoly_set_length(r, lenr, ctx);

   fmpz_mpoly_reverse(q, q, ctx);
   fmpz_mpoly_reverse(r, r, ctx);

   if (free2)
      flint_free(exp2);

   if (free3)
      flint_free(exp3);

cleanup:

   TMP_END;

   return res;
}
