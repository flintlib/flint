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

void _fmpz_mpoly_submul_array1_slong(ulong * poly1, 
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
                  sub_dddmmmsss(cy, c[1], c[0], 0, c[1], c[0], 0, p[1], p[0]);
                  c[2] += (0 <= (slong) p[1]) ? cy : cy - 1;
               }
            }
         }
      }
   }
}

void _fmpz_mpoly_submul_array1_slong2(ulong * poly1, 
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
                  sub_ddmmss(c[1], c[0], c[1], c[0], p[1], p[0]);
               }
            }
         }
      }
   }
}

void _fmpz_mpoly_submul_array1_fmpz(fmpz * poly1, 
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
                  fmpz_submul(c, poly2 + i, poly3 + i);
               }
            }
         }
      }
   }
}

void _fmpz_mpoly_submul_array1_slong_1(ulong * poly1, 
                          slong d, const ulong exp2,
                           const slong * poly3, const ulong * exp3, slong len3)
{
   slong j;
   ulong cy;
   ulong p[2]; /* for products of coefficients */
   ulong * c2, * c;

   c2 = poly1 + 3*((slong) exp2);

   if (d != 0)
   {
      for (j = 0; j < len3; j++)
      {
         c = c2 + 3*((slong) exp3[j]);

         smul_ppmm(p[1], p[0], d, poly3[j]);
         sub_dddmmmsss(cy, c[1], c[0], 0, c[1], c[0], 0, p[1], p[0]);
         c[2] += (0 <= (slong) p[1]) ? cy : cy - 1;
      }
   }
}

void _fmpz_mpoly_submul_array1_slong2_1(ulong * poly1, 
                           slong d, const ulong exp2,
                           const slong * poly3, const ulong * exp3, slong len3)
{
   slong j;
   ulong p[2]; /* for products of coefficients */
   ulong * c2, * c;

   c2 = poly1 + 2*((slong) exp2);

   if (d != 0)
   {
      for (j = 0; j < len3; j++)
      {
         c = c2 + 2*((slong) exp3[j]);

         smul_ppmm(p[1], p[0], d, poly3[j]);
         sub_ddmmss(c[1], c[0], c[1], c[0], p[1], p[0]);
      }
   }
}

void _fmpz_mpoly_submul_array1_fmpz_1(fmpz * poly1, 
                          const fmpz_t d, ulong exp2,
                           const fmpz * poly3, const ulong * exp3, slong len3)
{
   slong j;
   fmpz * c2, * c;

   c2 = poly1 + (slong) exp2;

   if (*d != 0)
   {
      for (j = 0; j < len3; j++)
      {
         c = c2 + (slong) exp3[j];
         fmpz_submul(c, d, poly3 + j);
      }
   }
}

void _fmpz_mpoly_to_ulong_array2(ulong * p, const fmpz * coeffs,
                                                 const ulong * exps, slong len)
{
   slong i, j;
   
   for (i = 0; i < len; i++)
   {
      ulong * ptr = p + 2*((slong) exps[i]);
      slong size = fmpz_size(coeffs + i);
      fmpz c = coeffs[i];

      if (!COEFF_IS_MPZ(c))
         ptr[0] = (ulong) c > 0 ? c : -c;
      else
      {
         __mpz_struct * m = COEFF_TO_PTR(c);

         for (j = 0; j < size; j++)
            ptr[j] = m->_mp_d[j];
      }
      
      if (fmpz_sgn(coeffs + i) < 0)
         mpn_neg(ptr, ptr, 2);
   }
}

void _fmpz_mpoly_to_ulong_array(ulong * p, const fmpz * coeffs,
                                                 const ulong * exps, slong len)
{
   slong i, j;
   
   for (i = 0; i < len; i++)
   {
      ulong * ptr = p + 3*((slong) exps[i]);
      slong size = fmpz_size(coeffs + i);
      fmpz c = coeffs[i];

      if (!COEFF_IS_MPZ(c))
         ptr[0] = (ulong) c > 0 ? c : -c;
      else
      {
         __mpz_struct * m = COEFF_TO_PTR(c);

         for (j = 0; j < size; j++)
            ptr[j] = m->_mp_d[j];
      }
      
      if (fmpz_sgn(coeffs + i) < 0)
         mpn_neg(ptr, ptr, 3);
   }
}

void _fmpz_mpoly_to_fmpz_array(fmpz * p, const fmpz * coeffs,
                                                 const ulong * exps, slong len)
{
   slong i;

   for (i = 0; i < len; i++)
      fmpz_set(p + (slong) exps[i], coeffs + i);
}

void _fmpz_mpoly_unpack_exponents_tight(ulong * e1, ulong * e2, slong len,
                               slong * mults, slong N, slong extra, slong bits)
{
   slong i, j;
   ulong exp;
   slong shift = FLINT_BITS - (N + extra)*bits;
   slong * prods;
   TMP_INIT;

   TMP_START;

   prods = (slong *) TMP_ALLOC((N + 1)*sizeof(slong));

   prods[0] = 1;
   for (i = 1; i <= N; i++)
      prods[i] = mults[i - 1]*prods[i - 1];

   for (i = 0; i < len; i++)
   {
      exp = 0;
         
      for (j = 0; j < N; j++)
         exp += (e2[i] % prods[j + 1])/prods[j] << bits*j;

      e1[i] = exp << shift;
   }

   TMP_END;
}

int _fmpz_mpoly_exponent_divides_tight(slong e1, slong e2,
                                                        slong * prods, slong N)
{
   slong j;

   for (j = 0; j < N; j++)
   {
      slong d1 = (e1 % prods[j + 1])/prods[j];
      slong d2 = (e2 % prods[j + 1])/prods[j];

      if (d1 < d2)
         return 0;
   }

   return 1;
}

slong _fmpz_mpoly_divexact_array_tight(fmpz ** poly1, ulong ** exp1,
                                           slong * alloc, slong len1, 
                      const fmpz * poly2, const ulong * exp2, slong len2,
                      const fmpz * poly3, const ulong * exp3, slong len3,
                                            slong * mults, slong bits, slong N)
{
   slong i, j, q, r, prod, bits1, bits2, bits3, len = len1;
   slong max3 = (slong) exp3[len3 - 1]; /* largest exponent in poly3 */
   slong min3 = (slong) exp3[0]; /* smallest exponent in poly3 */
   slong * prods;
   fmpz c3 = poly3[0];
   ulong u3 = ((ulong) FLINT_ABS(c3)) >> 1;
   fmpz * p1 = *poly1;
   ulong * e1 = *exp1;
   int small;
   TMP_INIT;

   TMP_START;

   prods = (slong *) TMP_ALLOC((N + 1)*sizeof(slong));

   prods[0] = 1;
   for (i = 1; i <= N; i++)
      prods[i] = mults[i - 1]*prods[i - 1];

   prod = prods[N];

   bits2 = _fmpz_vec_max_bits(poly2, len2);
   bits3 = _fmpz_vec_max_bits(poly3, len3);
   /* we assume a bound of FLINT_BITS - 2 for coefficients of the quotient */
   bits1 = FLINT_ABS(bits3) + FLINT_BITS + FLINT_BIT_COUNT(len3) - 2;

   small = FLINT_ABS(bits2) <= bits1 && FLINT_ABS(bits3) <= FLINT_BITS - 2;
   bits1 += 2; /* incr. so poly2 - q*poly3 doesn't overflow and for sign */

   if (small && bits1 <= 2*FLINT_BITS)
   {
      ulong * p2 = (ulong *) TMP_ALLOC(2*prod*sizeof(ulong));

      for (i = 0; i < 2*prod; i++)
         p2[i] = 0;

      _fmpz_mpoly_to_ulong_array2(p2, poly2, exp2, len2);

      for (i = 0; i + max3 < prod; i++)
      {
         ulong * ptr = p2 + 2*i;
         ulong p[2];

         if (ptr[0] != 0 || ptr[1] != 0)
         {
            if (0 > (slong) ptr[1])
               mpn_neg(p, ptr, 2);
            else
               flint_mpn_copyi(p, ptr, 2);

            if (u3 < p[1] || (u3 == 0 && 0 > (slong) p[0])) /* quotient too large */
            {
               for (j = len1; j < len; j++)
                  _fmpz_demote(p1 + j);
               len = len1;

               goto big;
            }

            sdiv_qrnnd(q, r, ptr[1], ptr[0], c3);

            if (COEFF_IS_MPZ(q))
            {
               for (j = len1; j < len; j++) /* quotient too large */
                  _fmpz_demote(p1 + j);
               len = len1;

               goto big;
            }

            if (r != 0 || /* not an exact division */
               !_fmpz_mpoly_exponent_divides_tight(i, min3, prods, N))
            {
               for (j = len1; j < len; j++)
                  _fmpz_demote(p1 + j);
               len = len1;

               goto cleanup;
            }

            _fmpz_mpoly_submul_array1_slong2_1(p2, q, i - min3,
                                                            poly3, exp3, len3);

            if (len >= *alloc)
            {
               p1 = (fmpz *) flint_realloc(p1, 2*(*alloc)*sizeof(fmpz));
               e1 = (ulong *) flint_realloc(e1, 2*(*alloc)*sizeof(ulong));
               flint_mpn_zero(p1 + *alloc, *alloc);
               (*alloc) *= 2;
            }
            
            p1[len] = q;
            e1[len++] = i - min3;         
         }
      }

      for ( ; i < prod; i++)
      {
         ulong * ptr = p2 + 2*i;

         if (ptr[0] != 0 || ptr[1] != 0)  /* not an exact division */
         {
            for (j = len1; j < len; j++)
               _fmpz_demote(p1 + j);
            len = len1;

            goto cleanup;
         }
      }
   }

   if (len == len1 && small)
   {
      ulong * p2 = (ulong *) TMP_ALLOC(3*prod*sizeof(ulong));

      for (i = 0; i < 3*prod; i++)
         p2[i] = 0;

      _fmpz_mpoly_to_ulong_array(p2, poly2, exp2, len2);

      for (i = 0; i + max3 < prod; i++)
      {
         ulong * ptr = p2 + 3*i;
         ulong p[3];

         if (ptr[0] != 0 || ptr[1] != 0 || ptr[2] != 0)
         {
            if (0 > (slong) ptr[2])
               mpn_neg(p, ptr, 3);
            else
               flint_mpn_copyi(p, ptr, 3);

            if (p[2] > 0 || u3 < p[1] || (u3 == 0 && 0 > (slong) p[0])) /* quotient too large */
            {
               for (j = len1; j < len; j++)
                  _fmpz_demote(p1 + j);
               len = len1;

               goto big;
            }

            sdiv_qrnnd(q, r, ptr[1], ptr[0], c3);

            if (COEFF_IS_MPZ(q))
            {
               for (j = len1; j < len; j++) /* quotient too large */
                  _fmpz_demote(p1 + j);
               len = len1;

               goto big;
            }

            if (r != 0 || /* not an exact division */
               !_fmpz_mpoly_exponent_divides_tight(i, min3, prods, N)) 
            {
               for (j = len1; j < len; j++)
                  _fmpz_demote(p1 + j);
               len = len1;

               goto cleanup;
            }

            _fmpz_mpoly_submul_array1_slong_1(p2, q, i - min3,
                                                            poly3, exp3, len3);

            if (len >= *alloc)
            {
               p1 = (fmpz *) flint_realloc(p1, 2*(*alloc)*sizeof(fmpz));
               e1 = (ulong *) flint_realloc(e1, 2*(*alloc)*sizeof(ulong));
               flint_mpn_zero(p1 + *alloc, *alloc);
               (*alloc) *= 2;
            }
            
            p1[len] = q;
            e1[len++] = i - min3;
         }
      }

      for ( ; i < prod; i++)
      {
         ulong * ptr = p2 + 3*i;

         if (ptr[0] != 0 || ptr[1] != 0 || ptr[2] != 0) 
         {
            /* not an exact division */

            for (j = len1; j < len; j++)
               _fmpz_demote(p1 + j);
            len = len1;

            goto cleanup;
         }
      }
   }

big:

   if (len == len1)
   {
      fmpz * p2 = (fmpz *) TMP_ALLOC(prod*sizeof(fmpz));
      fmpz_t fq, fr;

      fmpz_init(fq);
      fmpz_init(fr);

      for (i = 0; i < prod; i++)
         fmpz_init(p2 + i);

      _fmpz_mpoly_to_fmpz_array(p2, poly2, exp2, len2);
      
      for (i = 0; i + max3 < prod; i++)
      {
         if (!fmpz_is_zero(p2 + i))
         {
            fmpz_fdiv_qr(fq, fr, p2 + i, poly2 + 0);

            if (!fmpz_is_zero(fr) || /* not an exact division */
               !_fmpz_mpoly_exponent_divides_tight(i, min3, prods, N))
            {
               for (j = len1; j < len; j++)
                  _fmpz_demote(p1 + j);
               len = len1;

               goto cleanup2;
            }

            _fmpz_mpoly_submul_array1_fmpz_1(p2, fq, i - min3,
                                                            poly3, exp3, len3);

            if (len >= *alloc)
            {
               p1 = (fmpz *) flint_realloc(p1, 2*(*alloc)*sizeof(fmpz));
               e1 = (ulong *) flint_realloc(e1, 2*(*alloc)*sizeof(ulong));
               flint_mpn_zero(p1 + *alloc, *alloc);
               (*alloc) *= 2;
            }
            
            fmpz_set(p1 + len, fq);
            e1[len++] = i - min3;
         }
      }

      for ( ; i < prod; i++)
      {
         if (!fmpz_is_zero(p2 + i))
         {
            /* not an exact division */

            for (j = len1; j < len; j++)
               _fmpz_demote(p1 + j);
            len = len1;

            goto cleanup2;
         }
      }

cleanup2:

      fmpz_clear(fq);
      fmpz_clear(fr);

      for (i = 0; i < prod; i++)
         fmpz_clear(p2 + i);
   }
  
cleanup:

   (*poly1) = p1;
   (*exp1) = e1;

   TMP_END;

   return len - len1;
}

/*
   use array exact division to set poly1 to poly2/poly3 in N + 1 variables,
   given a list of multipliers to tightly pack exponents and a number of bits
   for the fields of the exponents of the result, assuming no aliasing
   classical exact division in main variable, array multiplication (submul)
   for multivariate coefficients
*/
slong _fmpz_mpoly_divexact_array_univariate(fmpz_mpoly_t poly1,
                    const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3,
                                            slong * mults, slong bits, slong N)
{
   slong i, j, k, prod, len = 0, len1, len2, len3, bits1, bits2, bits3, tlen, talloc;
   slong shift2 = FLINT_BITS - poly2->bits;
   slong shift3 = FLINT_BITS - poly3->bits;
   slong shift = FLINT_BITS - bits;
   slong * i1, * i2, * i3, * n1, * n2, * n3;
   ulong * exp2, * exp3, * texp;
   fmpz * temp;
   int small;
   TMP_INIT;

   prod = 1;
   for (i = 0; i < N; i++)
      prod *= mults[i];

   len2 = 1 + (slong) (poly2->exps[poly2->length - 1] >> shift2);
   len3 = 1 + (slong) (poly3->exps[poly3->length - 1] >> shift3);

   len1 = len2 - len3 + 1;

   TMP_START;

   /* compute indices and lengths of coefficients of polys in main variable */

   i1 = (slong *) TMP_ALLOC(len1*sizeof(slong));
   n1 = (slong *) TMP_ALLOC(len1*sizeof(slong));
   i2 = (slong *) TMP_ALLOC(len2*sizeof(slong));
   n2 = (slong *) TMP_ALLOC(len2*sizeof(slong));
   i3 = (slong *) TMP_ALLOC(len3*sizeof(slong));
   n3 = (slong *) TMP_ALLOC(len3*sizeof(slong));

   i2[0] = 0;
   j = 0;
   for (i = 0; i < len2 - 1; i++)
   {
      while (j < poly2->length && i == (slong) (poly2->exps[j] >> shift2))
         j++;

      i2[i + 1] = j;
      n2[i] = j - i2[i];
   }
   n2[len2 - 1] = poly2->length - j;

   i3[0] = 0;
   j = 0;
   for (i = 0; i < len3 - 1; i++)
   {
      while (j < poly3->length && i == (slong) (poly3->exps[j] >> shift3))
         j++;

      i3[i + 1] = j;
      n3[i] = j - i3[i];
   }
   n3[len3 - 1] = poly3->length - j;

   /* pack input coefficients tightly */

   exp2 = (ulong *) TMP_ALLOC(poly2->length*sizeof(ulong));
   exp3 = (ulong *) TMP_ALLOC(poly3->length*sizeof(ulong));

   _fmpz_mpoly_pack_exponents_tight(exp2, poly2->exps, poly2->length,
                                                        mults, N, 1, poly2->bits);

   _fmpz_mpoly_pack_exponents_tight(exp3, poly3->exps, poly3->length,
                                                        mults, N, 1, poly3->bits);

   bits2 = fmpz_mpoly_max_bits(poly2);
   bits3 = fmpz_mpoly_max_bits(poly3);
   /* we assume a bound of FLINT_BITS - 2 for coefficients of the quotient */
   bits1 = FLINT_ABS(bits3) + FLINT_BITS + FLINT_BIT_COUNT(poly3->length) - 2;

   small = FLINT_ABS(bits2) <= bits1 && FLINT_ABS(bits3) <= FLINT_BITS - 2;
   bits1 += 2; /* incr. so poly2 - q*poly3 doesn't overflow and for sign */

   /* make copy of first "coefficient" of poly2 */

   temp = (fmpz *) flint_calloc(n2[0] + 1, sizeof(fmpz));
   texp = (ulong *) flint_malloc((n2[0] + 1)*sizeof(ulong));
   talloc = n2[0] + 1;

   if (small && bits1 <= 2*FLINT_BITS)
   {
      ulong * p2 = (ulong *) TMP_ALLOC(2*prod*sizeof(ulong));

      for (i = 0; i < len2; i++)
      {
         for (j = 0; j < 2*prod; j++)
            p2[j] = 0;

         _fmpz_mpoly_to_ulong_array2(p2, poly2->coeffs + i2[i],
                                                          exp2 + i2[i], n2[i]);

         /* submuls */

         for (j = 0; j < i && j < len1; j++)
         {
            k = i - j;

            if (k < len3)
               _fmpz_mpoly_submul_array1_slong2(p2, poly1->coeffs + i1[j],
                          poly1->exps + i1[j], n1[j], poly3->coeffs + i3[k], 
                                                          exp3 + i3[k], n3[k]);
         }

         tlen = _fmpz_mpoly_from_ulong_array2(&temp, &texp, &talloc, 
                                                        p2, mults, N, bits, 0);

         if (i < len1)
         {

            _fmpz_mpoly_pack_exponents_tight(texp, texp, tlen, mults,
                                                                   N, 0, bits);

            i1[i] = len;
            
            if (tlen != 0)
            {
               n1[i] = _fmpz_mpoly_divexact_array_tight(&poly1->coeffs,
                        &poly1->exps, &poly1->alloc, len, temp, texp, tlen,
                                   poly3->coeffs, exp3, n3[0], mults, bits, N);

               if (n1[i] == 0) /* not an exact division */
               {
                  for (j = 0; j < len; j++)
                     _fmpz_demote(poly1->coeffs + j);
                  len = 0;

                  goto cleanup;
               }
            } else
               n1[i] = 0;

            /* check the quotient didn't have large coefficients */
            if (FLINT_ABS(_fmpz_vec_max_bits(poly1->coeffs + len, n1[i])) >
                                                                FLINT_BITS - 2)
            {
               for (j = 0; j < len; j++)
                  _fmpz_demote(poly1->coeffs + j);
               len = 0;

               goto big;
            }

            len += n1[i];
         } else /* check coefficient is zero */
         {
            for (j = 0; j < tlen; j++)
            {
               if (!fmpz_is_zero(temp + j)) /* not an exact division */
               {
                  for (j = 0; j < len; j++)
                     _fmpz_demote(poly1->coeffs + j);
                  len = 0;

                  goto cleanup;
               }
            }
         } 
      }
   }

   if (len == 0 && small)
   {
      ulong * p2 = (ulong *) TMP_ALLOC(3*prod*sizeof(ulong));

      for (i = 0; i < len2; i++)
      {
         for (j = 0; j < 3*prod; j++)
            p2[j] = 0;

         _fmpz_mpoly_to_ulong_array(p2, poly2->coeffs + i2[i],
                                                          exp2 + i2[i], n2[i]);

         /* submuls */

         for (j = 0; j < i && j < len1; j++)
         {
            k = i - j;

            if (k < len3)
               _fmpz_mpoly_submul_array1_slong(p2, poly1->coeffs + i1[j],
                          poly1->exps + i1[j], n1[j], poly3->coeffs + i3[k], 
                                                          exp3 + i3[k], n3[k]);
         }

         tlen = _fmpz_mpoly_from_ulong_array(&temp, &texp, &talloc, 
                                                        p2, mults, N, bits, 0);

         if (i < len1)
         {

            _fmpz_mpoly_pack_exponents_tight(texp, texp, tlen, mults,
                                                                   N, 0, bits);

            i1[i] = len;
            
            if (tlen != 0)
            {
               n1[i] = _fmpz_mpoly_divexact_array_tight(&poly1->coeffs,
                        &poly1->exps, &poly1->alloc, len, temp, texp, tlen,
                                   poly3->coeffs, exp3, n3[0], mults, bits, N);

               if (n1[i] == 0) /* not an exact division */
               {
                  for (j = 0; j < len; j++)
                     _fmpz_demote(poly1->coeffs + j);
                  len = 0;

                  goto cleanup;
               }
            } else
               n1[i] = 0;

            /* check the quotient didn't have large coefficients */
            if (FLINT_ABS(_fmpz_vec_max_bits(poly1->coeffs + len, n1[i])) >
                                                                FLINT_BITS - 2)
            {
               for (j = 0; j < len; j++)
                  _fmpz_demote(poly1->coeffs + j);
               len = 0;

               goto big;
            }

            len += n1[i];
         } else /* check coefficient is zero */
         {
            for (j = 0; j < tlen; j++)
            {
               if (!fmpz_is_zero(temp + j)) /* not an exact division */
               {
                  for (j = 0; j < len; j++)
                     _fmpz_demote(poly1->coeffs + j);
                  len = 0;

                  goto cleanup;
               }
            }
         } 
      }
   }

big:

   if (len == 0)
   {
      fmpz * p2 = (fmpz *) TMP_ALLOC(prod*sizeof(ulong));

      for (j = 0; j < prod; j++)
            fmpz_init(p2 + j);
      
      for (i = 0; i < len2; i++)
      {
         for (j = 0; j < prod; j++)
            fmpz_zero(p2 + j);

         _fmpz_mpoly_to_fmpz_array(p2, poly2->coeffs + i2[i],
                                                          exp2 + i2[i], n2[i]);

         /* submuls */

         for (j = 0; j < i && j < len1; j++)
         {
            k = i - j;

            if (k < len3)
               _fmpz_mpoly_submul_array1_fmpz(p2, poly1->coeffs + i1[j],
                          poly1->exps + i1[j], n1[j], poly3->coeffs + i3[k], 
                                                          exp3 + i3[k], n3[k]);
         }

         tlen = _fmpz_mpoly_from_fmpz_array(&temp, &texp, &talloc, 
                                                        p2, mults, N, bits, 0);

         if (i < len1)
         {

            _fmpz_mpoly_pack_exponents_tight(texp, texp, tlen, mults,
                                                                   N, 0, bits);

            i1[i] = len;
            
            if (tlen != 0)
            {
               n1[i] = _fmpz_mpoly_divexact_array_tight(&poly1->coeffs,
                        &poly1->exps, &poly1->alloc, len, temp, texp, tlen,
                                   poly3->coeffs, exp3, n3[0], mults, bits, N);

               if (n1[i] == 0) /* not an exact division */
               {
                  for (j = 0; j < len; j++)
                     _fmpz_demote(poly1->coeffs + j);
                  len = 0;

                  goto cleanup2;
               }
            } else
               n1[i] = 0;

            len += n1[i];
         } else /* check coefficient is zero */
         {
            for (j = 0; j < tlen; j++)
            {
               if (!fmpz_is_zero(temp + j)) /* not an exact division */
               {
                  for (j = 0; j < len; j++)
                     _fmpz_demote(poly1->coeffs + j);
                  len = 0;

                  goto cleanup2;
               }
            }
         } 
      }

cleanup2:

      for (j = 0; j < prod; j++)
            fmpz_clear(p2 + j);
   }

   if (len != 0)
   {
      _fmpz_mpoly_unpack_exponents_tight(poly1->exps, poly1->exps, len,
                                                            mults, N, 1, bits);

      /* put main variable back in */
      for (i = 0; i < len1; i++)
      {
         for (j = 0; j < n1[i]; j++)
            poly1->exps[i1[i] + j] += (i << shift);
      }
   }

cleanup:

   flint_free(temp);
   flint_free(texp);

   TMP_END;

   return len;
}

/*
   use array exact division to set poly1 to poly2/poly3 in N variables, given
   a list of multipliers to tightly pack exponents and a number of bits for the
   fields of the exponents of the result, assuming no aliasing
*/
slong _fmpz_mpoly_divexact_array(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                  const fmpz_mpoly_t poly3, slong * mults, slong bits, slong N)
{
   slong i;
   ulong * exp2, * exp3;
   slong len, prod;
   TMP_INIT;

   prod = 1;
   for (i = 0; i < N; i++)
      prod *= mults[i];

   if (prod > MAX_ARRAY_SIZE)
      return _fmpz_mpoly_divexact_array_univariate(poly1, poly2, poly3,
                                                               mults, bits, N - 1);

   TMP_START;

   exp2 = (ulong *) TMP_ALLOC(poly2->length*sizeof(ulong));
   exp3 = (ulong *) TMP_ALLOC(poly3->length*sizeof(ulong));

   _fmpz_mpoly_pack_exponents_tight(exp2, poly2->exps, poly2->length,
                                                     mults, N, 0, poly2->bits);

   _fmpz_mpoly_pack_exponents_tight(exp3, poly3->exps, poly3->length,
                                                     mults, N, 0, poly3->bits);

   len = _fmpz_mpoly_divexact_array_tight(&poly1->coeffs, &poly1->exps,
                     &poly1->alloc, 0,  poly2->coeffs, exp2, poly2->length,
                           poly3->coeffs, exp3, poly3->length, mults, bits, N);

   _fmpz_mpoly_unpack_exponents_tight(poly1->exps, poly1->exps, len,
                                                            mults, N, 0, bits);

   TMP_END;

   return len;
}

int fmpz_mpoly_divexact_array(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                          const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx)
{
   slong i, bits, exp_bits, N, len = 0, array_size;
   ulong * max_degs2;
   ulong max = 0;
   int res = 1;

   TMP_INIT;

   if (poly3->length == 0)
      flint_throw(FLINT_DIVZERO, "Divide by zero in fmpz_mpoly_divexact_monagan_pearce");

   if (poly2->length == 0)
   {
      fmpz_mpoly_zero(poly1, ctx);

      return 1;
   }

   TMP_START;

   max_degs2 = (ulong *) TMP_ALLOC(ctx->n*sizeof(ulong));

   fmpz_mpoly_max_degrees(max_degs2, poly2, ctx);

   for (i = 0; i < ctx->n; i++)
   {
      if (max_degs2[i] > max)
         max = max_degs2[i];
   }

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
      max_degs2[i]++;
      array_size *= max_degs2[i];
   }  
   max_degs2[ctx->n - 1]++;
   
   if (array_size > MAX_ARRAY_SIZE)
   {
      res = 0;

      goto cleanup;
   }

   if (poly2->exps[poly2->length - 1] < poly3->exps[poly3->length - 1])
      goto inexact;

   if (poly1 == poly2 || poly1 == poly3)
   {
      fmpz_mpoly_t temp;

      fmpz_mpoly_init2(temp, poly2->length/poly3->length + 1, ctx);
      fmpz_mpoly_fit_bits(temp, exp_bits, ctx);

      len = _fmpz_mpoly_divexact_array(temp, poly2, poly3,
                                        (slong *) max_degs2, exp_bits, ctx->n);

      fmpz_mpoly_swap(temp, poly1, ctx);

      fmpz_mpoly_clear(temp, ctx);
   } else
   {
      fmpz_mpoly_fit_length(poly1, poly2->length/poly3->length + 1, ctx);
      fmpz_mpoly_fit_bits(poly1, exp_bits, ctx);

      len = _fmpz_mpoly_divexact_array(poly1, poly2, poly3,
                                        (slong *) max_degs2, exp_bits, ctx->n);
   }

   _fmpz_mpoly_set_length(poly1, len, ctx);

inexact:

   if (len == 0)
   {
      TMP_END;

      flint_throw(FLINT_INEXACT, "Not an exact division in fmpz_mpoly_divexact_array");
   }

cleanup:

   TMP_END;

   return res;
}
