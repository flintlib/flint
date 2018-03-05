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

/* improve locality */
#define BLOCK 128
#define MAX_ARRAY_SIZE (WORD(300000))

/*
   Addmul into a dense array poly1, given polys with coefficients
   fitting into a word, and exponents tightly packed with mixed
   bases equal to the largest exponent for each variable, e.g.
   the input polys have exponents of the form
   a_0 + a_1*b1 + a_2*b_2*b_2 + .... where b_0, b_1, b_2, etc, are
   the bases, which are equal to the largest possible exponent for
   each of the respective variables in the exponent. These exponents
   are use as array indices in the output polynomial. The
   output poly is assumed to fit into one word per coefficient.
   The input polynomials are broken into blocks to improve
   cache efficiency.
*/
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

/*
   Addmul into a dense array poly1, given polys with coefficients
   fitting into a word, and exponents tightly packed with mixed
   bases equal to the largest exponent for each variable, e.g.
   the input polys have exponents of the form
   a_0 + a_1*b1 + a_2*b_2*b_2 + .... where b_0, b_1, b_2, etc, are
   the bases, which are equal to the largest possible exponent for
   each of the respective variables in the exponent. These exponents
   are use as array indices in the output polynomial. The
   output poly is assumed to fit into three words per coefficient.
   The input polynomials are broken into blocks to improve
   cache efficiency.
*/
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

/*
   Addmul into a dense array poly1, given polys with coefficients
   fitting into a word, and exponents tightly packed with mixed
   bases equal to the largest exponent for each variable, e.g.
   the input polys have exponents of the form
   a_0 + a_1*b1 + a_2*b_2*b_2 + .... where b_0, b_1, b_2, etc, are
   the bases, which are equal to the largest possible exponent for
   each of the respective variables in the exponent. These exponents
   are use as array indices in the output polynomial. The
   output poly is assumed to fit into two words per coefficient.
   The input polynomials are broken into blocks to improve
   cache efficiency.
*/
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

/*
   Addmul into a dense array poly1, given polys with coefficients
   fitting into a word, and exponents tightly packed with mixed
   bases equal to the largest exponent for each variable, e.g.
   the input polys have exponents of the form
   a_0 + a_1*b1 + a_2*b_2*b_2 + .... where b_0, b_1, b_2, etc, are
   the bases, which are equal to the largest possible exponent for
   each of the respective variables in the exponent. These exponents
   are use as array indices in the output polynomial. The
   output poly is unrestricted, having multiprecision coefficients.
   The input polynomials are broken into blocks to improve
   cache efficiency.
*/
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
                  fmpz_addmul(c, poly2 + i, poly3 + j);
               }
            }
         }
      }
   }
}

/* 
   Convert a polynomial in dense array format to an fmpz_mpoly with
   the given number of bits per exponent field. This function destroys
   poly2 and starts writing poly1 at index k. The function may reallocate
   its output. The array "mults" is a list of bases uses in encoding
   the array indices from the exponents. The value num is the number of
   fields in the output exponent vectors, also the number of entries in
   mults. Exponents are assumed to be packed into a single word. The
   array, poly2 is assumed to have three words per coefficient.
*/
slong _fmpz_mpoly_from_ulong_array(fmpz ** poly1, ulong ** exp1, slong * alloc, 
              ulong * poly2, const slong * mults, slong num, slong bits, slong k)
{
   slong i, j;
   ulong exp;
   ulong * c;
   slong * prods;
   fmpz * p1 = *poly1;
   ulong * e1 = *exp1;
   TMP_INIT;

   TMP_START;

   prods = (slong *) TMP_ALLOC((num + 1)*sizeof(slong));

   /*
      compute products 1, b0, b0*b1, b0*b1*b2 ...
      from list of bases b0, b1, b2, ...
   */
   prods[0] = 1;
   for (i = 1; i <= num; i++)
     prods[i] = mults[i - 1]*prods[i - 1];
   
   /* for each coeff in array */
   for (i = prods[num] - 1; i >= 0; i--)
   {
      c = poly2 + i*3;

      /* if coeff is nonzero */
      if (c[0] != 0 || c[1] != 0 || c[2] != 0)
      {
         _fmpz_mpoly_fit_length(&p1, &e1, alloc, k + 1, 1);

         exp = 0;
         
         /* compute exponent from index */
         for (j = 0; j < num; j++)
            exp += (i % prods[j + 1])/prods[j] << bits*j;

         /* shift exponent vector into place */
         e1[k] = exp;

         /* set coefficient */
         fmpz_set_signed_uiuiui(p1 + k, c[2], c[1], c[0]);
         
         k++;
      }
   }

   *poly1 = p1;
   *exp1 = e1;

   TMP_END;

   return k;
}

/*
   Convert a polynomial in dense array format to an fmpz_mpoly with
   the given number of bits per exponent field. This function destroys
   poly2 and starts writing poly1 at index k. The function may reallocate
   its output. The array "mults" is a list of bases uses in encoding
   the array indices from the exponents. The value num is the number of
   fields in the output exponent vectors, also the number of entries in
   mults. Exponents are assumed to be packed into a single word. The
   array, poly2 is assumed to have two words per coefficient.
*/
slong _fmpz_mpoly_from_ulong_array2(fmpz ** poly1, ulong ** exp1, slong * alloc, 
              ulong * poly2, const slong * mults, slong num, slong bits, slong k)
{
   slong i, j;
   ulong exp;
   ulong * c;
   slong * prods;
   fmpz * p1 = *poly1;
   ulong * e1 = *exp1;
   TMP_INIT;

   TMP_START;

   prods = (slong *) TMP_ALLOC((num + 1)*sizeof(slong));

   /*
      compute products 1, b0, b0*b1, b0*b1*b2 ...
      from list of bases b0, b1, b2, ...
   */
   prods[0] = 1;
   for (i = 1; i <= num; i++)
     prods[i] = mults[i - 1]*prods[i - 1];
   
   /* for each coeff in array */
   for (i = prods[num] - 1; i >= 0; i--)
   {
      c = poly2 + i*2;

      /* if coeff is nonzero */
      if (c[0] != 0 || c[1] != 0)
      {
         _fmpz_mpoly_fit_length(&p1, &e1, alloc, k + 1, 1);

         exp = 0;
         
         /* compute exponent from index */
         for (j = 0; j < num; j++)
            exp += (i % prods[j + 1])/prods[j] << bits*j;

         /* shift exponent vector into place */
         e1[k] = exp;

         /* set coefficient */
         fmpz_set_signed_uiui(p1 + k, c[1], c[0]);
         
         k++;
      }
   }

   *poly1 = p1;
   *exp1 = e1;

   TMP_END;

   return k;
}

/*
   Convert a polynomial in dense array format to an fmpz_mpoly with
   the given number of bits per exponent field. This function destroys
   poly2 and starts writing poly1 at index k. The function may reallocate
   its output. The array "mults" is a list of bases uses in encoding
   the array indices from the exponents. The value num is the number of
   fields in the output exponent vectors, also the number of entries in
   mults. Exponents are assumed to be packed into a single word. The
   array, poly2 is assumed to have three words per coefficient.
*/
slong _fmpz_mpoly_from_ulong_array1(fmpz ** poly1, ulong ** exp1, slong * alloc, 
              ulong * poly2, const slong * mults, slong num, slong bits, slong k)
{
   slong i, j;
   ulong exp;
   ulong * c;
   slong * prods;
   fmpz * p1 = *poly1;
   ulong * e1 = *exp1;
   TMP_INIT;

   TMP_START;

   prods = (slong *) TMP_ALLOC((num + 1)*sizeof(slong));

   /*
      compute products 1, b0, b0*b1, b0*b1*b2 ...
      from list of bases b0, b1, b2, ...
   */
   prods[0] = 1;
   for (i = 1; i <= num; i++)
     prods[i] = mults[i - 1]*prods[i - 1];

   /* for each coeff in array */
   for (i = prods[num] - 1; i >= 0; i--)
   {
      c = poly2 + i;

      /* if coeff is nonzero */
      if (c[0] != 0)
      {
         _fmpz_mpoly_fit_length(&p1, &e1, alloc, k + 1, 1);

         exp = 0;
         
         /* compute exponent from index */
         for (j = 0; j < num; j++)
            exp += (i % prods[j + 1])/prods[j] << bits*j;

         /* shift exponent vector into place */
         e1[k] = exp;
         
         /* set coefficient */
         fmpz_set_si(p1 + k, c[0]);
         
         k++;
      }
   }

   *poly1 = p1;
   *exp1 = e1;

   TMP_END;

   return k;
}

/*
   Convert a polynomial in dense array format to an fmpz_mpoly with
   the given number of bits per exponent field. This function destroys
   poly2 and starts writing poly1 at index k. The function may reallocate
   its output. The array "mults" is a list of bases uses in encoding
   the array indices from the exponents. The value num is the number of
   fields in the output exponent vectors, also the number of entries in
   mults. Exponents are assumed to be packed into a single word. The
   array, poly2 has no restrictions with respect to coefficients; they
   may be multiprecision integers.
*/
slong _fmpz_mpoly_from_fmpz_array(fmpz ** poly1, ulong ** exp1, slong * alloc, 
               fmpz * poly2, const slong * mults, slong num, slong bits, slong k)
{
   slong i, j;
   ulong exp;
   fmpz * c;
   slong * prods;
   fmpz * p1 = *poly1;
   ulong * e1 = *exp1;
   TMP_INIT;

   TMP_START;

   prods = (slong *) TMP_ALLOC((num + 1)*sizeof(slong));

   /*
      compute products 1, b0, b0*b1, b0*b1*b2 ...
      from list of bases b0, b1, b2, ...
   */
   prods[0] = 1;
   for (i = 1; i <= num; i++)
      prods[i] = mults[i - 1]*prods[i - 1];

   /* for each coeff in array */
   for (i = prods[num] - 1; i >= 0; i--)
   {
      c = poly2 + i;

      /* if coeff is nonzero */
      if (!fmpz_is_zero(c))
      {
         _fmpz_mpoly_fit_length(&p1, &e1, alloc, k + 1, 1);

         exp = 0;
         
         /* compute exponent from index */
         for (j = 0; j < num; j++)
            exp += (i % prods[j + 1])/prods[j] << bits*j;

         /* shift exponent vector into place */
         e1[k] = exp;
         
         /* set coefficient */
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
   Compute two bounds for the i-th multivariate "coefficient" (chunk) of
   the chunked polynomial. The first bound, stored at the i-th location
   in the  array b1 is the number of bits of the sum of the absolute value
   of the coefficients of the chunk. The second, stored in maxb1 is the
   maximum absolute value of the coefficients in the chunk. The array i1
   contains the starting indices of the chunks within the polynomial poly1,
   and the array n1 contains the corresponding lengths. Thus the i-th chunk
   starts at index i1[i] and has length n1[i].
*/
void _fmpz_mpoly_chunk_max_bits(slong * b1, slong * maxb1,
                           const fmpz * poly1, slong * i1, slong * n1, slong i)
{
   slong j;
   ulong hi = 0, lo = 0;

   maxb1[i] = 0;

   /* for each coeff in the chunk */
   for (j = 0; j < n1[i] && fmpz_fits_si(poly1 + i1[i] + j); j++)
   {
      slong c = fmpz_get_si(poly1 + i1[i] + j);
      ulong uc = (ulong) FLINT_ABS(c);

      /* compute max abs value of the coeff */
      if (FLINT_BIT_COUNT(uc) > maxb1[i])
         maxb1[i] = FLINT_BIT_COUNT(uc);

      /* sum the absolute values of the coeffs */
      add_ssaaaa(hi, lo, hi, lo, UWORD(0), uc);
   }

   if (j == n1[i]) /* no large coeffs encountered */
   {
      /* write out the number of bits */
      if (hi != 0)
         b1[i] = FLINT_BIT_COUNT(hi) + FLINT_BITS;
      else
         b1[i] = FLINT_BIT_COUNT(lo);
   } else /* restart with multiprecision routine */
   {
      fmpz_t sum;

      fmpz_init(sum);

      for (j = 0; j < n1[i]; j++)
      {
         slong size;

         if (fmpz_sgn(poly1 + i1[i] + j) < 0)
            fmpz_sub(sum, sum, poly1 + i1[i] + j);
         else
            fmpz_add(sum, sum, poly1 + i1[i] + j);

         size = fmpz_sizeinbase(poly1 + i1[i] + j, 2);

         if (size > maxb1[i])
            maxb1[i] = size;
      }

      b1[i] = fmpz_sizeinbase(sum, 2);

      fmpz_clear(sum);
   }
}

/*
   Use dense array multiplication to set poly1 to poly2*poly3 in num + 1
   variables, given a list of multipliers to tightly pack exponents and a
   number of bits for the fields of the exponents of the result, assuming
   no aliasing. Classical multiplication in main variable, array
   multiplication for multivariate coefficients in remaining num variables.
   The array "mults" is a list of bases to be used in encoding the array
   indices from the exponents. The function reallocates its output. 
*/
slong _fmpz_mpoly_mul_array_chunked(fmpz ** poly1, ulong ** exp1,
        slong * alloc, const fmpz * poly2, const ulong * exp2, slong len2,
                       const fmpz * poly3, const ulong * exp3, slong len3, 
                                          slong * mults, slong num, slong bits)
{
   slong i, j, k = 0, len, l1, l2, l3, prod, bits1, bits2 = 0, bits3 = 0;
   slong shift = bits*(num);
   slong * i2, * i3, * n2, * n3, * b2, * maxb2, * b3, * maxb3;
   ulong * e2, * e3, * p1;
   int small;
   TMP_INIT;

   /*
      compute products 1, b0, b0*b1, b0*b1*b2 ...
      from list of bases b0, b1, b2, ...
   */
   prod = 1;
   for (i = 0; i < num; i++)
      prod *= mults[i];

   /* compute lengths of poly2 and poly3 in chunks */
   l2 = 1 + (slong) (exp2[0] >> shift);
   l3 = 1 + (slong) (exp3[0] >> shift);

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
   
   /* compute chunks of the input polys with respect to the main variable */
   mpoly_main_variable_terms1(i2, n2, exp2, l2, len2, num + 1, num + 1, bits);
   mpoly_main_variable_terms1(i3, n3, exp3, l3, len3, num + 1, num + 1, bits);

   /* pack input exponents tightly with mixed bases specified by "mults" */
   e2 = (ulong *) TMP_ALLOC(len2*sizeof(ulong));
   e3 = (ulong *) TMP_ALLOC(len3*sizeof(ulong));

   mpoly_pack_monomials_tight(e2, exp2, len2, mults, num, bits);
   mpoly_pack_monomials_tight(e3, exp3, len3, mults, num, bits);

   /* work out max bits for each chunk and optimal bits */

   for (i = 0; i < l2; i++)
   {
      _fmpz_mpoly_chunk_max_bits(b2, maxb2, poly2, i2, n2, i);

      if (bits2 < maxb2[i])
         bits2 = maxb2[i];
   }

   for (i = 0; i < l3; i++)
   {
      _fmpz_mpoly_chunk_max_bits(b3, maxb3, poly3, i3, n3, i);

      if (bits3 < maxb3[i])
         bits3 = maxb3[i];
   }

   /* whether the output coefficients are "small" */
   small = bits2 <= (FLINT_BITS - 2) &&
           bits3 <= (FLINT_BITS - 2);

   /* classical multiplication one output chunk at a time */

   l1 = l2 + l3 - 1; /* length of output in chunks */

   if (small)
   {
      p1 = (ulong *) TMP_ALLOC(3*prod*sizeof(ulong));

      /* for each output chunk */
      for (i = 0; i < l1; i++)
      {
         slong num1 = 0;
         bits1 = 0;

         /* compute bound on coeffs of output chunk */
         
         for (j = 0; j < l2 && j <= i; j++)
         {
            /* for each cross product of chunks */
            if (i - j < l3)
            {
               bits1 = FLINT_MAX(bits1, FLINT_MIN(b2[j] +
                           maxb3[i - j], maxb2[j] + b3[i - j]));
               num1++;
            }
         }

         bits1 += FLINT_BIT_COUNT(num1) + 1; /* includes one bit for sign */

         /* output coeffs fit in one word */
         if (bits1 <= FLINT_BITS)
         {
            for (j = 0; j < prod; j++)
               p1[j] = 0;

            /* addmuls for each cross product of chunks */
            for (j = 0; j < l2 && j <= i; j++)
            {
               if (i - j < l3)
               {
                  _fmpz_mpoly_addmul_array1_slong1(p1, 
                     (slong *) poly2 + i2[j], e2 + i2[j], n2[j],
                     (slong *) poly3 + i3[i - j], e3 + i3[i - j], n3[i - j]);
               }
            }

            /* convert array to fmpz_poly */
            len = _fmpz_mpoly_from_ulong_array1(poly1, exp1, alloc, 
                                                      p1, mults, num, bits, k) - k;

            /* insert main variable into exponents */
            for (j = 0; j < len; j++)
               (*exp1)[k + j] += ((l1 - i - 1) << shift);

            k += len;
         } else if (bits1 <= 2*FLINT_BITS) /* output coeffs fit in two words */        
         {
            for (j = 0; j < 2*prod; j++)
               p1[j] = 0;

            /* addmuls for each cross product of chunks */
            for (j = 0; j < l2 && j <= i; j++)
            {
               if (i - j < l3)
               {
                  _fmpz_mpoly_addmul_array1_slong2(p1, 
                     (slong *) poly2 + i2[j], e2 + i2[j], n2[j],
                     (slong *) poly3 + i3[i - j], e3 + i3[i - j], n3[i - j]);

               }
            }

            /* convert array to fmpz_poly */
            len = _fmpz_mpoly_from_ulong_array2(poly1, exp1, alloc, 
                                                      p1, mults, num, bits, k) - k;

            /* insert main variable into exponents */
            for (j = 0; j < len; j++)
               (*exp1)[k + j] += ((l1 - i - 1) << shift);

            k += len;
         } else /* output coeffs fit in three words */
         {
            for (j = 0; j < 3*prod; j++)
               p1[j] = 0;

            /* addmuls for each cross product of chunks */
            for (j = 0; j < l2 && j <= i; j++)
            {
               if (i - j < l3)
               {
                  _fmpz_mpoly_addmul_array1_slong(p1, 
                     (slong *) poly2 + i2[j], e2 + i2[j], n2[j],
                     (slong *) poly3 + i3[i - j], e3 + i3[i - j], n3[i - j]);

               }
            }

            /* convert array to fmpz_poly */
            len = _fmpz_mpoly_from_ulong_array(poly1, exp1, alloc, 
                                                      p1, mults, num, bits, k) - k;

            /* insert main variable into exponents */
            for (j = 0; j < len; j++)
               (*exp1)[k + j] += ((l1 - i - 1) << shift);

            k += len;
         }
      }
   } else /* output coeffs may be arbitrary size */
   {
      fmpz * p1 = (fmpz *) TMP_ALLOC(prod*sizeof(fmpz));
     
      /* for each output chunk */
      for (i = 0; i < l1; i++)
      {
         for (j = 0; j < prod; j++)
            p1[j] = 0;

         /* addmuls for each cross product of chunks */
         for (j = 0; j < l2 && j <= i; j++)
         {
            if (i - j < l3)
            {
               _fmpz_mpoly_addmul_array1_fmpz(p1, 
                   poly2 + i2[j], e2 + i2[j], n2[j],
                   poly3 + i3[i - j], e3 + i3[i - j], n3[i - j]);

             }
          }

          /* convert array to fmpz_poly */
          len = _fmpz_mpoly_from_fmpz_array(poly1, exp1, alloc, 
                                                      p1, mults, num, bits, k) - k;

          /* insert main variable into exponents */
          for (j = 0; j < len; j++)
             (*exp1)[k + j] += ((l1 - i - 1) << shift);

          for (j = 0; j < prod; j++)
             _fmpz_demote(p1 + j);

          k += len;
       }
   }

   TMP_END;

   return k;
}

/*
   Use array multiplication to set poly1 to poly2*poly3 in num variables, given
   a list of multipliers to tightly pack exponents and a number of bits for the
   fields of the exponents of the result, assuming no aliasing. The array "mults"
   is a list of bases to be used in encoding the array indices from the exponents.
   The function reallocates its output.
*/
slong _fmpz_mpoly_mul_array(fmpz ** poly1, ulong ** exp1, slong * alloc,
                         const fmpz * poly2, const ulong * exp2, slong len2, 
                         const fmpz * poly3, const ulong * exp3, slong len3, 
                                          slong * mults, slong num, slong bits)
{
   slong i, bits1, bits2, bits3;
   ulong * e2, * e3;
   slong prod, len;
   int small;
   ulong hi = 0, lo = 0; /* two words to count bits */
   TMP_INIT;

   /*
      compute products 1, b0, b0*b1, b0*b1*b2 ...
      from list of bases b0, b1, b2, ...
   */
   prod = 1;
   for (i = 0; i < num; i++)
      prod *= mults[i];

   /* if array size will be too large, chunk the polynomials */
   if (prod > MAX_ARRAY_SIZE)
      return _fmpz_mpoly_mul_array_chunked(poly1, exp1, alloc,
                   poly2, exp2, len2, poly3, exp3, len3, mults, num - 1, bits);

   TMP_START;

   /* pack input exponents tightly with mixed bases specified by "mults" */
   
   e2 = (ulong *) TMP_ALLOC(len2*sizeof(ulong));
   e3 = (ulong *) TMP_ALLOC(len3*sizeof(ulong));

   mpoly_pack_monomials_tight(e2, exp2, len2, mults, num, bits);
   mpoly_pack_monomials_tight(e3, exp3, len3, mults, num, bits);

   /* compute bound on output bits and whether they are "small" */
   bits2 = _fmpz_vec_max_bits(poly2, len2);
   bits3 = _fmpz_vec_max_bits(poly3, len3);

   small = FLINT_ABS(bits2) <= (FLINT_BITS - 2) &&
           FLINT_ABS(bits3) <= (FLINT_BITS - 2);

   bits1 = -1; /* not used in large case */

   /* if output coeffs are "small" */
   if (small)
   {
      /* compute bound : sum of absolute value of coeffs of shorter poly */
      if (len2 < len3)
      {
         for (i = 0; i < len2; i++)
         {
            slong b2 = fmpz_get_si(poly2 + i);

            add_ssaaaa(hi, lo, hi, lo, UWORD(0), (ulong) FLINT_ABS(b2));
         }

         bits1 = FLINT_ABS(bits3) + 1; /* one bit for sign */
      } else
      {
         for (i = 0; i < len3; i++)
         {
            slong b3 = fmpz_get_si(poly3 + i);

            add_ssaaaa(hi, lo, hi, lo, UWORD(0), (ulong) FLINT_ABS(b3));
         }

         bits1 = FLINT_ABS(bits2) + 1; /* one bit for sign */
      }

      /* compute number of bits of sum of absolute values */
      if (hi != 0)
         bits1 += FLINT_BIT_COUNT(hi) + FLINT_BITS;
      else
         bits1 += FLINT_BIT_COUNT(lo);
   } 

   /* output coeffs fit in one word */
   if (small && bits1 <= FLINT_BITS)
   {
      ulong * p1 = (ulong *) TMP_ALLOC(prod*sizeof(ulong));

      for (i = 0; i < prod; i++)
         p1[i] = 0;

      /* array multiplication */
      _fmpz_mpoly_addmul_array1_slong1(p1, 
                         (slong *) poly2, e2, len2, (slong *) poly3, e3, len3);

      /* convert to fmpz_poly */
      len = _fmpz_mpoly_from_ulong_array1(poly1, exp1, alloc, 
                                                      p1, mults, num, bits, 0);

   } else if (small && bits1 <= 2*FLINT_BITS) /* output coeffs in two words */
   {
      ulong * p1 = (ulong *) TMP_ALLOC(2*prod*sizeof(ulong));

      for (i = 0; i < 2*prod; i++)
         p1[i] = 0;

      /* array multiplication */
      _fmpz_mpoly_addmul_array1_slong2(p1, 
                         (slong *) poly2, e2, len2, (slong *) poly3, e3, len3);

      /* convert to fmpz_poly */
      len = _fmpz_mpoly_from_ulong_array2(poly1, exp1, alloc, 
                                                      p1, mults, num, bits, 0);
   } else if (small) /* three words per output coeff */
   {
      ulong * p1 = (ulong *) TMP_ALLOC(3*prod*sizeof(ulong));

      for (i = 0; i < 3*prod; i++)
         p1[i] = 0;

      /* array multiplication */
      _fmpz_mpoly_addmul_array1_slong(p1, 
                         (slong *) poly2, e2, len2, (slong *) poly3, e3, len3);

      /* convert to fmpz_poly */
      len = _fmpz_mpoly_from_ulong_array(poly1, exp1, alloc, 
                                                      p1, mults, num, bits, 0);      
   } else /* multiprecision output coeffs */
   {
      fmpz * p1 = (fmpz *) TMP_ALLOC(prod*sizeof(fmpz));

      for (i = 0; i < prod; i++)
         p1[i] = 0;

      /* array multiplication */
      _fmpz_mpoly_addmul_array1_fmpz(p1, poly2, e2, len2, poly3, e3, len3);

      /* convert to fmpz_poly */
      len = _fmpz_mpoly_from_fmpz_array(poly1, exp1, alloc, 
                                                      p1, mults, num, bits, 0);

      for (i = 0; i < prod; i++)
         _fmpz_demote(p1 + i);
   }
  
   TMP_END;

   return len;
}

int fmpz_mpoly_mul_array(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                          const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx)
{
   slong i, bits, exp_bits, N, len = 0, array_size;
   ulong max, max2, max3, * max_fields2, * max_fields3;
   ulong * exp2 = poly2->exps, * exp3 = poly3->exps;
   int free2 = 0, free3 = 0;
   int res = 1;

   TMP_INIT;

   /* input poly is zero */
   if (poly2->length == 0 || poly3->length == 0)
   {
      fmpz_mpoly_zero(poly1, ctx);

      return 1;
   }

   TMP_START;

    /* compute maximum exponents for each variable */
    max_fields2 = (ulong *) TMP_ALLOC(ctx->minfo->nfields*sizeof(ulong));
    max_fields3 = (ulong *) TMP_ALLOC(ctx->minfo->nfields*sizeof(ulong));
    mpoly_max_fields_ui(max_fields2, poly2->exps, poly2->length,
                                                      poly2->bits, ctx->minfo);
    mpoly_max_fields_ui(max_fields3, poly3->exps, poly3->length,
                                                      poly3->bits, ctx->minfo);
    max2 = max3 = 0;
    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        if (max_fields2[i] > max2)
            max2 = max_fields2[i];

        if (max_fields3[i] > max3)
            max3 = max_fields3[i];
    }

   /* check that exponents won't overflow a word */
   max = max2 + max3;
   if (max < max2 || 0 > (slong) max)
      flint_throw(FLINT_EXPOF, "Exponent overflow in fmpz_mpoly_mul_array");

    /* compute number of bits required for output exponents */
    bits = FLINT_BIT_COUNT(max);
    exp_bits = FLINT_MAX(WORD(8), bits + 1);
    exp_bits = FLINT_MAX(exp_bits, poly2->bits);
    exp_bits = FLINT_MAX(exp_bits, poly3->bits);

   N = mpoly_words_per_exp(exp_bits, ctx->minfo);

   /* array multiplication expects each exponent vector in one word */
   /* current code is wrong for reversed orderings */
   if (N != 1 || mpoly_ordering_isrev(ctx->minfo))
   {
      res = 0;
      goto cleanup;
   }

   /* compute bounds on output exps, used as mixed bases for packing exps */
   array_size = 1;
   for (i = 0; i < ctx->minfo->nfields - 1; i++)
   {
      max_fields2[i] += max_fields3[i] + 1;
      array_size *= max_fields2[i];
   }
   max_fields2[ctx->minfo->nfields - 1] += max_fields3[ctx->minfo->nfields - 1] + 1;

   /* if exponents too large for array multiplication, exit silently */
   if (array_size > MAX_ARRAY_SIZE)
   {
      res = 0;
      goto cleanup;
   }

   /* expand input exponents to same number of bits as output */
   if (exp_bits > poly2->bits)
   {
      free2 = 1;
      exp2 = (ulong *) flint_malloc(N*poly2->length*sizeof(ulong));
      mpoly_repack_monomials(exp2, exp_bits, poly2->exps, poly2->bits,
                                                    poly2->length, ctx->minfo);
   }

   if (exp_bits > poly3->bits)
   {
      free3 = 1;
      exp3 = (ulong *) flint_malloc(N*poly3->length*sizeof(ulong));
      mpoly_repack_monomials(exp3, exp_bits, poly3->exps, poly3->bits,
                                                    poly3->length, ctx->minfo);
   }

   /* handle aliasing and do array multiplication */
   if (poly1 == poly2 || poly1 == poly3)
   {
      fmpz_mpoly_t temp;

      fmpz_mpoly_init2(temp, poly2->length + poly3->length - 1, ctx);
      fmpz_mpoly_fit_bits(temp, exp_bits, ctx);
      temp->bits = exp_bits;

      if (poly2->length >= poly3->length)
         len = _fmpz_mpoly_mul_array(&temp->coeffs, &temp->exps, &temp->alloc, 
                                           poly3->coeffs, exp3, poly3->length,
                                           poly2->coeffs, exp2, poly2->length,
                         (slong *) max_fields2, ctx->minfo->nfields, exp_bits);
      else
         len = _fmpz_mpoly_mul_array(&temp->coeffs, &temp->exps, &temp->alloc, 
                                           poly2->coeffs, exp2, poly2->length,
                                           poly3->coeffs, exp3, poly3->length,
                         (slong *) max_fields2, ctx->minfo->nfields, exp_bits);

      fmpz_mpoly_swap(temp, poly1, ctx);

      fmpz_mpoly_clear(temp, ctx);
   } else
   {
      fmpz_mpoly_fit_length(poly1, poly2->length + poly3->length - 1, ctx);
      fmpz_mpoly_fit_bits(poly1, exp_bits, ctx);
      poly1->bits = exp_bits;

      if (poly2->length >= poly3->length)
         len = _fmpz_mpoly_mul_array(&poly1->coeffs, &poly1->exps, &poly1->alloc,
                                            poly3->coeffs, exp3, poly3->length,
                                            poly2->coeffs, exp2, poly2->length,
                         (slong *) max_fields2, ctx->minfo->nfields, exp_bits);
      else
         len = _fmpz_mpoly_mul_array(&poly1->coeffs, &poly1->exps, &poly1->alloc,
                                            poly2->coeffs, exp2, poly2->length, 
                                            poly3->coeffs, exp3, poly3->length,
                         (slong *) max_fields2, ctx->minfo->nfields, exp_bits);
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
