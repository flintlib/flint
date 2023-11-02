/*
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mpoly.h"

/* improve locality */
#define BLOCK 128
#define MAX_ARRAY_SIZE (WORD(300000))

/*
   Submul into a dense array poly1, given polys with coefficients
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
                  c[2] += (0 <= (slong) p[1]) ? cy : cy + 1;
               }
            }
         }
      }
   }
}

/*
   Submul into a dense array poly1, given polys with coefficients
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

/*
   Submul into a dense array poly1, given polys with coefficients
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
void _fmpz_mpoly_submul_array1_slong1(ulong * poly1,
                 const slong * poly2, const ulong * exp2, slong len2,
                           const slong * poly3, const ulong * exp3, slong len3)
{
   slong ii, i, jj, j;
   ulong p; /* for products of coefficients */
   ulong * c2;

   for (ii = 0; ii < len2 + BLOCK; ii += BLOCK)
   {
      for (jj = 0; jj < len3 + BLOCK; jj += BLOCK)
      {
         for (i = ii; i < FLINT_MIN(ii + BLOCK, len2); i++)
         {
            c2 = poly1 + ((slong) exp2[i]);

            if (poly2[i] != 0)
            {
               for (j = jj; j < FLINT_MIN(jj + BLOCK, len3); j++)
               {
                  p = (ulong) ((slong) poly2[i])*((slong) poly3[j]);
                  c2[(slong) exp3[j]] -= p;
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
                  fmpz_submul(c, poly2 + i, poly3 + j);
               }
            }
         }
      }
   }
}

/*
   Polynomial by monomial submul into a dense array poly1, given
   an input poly and monomial with coefficients fitting into a word,
   and exponents tightly packed with mixed bases equal to the
   largest exponent for each variable, e.g. the input poly and monomial
   have exponents of the form a_0 + a_1*b1 + a_2*b_2*b_2 + ....
   where b_0, b_1, b_2, etc, are the bases, which are equal to the
   largest possible exponent for each of the respective variables in
   the exponent. These exponents are use as array indices in the output
   polynomial. The output poly is assumed to fit into three words per
   coefficient.
*/
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
         c[2] += (0 <= (slong) p[1]) ? cy : cy + 1;
      }
   }
}

/*
   Polynomial by monomial submul into a dense array poly1, given
   an input poly and monomial with coefficients fitting into a word,
   and exponents tightly packed with mixed bases equal to the
   largest exponent for each variable, e.g. the input poly and monomial
   have exponents of the form a_0 + a_1*b1 + a_2*b_2*b_2 + ....
   where b_0, b_1, b_2, etc, are the bases, which are equal to the
   largest possible exponent for each of the respective variables in
   the exponent. These exponents are use as array indices in the output
   polynomial. The output poly is assumed to fit into two words per
   coefficient.
*/
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

/*
   Polynomial by monomial submul into a dense array poly1, given
   an input poly and monomial with coefficients fitting into a word,
   and exponents tightly packed with mixed bases equal to the
   largest exponent for each variable, e.g. the input poly and monomial
   have exponents of the form a_0 + a_1*b1 + a_2*b_2*b_2 + ....
   where b_0, b_1, b_2, etc, are the bases, which are equal to the
   largest possible exponent for each of the respective variables in
   the exponent. These exponents are use as array indices in the output
   polynomial. The output poly is unrestricted, having multiprecision
   coefficients.
*/
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

/*
   Convert an fmpz_mpoly to dense array format, where the exponents
   of the input poly are tightly packed with mixed bases equal to the
   largest exponent for each variable, e.g. the input polys have
   exponents of the form a_0 + a_1*b1 + a_2*b_2*b_2 + .... where
   b_0, b_1, b_2, etc, are the bases, which are equal to the largest
   possible exponent for each of the respective variables in the
   exponent. The output array is assumed to have two words per
   coefficient.
*/
void _fmpz_mpoly_to_ulong_array2(ulong * p, const fmpz * coeffs,
                                                 const ulong * exps, slong len)
{
   slong i, j;

   /* for each term of the input poly */
   for (i = 0; i < len; i++)
   {
      ulong * ptr = p + 2*((slong) exps[i]);
      slong size = fmpz_size(coeffs + i);
      fmpz c = coeffs[i];

      if (!COEFF_IS_MPZ(c))
      {
         ptr[0] = (ulong) c;
         ptr[1] = (ulong) (c > 0 ? 0 : -WORD(1));
      }
      else
      {
         __mpz_struct * m = COEFF_TO_PTR(c);

         for (j = 0; j < size; j++)
            ptr[j] = m->_mp_d[j];

         if (fmpz_sgn(coeffs + i) < 0)
            mpn_neg(ptr, ptr, 2);
      }
   }
}

/*
   Convert an fmpz_mpoly to dense array format, where the exponents
   of the input poly are tightly packed with mixed bases equal to the
   largest exponent for each variable, e.g. the input polys have
   exponents of the form a_0 + a_1*b1 + a_2*b_2*b_2 + .... where
   b_0, b_1, b_2, etc, are the bases, which are equal to the largest
   possible exponent for each of the respective variables in the
   exponent. The output array is assumed to have one word per
   coefficient.
*/
void _fmpz_mpoly_to_ulong_array1(ulong * p, const fmpz * coeffs,
                                                 const ulong * exps, slong len)
{
   slong i;

   for (i = 0; i < len; i++)
   {
      ulong * ptr = p + ((slong) exps[i]);
      slong size = fmpz_size(coeffs + i);
      fmpz c = coeffs[i];

      if (!COEFF_IS_MPZ(c))
         ptr[0] = c;
      else
      {
         __mpz_struct * m = COEFF_TO_PTR(c);

         if (size != 0)
         {
            if (fmpz_sgn(coeffs + i) > 0)
               ptr[0] = m->_mp_d[0];
            else
               ptr[0] = -m->_mp_d[0];
         }
      }
   }
}

/*
   Convert an fmpz_mpoly to dense array format, where the exponents
   of the input poly are tightly packed with mixed bases equal to the
   largest exponent for each variable, e.g. the input polys have
   exponents of the form a_0 + a_1*b1 + a_2*b_2*b_2 + .... where
   b_0, b_1, b_2, etc, are the bases, which are equal to the largest
   possible exponent for each of the respective variables in the
   exponent. The output array is assumed to have three words per
   coefficient.
*/
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
      {
         ptr[0] = (ulong) c;
         if (c > 0)
         {
            ptr[1] = 0;
            ptr[2] = 0;
         } else
         {
            ptr[1] = (ulong) -WORD(1);
            ptr[2] = (ulong) -WORD(1);
         }
      }
      else
      {
         __mpz_struct * m = COEFF_TO_PTR(c);

         for (j = 0; j < size; j++)
            ptr[j] = m->_mp_d[j];

         if (fmpz_sgn(coeffs + i) < 0)
            mpn_neg(ptr, ptr, 3);
      }
   }
}

/*
   Convert an fmpz_mpoly to dense array format, where the exponents
   of the input poly are tightly packed with mixed bases equal to the
   largest exponent for each variable, e.g. the input polys have
   exponents of the form a_0 + a_1*b1 + a_2*b_2*b_2 + .... where
   b_0, b_1, b_2, etc, are the bases, which are equal to the largest
   possible exponent for each of the respective variables in the
   exponent. The output array is unrestricted, having multiprecision
   coefficients.
*/
void _fmpz_mpoly_to_fmpz_array(fmpz * p, const fmpz * coeffs,
                                                 const ulong * exps, slong len)
{
   slong i;

   for (i = 0; i < len; i++)
      fmpz_set(p + (slong) exps[i], coeffs + i);
}

/*
   Set poly1 to the quotient of poly2 by poly3 if the quotient is exact,
   and return the length of the quotient. If the quotient is not exact,
   return 0. The algorithm aborts as early as possible if the quotient is
   not exact. It is assumed that poly2 is nonzero so that the quotient is
   either exact and nonzero, or inexact. The polynomials have their
   exponents tightly packed,  with mixed bases equal to the largest
   exponent for each variable, e.g. the input polys have exponents of the
   form a_0 + a_1*b1 + a_2*b_2*b_2 + .... where b_0, b_1, b_2, etc, are
   the bases, which are equal to the largest possible exponent for
   each of the respective variables in the exponent. There are assumed to
   be "num" variables and the bases b_i are passed in the array "mults".
   The function reallocates its output. The quotient poly is not assumed
   to be zero on input. The quotient are appended to the existing terms in
   that polys.
*/
slong _fmpz_mpoly_divides_array_tight(fmpz ** poly1, ulong ** exp1,
                                           slong * alloc, slong len1,
                      const fmpz * poly2, const ulong * exp2, slong len2,
                      const fmpz * poly3, const ulong * exp3, slong len3,
                                                      slong * mults, slong num)
{
   slong i, j, q, r, prod, bits1, bits2, bits3, len = len1;
   slong max3 = (slong) exp3[0]; /* largest exponent in poly3 */
   slong min3 = (slong) exp3[len3 - 1]; /* smallest exponent in poly3 */
   slong * prods;
   fmpz c3 = poly3[0];
   /* abs val of trailing coeff of poly3 */
   ulong u3 = ((ulong) FLINT_ABS(c3)) >> 1;
   fmpz * p1 = *poly1;
   ulong * e1 = *exp1;
   int small;
   TMP_INIT;

   TMP_START;

   /* check there are at least as many zero exponents in dividend as divisor */
   if (exp2[len2 - 1] < min3)
      goto cleanup; /* not an exact quotient */

   prods = (slong *) TMP_ALLOC((num + 1)*sizeof(slong));

   /*
      compute products 1, b0, b0*b1, b0*b1*b2 ...
      from list of bases b0, b1, b2, ...
   */
   prods[0] = 1;
   for (i = 1; i <= num; i++)
      prods[i] = mults[i - 1]*prods[i - 1];

   prod = prods[num];

   /* quick check leading terms divide */
   if (!mpoly_monomial_divides_tight(exp2[0], exp3[0], prods, num))
      goto cleanup;

   /* compute bound on poly2 - q*poly3 assuming quotient remains small */
   bits2 = _fmpz_vec_max_bits(poly2, len2);
   bits3 = _fmpz_vec_max_bits(poly3, len3);
   /* we assume a bound of SMALL_FMPZ_BITCOUNT_MAX for coefficients of the quotient */
   bits1 = FLINT_ABS(bits3) + FLINT_BITS + FLINT_BIT_COUNT(len3) - 2;

   small = FLINT_ABS(bits2) <= bits1 && FLINT_ABS(bits3) <= SMALL_FMPZ_BITCOUNT_MAX;
   bits1 += 2; /* incr. so poly2 - q*poly3 doesn't overflow and for sign */

   /* input coeffs small and intermediate computations fit two words */
   if (small && bits1 <= 2*FLINT_BITS)
   {
      ulong * p2 = (ulong *) TMP_ALLOC(2*prod*sizeof(ulong));

      for (i = 0; i < 2*prod; i++)
         p2[i] = 0;

      /* poly2 to array format */
      _fmpz_mpoly_to_ulong_array2(p2, poly2, exp2, len2);

      /* for each term of poly2 array relevant to exact quotient */
      for (i = prod - 1; i >= max3; i--)
      {
         ulong * ptr = p2 + 2*i;
         ulong p[2];

         /* if coeff is nonzero */
         if (ptr[0] != 0 || ptr[1] != 0)
         {
            if (0 > (slong) ptr[1])
               mpn_neg(p, ptr, 2);
            else
               flint_mpn_copyi(p, ptr, 2);

            /* check quotient won't overflow a word */
            if (u3 <= p[1] || (u3 == 0 && 0 > (slong) p[0])) /* quotient too large */
            {
               for (j = len1; j < len; j++)
                  _fmpz_demote(p1 + j);
               len = len1;

               goto big;
            }

            /* quotient and remainder of coeffs */
            sdiv_qrnnd(q, r, ptr[1], ptr[0], c3);

            /* check coefficient is small, else restart with multiprec code */
            if (COEFF_IS_MPZ(FLINT_ABS(q)))
            {
               for (j = len1; j < len; j++) /* quotient too large */
                  _fmpz_demote(p1 + j);
               len = len1;

               goto big;
            }

            /* check coeff and monomial quotients were exact */
            if (r != 0 || /* not an exact division */
               !mpoly_monomial_divides_tight(i, max3, prods, num))
            {
               for (j = len1; j < len; j++)
                  _fmpz_demote(p1 + j);
               len = len1;

               goto cleanup;
            }

            /* submul q*poly3 */
            _fmpz_mpoly_submul_array1_slong2_1(p2, q, i - max3,
                                                            poly3, exp3, len3);

            /* reallocate quotient */
            _fmpz_mpoly_fit_length(&p1, &e1, alloc, len + 1, 1);

            /* write quotient term */
            fmpz_set_si(p1 + len, q);
            e1[len++] = i - max3;
         }
      }

      /* check there are no nonzero terms left in array */
      for ( ; i >= min3; i--)
      {
         ulong * ptr = p2 + 2*i;

         /* if coeff nonzero */
         if (ptr[0] != 0 || ptr[1] != 0)  /* not an exact division */
         {
            for (j = len1; j < len; j++)
               _fmpz_demote(p1 + j);
            len = len1;

            goto cleanup;
         }
      }
   }

   /* not done, coeffs small, intermediate computations fit three words */
   if (len == len1 && small)
   {
      ulong * p2 = (ulong *) TMP_ALLOC(3*prod*sizeof(ulong));

      for (i = 0; i < 3*prod; i++)
         p2[i] = 0;

      /* poly2 to array format */
      _fmpz_mpoly_to_ulong_array(p2, poly2, exp2, len2);

      /* for each term of poly2 array relevant to exact quotient */
      for (i = prod - 1; i >= max3; i--)
      {
         ulong * ptr = p2 + 3*i;
         ulong p[3];

         /* if coeff is nonzero */
         if (ptr[0] != 0 || ptr[1] != 0 || ptr[2] != 0)
         {
            if (0 > (slong) ptr[2])
               mpn_neg(p, ptr, 3);
            else
               flint_mpn_copyi(p, ptr, 3);

            /* check quotient won't overflow a word */
            if (p[2] > 0 || u3 <= p[1] || (u3 == 0 && 0 > (slong) p[0])) /* quotient too large */
            {
               for (j = len1; j < len; j++)
                  _fmpz_demote(p1 + j);
               len = len1;

               goto big;
            }

            /* quotient and remainder of coeffs */
            sdiv_qrnnd(q, r, ptr[1], ptr[0], c3);

            /* check coefficient is small, else restart with multiprec code */
            if (COEFF_IS_MPZ(FLINT_ABS(q)))
            {
               for (j = len1; j < len; j++) /* quotient too large */
                  _fmpz_demote(p1 + j);
               len = len1;

               goto big;
            }

            /* check coeff and monomial quotients were exact */
            if (r != 0 || /* not an exact division */
               !mpoly_monomial_divides_tight(i, max3, prods, num))
            {
               for (j = len1; j < len; j++)
                  _fmpz_demote(p1 + j);
               len = len1;

               goto cleanup;
            }

            /* submul q*poly3 */
            _fmpz_mpoly_submul_array1_slong_1(p2, q, i - max3,
                                                            poly3, exp3, len3);

            /* reallocate quotient */
            _fmpz_mpoly_fit_length(&p1, &e1, alloc, len + 1, 1);

            /* write quotient term */
            fmpz_set_si(p1 + len, q);
            e1[len++] = i - max3;
         }
      }

      /* check there are no nonzero terms left in array */
      for ( ; i >= min3; i--)
      {
         ulong * ptr = p2 + 3*i;

         /* if coeff nonzero */
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

   /* if not done, use multiprecision coeffs instead */
   if (len == len1)
   {
      fmpz * p2 = (fmpz *) TMP_ALLOC(prod*sizeof(fmpz));
      fmpz_t fq, fr;

      fmpz_init(fq);
      fmpz_init(fr);

      for (i = 0; i < prod; i++)
         fmpz_init(p2 + i);

      /* poly2 to array format */
      _fmpz_mpoly_to_fmpz_array(p2, poly2, exp2, len2);

      /* for each term of poly2 array relevant to exact quotient */
      for (i = prod - 1; i >= max3; i--)
      {
         /* if coeff is nonzero */
         if (!fmpz_is_zero(p2 + i))
         {
            /* quotient and remainder of coeffs */
            fmpz_fdiv_qr(fq, fr, p2 + i, poly3);

            /* check coeff and monomial quotients were exact */
            if (!fmpz_is_zero(fr) || /* not an exact division */
               !mpoly_monomial_divides_tight(i, max3, prods, num))
            {
               for (j = len1; j < len; j++)
                  _fmpz_demote(p1 + j);
               len = len1;

               goto cleanup2;
            }

            /* submul q*poly3 */
            _fmpz_mpoly_submul_array1_fmpz_1(p2, fq, i - max3,
                                                            poly3, exp3, len3);

            /* reallocate quotient */
            _fmpz_mpoly_fit_length(&p1, &e1, alloc, len + 1, 1);

            /* write quotient term */
            fmpz_set(p1 + len, fq);

            e1[len++] = i - max3;
         }
      }

      /* check there are no nonzero terms left in array */
      for ( ; i >= min3; i--)
      {
         /* if coeff nonzero */
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
   Use dense array exact division to set poly1 to poly2/poly3 in num + 1
   variables, given a list of multipliers to tightly pack exponents and a
   number of bits for the fields of the exponents of the result, assuming no
   aliasing. Classical exact division in main variable, array multiplication
   (submul) for multivariate coefficients in remaining num variables.
   The array "mults" is a list of bases to be used in encoding the array
   indices from the exponents. The function reallocates its output and returns
   the length of its output if the quotient is exact, or zero if not. It is
   assumed that poly2 is not zero.
*/
slong _fmpz_mpoly_divides_array_chunked(fmpz ** poly1, ulong ** exp1,
         slong * alloc, const fmpz * poly2, const ulong * exp2, slong len2,
                        const fmpz * poly3, const ulong * exp3, slong len3,
                                          slong * mults, slong num, slong bits)
{
   slong i, j, k, prod, len = 0, l1, l2, l3;
   slong bits1, bits2, bits3 = 0, tlen, talloc, skip, max_exp;
   slong shift = bits*(num);
   slong * i1, * i2, * i3, * n1, * n2, * n3;
   slong * b1, * b3, * maxb1, * maxb3, * max_exp1, * max_exp3;
   ulong * e2, * e3, * texp, * p2;
   fmpz * temp;
   int small;
   TMP_INIT;

   /*
      compute products 1, b0, b0*b1, b0*b1*b2 ...
      from list of bases b0, b1, b2, ...
   */
   prod = 1;
   for (i = 0; i < num; i++)
      prod *= mults[i];

   /* lengths of poly2, poly3 in chunks, and poly1 assuming exact division */
   l2 = 1 + (slong) (exp2[0] >> shift);
   l3 = 1 + (slong) (exp3[0] >> shift);

   l1 = l2 - l3 + 1;

   TMP_START;

   /* indices and lengths of coefficients/chunks of polys in main variable */

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
   max_exp1 = (slong *) TMP_ALLOC(l1*sizeof(slong));
   max_exp3 = (slong *) TMP_ALLOC(l3*sizeof(slong));

   mpoly_main_variable_terms1(i2, n2, exp2, l2, len2, num + 1, num + 1, bits);
   mpoly_main_variable_terms1(i3, n3, exp3, l3, len3, num + 1, num + 1, bits);

   /* alloc space for copy of coeff/chunk of poly2 */

   temp = (fmpz *) flint_calloc(n2[0] + 1, sizeof(fmpz));
   texp = (ulong *) flint_malloc((n2[0] + 1)*sizeof(ulong));
   talloc = n2[0] + 1;

   /* figure out how many trailing zero chunks poly3 has */
   skip = 0;
   while (n3[l3 - skip - 1] == 0)
   {
      /* check poly2 has at least as many trailing zero chunks */
      if (n2[l2 - skip - 1] != 0) /* not an exact quotient */
      {
         goto cleanup;
      }
      skip++;
   }

   /* work out max bits for each coeff/chunk, optimal bits */

   for (i = 0; i < l3; i++)
   {
      _fmpz_vec_sum_max_bits(&b3[i], &maxb3[i], poly3+i3[i], n3[i]);

      if (bits3 < maxb3[i])
         bits3 = maxb3[i];
   }

   /* pack input coefficients tightly */

   e2 = (ulong *) TMP_ALLOC(len2*sizeof(ulong));
   e3 = (ulong *) TMP_ALLOC(len3*sizeof(ulong));

   mpoly_pack_monomials_tight(e2, exp2, len2, mults, num, bits);
   mpoly_pack_monomials_tight(e3, exp3, len3, mults, num, bits);

   /* work out maximum packed exponent for each chunk */
   for (i = 0; i < l3; i++)
   {
      max_exp = 0;

      for (j = 0; j < n3[i]; j++)
      {
         if (e3[i3[i] + j] > max_exp)
            max_exp = e3[i3[i] + j];
      }

      max_exp3[i] = max_exp;
   }

   /* bound poly2 coeffs and check input/output coeffs likely small */
   bits2 = _fmpz_vec_max_bits(poly2, len2);
   /* we assume a bound of SMALL_FMPZ_BITCOUNT_MAX for coefficients of the quotient */
   bits1 = FLINT_ABS(bits3) + FLINT_BITS + FLINT_BIT_COUNT(len3) - 2;

   small = FLINT_ABS(bits2) <= bits1 && FLINT_ABS(bits3) <= SMALL_FMPZ_BITCOUNT_MAX;

   /* enough space for three words per coeff, even if only one or two needed */
   p2 = (ulong *) TMP_ALLOC(3*prod*sizeof(ulong));

   /* coefficients likely to be small */
   if (small)
   {
      /* for each chunk of poly2 */
      for (i = 0; i < l2 - skip; i++)
      {
         slong num1 = 0;
         bits1 = 0;

         /* if there are already quotient terms */
         if (i != 0)
         {
            /* compute bound on intermediate computations a - q*b */
            for (j = 0; j < i && j < l1; j++)
            {
               k = i - j;

               if (k < l3)
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

         /* intermediate computations fit in one word */
         if (bits1 <= FLINT_BITS)
         {
            for (j = 0; j < prod; j++)
               p2[j] = 0;

            /* convert relevant coeff/chunk of poly2 to array format */
            _fmpz_mpoly_to_ulong_array1(p2, poly2 + i2[i], e2 + i2[i], n2[i]);

            /* submuls */

            for (j = 0; j < i && j < l1; j++)
            {
               k = i - j;

               if (k < l3)
               {
                  /* check for exponent overflow: not exact division */
                  if (max_exp1[j] + max_exp3[k] >= prod)
                  {
                     for (j = 0; j < len; j++)
                        _fmpz_demote((*poly1) + j);
                     len = 0;

                     goto cleanup;
                  }

                  _fmpz_mpoly_submul_array1_slong1(p2, (*poly1) + i1[j],
                     (*exp1) + i1[j], n1[j], poly3 + i3[k], e3 + i3[k], n3[k]);
               }
            }

            /* convert chunk from array format */
            tlen = _fmpz_mpoly_from_ulong_array1(&temp, &texp, &talloc,
                                                      p2, mults, num, bits, 0);
         } else if (bits1 <= 2*FLINT_BITS) /* intermed comps fit two words */
         {
            for (j = 0; j < 2*prod; j++)
               p2[j] = 0;

            /* convert relevant coeff/chunk of poly2 to array format */
            _fmpz_mpoly_to_ulong_array2(p2, poly2 + i2[i], e2 + i2[i], n2[i]);

            /* submuls */

            for (j = 0; j < i && j < l1; j++)
            {
               k = i - j;

               if (k < l3)
               {
                  /* check for exponent overflow: not exact division */
                  if (max_exp1[j] + max_exp3[k] >= prod)
                  {
                     for (j = 0; j < len; j++)
                        _fmpz_demote((*poly1) + j);
                     len = 0;

                     goto cleanup;
                  }

                  _fmpz_mpoly_submul_array1_slong2(p2, (*poly1) + i1[j],
                     (*exp1) + i1[j], n1[j], poly3 + i3[k], e3 + i3[k], n3[k]);
               }
            }

            /* convert chunk from array format */
            tlen = _fmpz_mpoly_from_ulong_array2(&temp, &texp, &talloc,
                                                      p2, mults, num, bits, 0);
         } else /* intermed comps fit three words */
         {
            for (j = 0; j < 3*prod; j++)
               p2[j] = 0;

            /* convert relevant coeff/chunk of poly2 to array format */
            _fmpz_mpoly_to_ulong_array(p2, poly2 + i2[i], e2 + i2[i], n2[i]);

            /* submuls */

            for (j = 0; j < i && j < l1; j++)
            {
               k = i - j;

               if (k < l3)
               {
                  /* check for exponent overflow: not exact division */
                  if (max_exp1[j] + max_exp3[k] >= prod)
                  {
                     for (j = 0; j < len; j++)
                        _fmpz_demote((*poly1) + j);
                     len = 0;

                     goto cleanup;
                  }

                  _fmpz_mpoly_submul_array1_slong(p2, (*poly1) + i1[j],
                     (*exp1) + i1[j], n1[j], poly3 + i3[k], e3 + i3[k], n3[k]);
               }
            }

            /* convert chunk from array format */
            tlen = _fmpz_mpoly_from_ulong_array(&temp, &texp, &talloc,
                                                      p2, mults, num, bits, 0);
         }

         /* for terms where there may be a nonzero quotient if exact */
         if (i < l1)
         {
            /* tightly pack chunk exponents */
            mpoly_pack_monomials_tight(texp, texp, tlen, mults, num, bits);

            /* set starting index for quotient chunk we are about to compute */
            i1[i] = len;

            /* if chunk is nonzero */
            if (tlen != 0)
            {
               /* compute quotient chunk and set length of quotient chunk */
               n1[i] = _fmpz_mpoly_divides_array_tight(poly1,
                                   exp1, alloc, len, temp, texp, tlen,
                        poly3 + i3[0], e3 + i3[0], n3[0], mults, num);

               /* check quotient was exact */
               if (n1[i] == 0) /* not an exact division */
               {
                  for (j = 0; j < len; j++)
                     _fmpz_demote((*poly1) + j);
                  len = 0;

                  goto cleanup;
               }
            } else
               n1[i] = 0;

            /* compute maximum packed exponent for chunk */
            max_exp = 0;

            for (j = 0; j < n1[i]; j++)
            {
               if ((*exp1)[i1[i] + j] > max_exp)
                  max_exp = (*exp1)[i1[i] + j];
            }

            max_exp1[i] = max_exp;

            /* check the quotient chunk didn't have large coefficients */
            if (FLINT_ABS(_fmpz_vec_max_bits((*poly1) + len, n1[i])) >
                                                             SMALL_FMPZ_BITCOUNT_MAX)
            {
               for (j = 0; j < len; j++)
                  _fmpz_demote((*poly1) + j);
               len = 0;

               goto big;
            }

            /* abs bound and sum of abs vals of coeffs of quotient chunk */
            _fmpz_vec_sum_max_bits(&b1[i], &maxb1[i], (*poly1)+i1[i], n1[i]);

            /* update length of output quotient poly */
            len += n1[i];
         } else /* should be zero quotient, check coefficient is zero */
         {
            /* for each coeff in chunk */
            for (j = 0; j < tlen; j++)
            {
               /* if coeff is nonzero */
               if (!fmpz_is_zero(temp + j)) /* not an exact division */
               {
                  for (j = 0; j < len; j++)
                     _fmpz_demote((*poly1) + j);
                  len = 0;

                  goto cleanup;
               }
            }
         }
      }
   }

big:

   /* if not done, use multiprecision coeffs instead */
   if (len == 0)
   {
      fmpz * p2 = (fmpz *) TMP_ALLOC(prod*sizeof(fmpz));

      for (j = 0; j < prod; j++)
            fmpz_init(p2 + j);

      /* for each chunk of poly2 */
      for (i = 0; i < l2 - skip; i++)
      {
         for (j = 0; j < prod; j++)
            fmpz_zero(p2 + j);

         /* convert relevant coeff/chunk of poly2 to array format */
         _fmpz_mpoly_to_fmpz_array(p2, poly2 + i2[i], e2 + i2[i], n2[i]);

         /* submuls */

         for (j = 0; j < i && j < l1; j++)
         {
            k = i - j;

            if (k < l3)
            {
               /* check for exponent overflow: not exact division */
               if (max_exp1[j] + max_exp3[k] >= prod)
               {
                  for (j = 0; j < len; j++)
                     _fmpz_demote((*poly1) + j);
                  len = 0;

                  goto cleanup;
               }

               _fmpz_mpoly_submul_array1_fmpz(p2, (*poly1) + i1[j],
                     (*exp1) + i1[j], n1[j], poly3 + i3[k], e3 + i3[k], n3[k]);
            }
         }

         /* convert chunk from array format */
         tlen = _fmpz_mpoly_from_fmpz_array(&temp, &texp, &talloc,
                                                      p2, mults, num, bits, 0);

         /* for terms where there may be a nonzero quotient if exact */
         if (i < l1)
         {
            /* tightly pack chunk exponents */
            mpoly_pack_monomials_tight(texp, texp, tlen, mults, num, bits);

            /* set starting index of quotient chunk we are about to compute */
            i1[i] = len;

            /* if chunk is nonzero */
            if (tlen != 0)
            {
               /* compute quotient chunk and set length of quotient chunk */
               n1[i] = _fmpz_mpoly_divides_array_tight(poly1,
                                   exp1, alloc, len, temp, texp, tlen,
                        poly3 + i3[0], e3 + i3[0], n3[0], mults, num);

               /* check quotient was exact */
               if (n1[i] == 0) /* not an exact division */
               {
                  for (j = 0; j < len; j++)
                     _fmpz_demote((*poly1) + j);
                  len = 0;

                  goto cleanup2;
               }
            } else
               n1[i] = 0;

            /* compute maximum packed exponent for chunk */
            max_exp = 0;

            for (j = 0; j < n1[i]; j++)
            {
               if ((*exp1)[i1[i] + j] > max_exp)
                  max_exp = (*exp1)[i1[i] + j];
            }

            max_exp1[i] = max_exp;

            /* update length of output quotient poly */
            len += n1[i];
         } else /* should be zero quotient, check coefficient is zero */
         {
            /* for each coeff in chunk */
            for (j = 0; j < tlen; j++)
            {
               if (!fmpz_is_zero(temp + j)) /* not an exact division */
               {
                  for (j = 0; j < len; j++)
                     _fmpz_demote((*poly1) + j);
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

   /* if quotient was exact */
   if (len != 0)
   {
      /* unpack monomials of quotient */
      mpoly_unpack_monomials_tight((*exp1), (*exp1), len, mults, num, bits);

      /* put main variable back in quotient */
      for (i = 0; i < l1; i++)
      {
         for (j = 0; j < n1[i]; j++)
            (*exp1)[i1[i] + j] += ((l1 - i - 1) << shift);
      }
   }

cleanup:

   for (j = 0; j < talloc; j++)
      fmpz_clear(temp + j);

   flint_free(temp);
   flint_free(texp);

   TMP_END;

   return len;
}

/*
   Use dense array exact division to set poly1 to poly2/poly3 in num variables,
   given a list of multipliers to tightly pack exponents and a number of bits
   for the fields of the exponents of the result, assuming no aliasing. The
   array "mults" is a list of bases to be used in encoding the array indices
   from the exponents. The function reallocates its output and returns the
   length of its output if the quotient is exact, or zero if not. It is assumed
   that poly2 is not zero.
*/
slong _fmpz_mpoly_divides_array(fmpz ** poly1, ulong ** exp1, slong * alloc,
                           const fmpz * poly2, const ulong * exp2, slong len2,
                           const fmpz * poly3, const ulong * exp3, slong len3,
                                          slong * mults, slong num, slong bits)
{
   slong i;
   ulong * e2, * e3;
   slong len, prod;
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
      return _fmpz_mpoly_divides_array_chunked(poly1, exp1, alloc,
                   poly2, exp2, len2, poly3, exp3, len3, mults, num - 1, bits);

   TMP_START;

   e2 = (ulong *) TMP_ALLOC(len2*sizeof(ulong));
   e3 = (ulong *) TMP_ALLOC(len3*sizeof(ulong));

   /* pack input exponents tightly with mixed bases specified by "mults" */

   mpoly_pack_monomials_tight(e2, exp2, len2, mults, num, bits);
   mpoly_pack_monomials_tight(e3, exp3, len3, mults, num, bits);

   /* do exact quotient with divisibility test on tightly packed polys */
   len = _fmpz_mpoly_divides_array_tight(poly1, exp1,
                      alloc, 0,  poly2, e2, len2, poly3, e3, len3, mults, num);

   /* unpack output quotient exponents */
   mpoly_unpack_monomials_tight((*exp1), (*exp1), len, mults, num, bits);

   TMP_END;

   return len;
}

/*
   Return 1 if exact quotient and set poly1 to poly2/poly3, else return 0
   if not exact quotient or -1 if array method not suitable.
*/
int fmpz_mpoly_divides_array(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                          const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx)
{
   slong i, bits, exp_bits, N, len = 0, array_size;
   ulong max, * max_fields1, * max_fields2, * max_fields3;
   ulong * exp2 = poly2->exps, * exp3 = poly3->exps;
   int free2 = 0, free3 = 0;
   int res = -1;

   TMP_INIT;

   /* check divisor is not zero */
   if (poly3->length == 0)
      flint_throw(FLINT_DIVZERO, "Divide by zero in fmpz_mpoly_divides_array");

   /* dividend is zero */
   if (poly2->length == 0)
   {
      fmpz_mpoly_zero(poly1, ctx);

      return 1;
   }

   TMP_START;

    /* compute maximum fields */
    max_fields1 = (ulong *) TMP_ALLOC(ctx->minfo->nfields*sizeof(ulong));
    max_fields2 = (ulong *) TMP_ALLOC(ctx->minfo->nfields*sizeof(ulong));
    max_fields3 = (ulong *) TMP_ALLOC(ctx->minfo->nfields*sizeof(ulong));
    mpoly_max_fields_ui_sp(max_fields2, poly2->exps, poly2->length,
                                                      poly2->bits, ctx->minfo);
    mpoly_max_fields_ui_sp(max_fields3, poly3->exps, poly3->length,
                                                      poly3->bits, ctx->minfo);
    max = 0;
    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        if (max_fields2[i] > max)
            max = max_fields2[i];
        /*
            cannot be an exact division if variable in dividend has smaller degree
            than corresponding variable in divisor
        */
        if (max_fields2[i] < max_fields3[i])
        {
            res = 0;
            goto cleanup;
        }
    }

    /* compute number of bits required for output exponents */
    bits = FLINT_BIT_COUNT(max);
    exp_bits = FLINT_MAX(WORD(8), bits + 1);
    exp_bits = FLINT_MAX(exp_bits, poly2->bits);
    exp_bits = FLINT_MAX(exp_bits, poly3->bits);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);

    N = mpoly_words_per_exp(exp_bits, ctx->minfo);

   /* array division expects each exponent vector in one word */
   /* current code is wrong for reversed orderings */
   if (N != 1 || mpoly_ordering_isrev(ctx->minfo))
      goto cleanup;

   /* compute bounds on output exps, used as mixed bases for packing exps */
   array_size = 1;
   for (i = 0; i < ctx->minfo->nfields - 1; i++)
   {
      max_fields2[i]++;
      array_size *= max_fields2[i];
   }
   max_fields2[ctx->minfo->nfields - 1]++;

   /* if exponents too large for array multiplication, exit silently */
   if (array_size > MAX_ARRAY_SIZE)
      goto cleanup;

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

   /* handle aliasing and do array division */
   if (poly1 == poly2 || poly1 == poly3)
   {
      fmpz_mpoly_t temp;

      fmpz_mpoly_init2(temp, poly2->length/poly3->length + 1, ctx);
      fmpz_mpoly_fit_bits(temp, exp_bits, ctx);
      temp->bits = exp_bits;

      len = _fmpz_mpoly_divides_array(&temp->coeffs, &temp->exps,
                     &temp->alloc, poly2->coeffs, exp2, poly2->length,
                                       poly3->coeffs, exp3, poly3->length,
                         (slong *) max_fields2, ctx->minfo->nfields, exp_bits);

      fmpz_mpoly_swap(temp, poly1, ctx);

      fmpz_mpoly_clear(temp, ctx);
   } else
   {
      fmpz_mpoly_fit_length(poly1, poly2->length/poly3->length + 1, ctx);
      fmpz_mpoly_fit_bits(poly1, exp_bits, ctx);
      poly1->bits = exp_bits;

      len = _fmpz_mpoly_divides_array(&poly1->coeffs, &poly1->exps,
                     &poly1->alloc, poly2->coeffs, exp2, poly2->length,
                                        poly3->coeffs, exp3, poly3->length,
                         (slong *) max_fields2, ctx->minfo->nfields, exp_bits);
   }

   _fmpz_mpoly_set_length(poly1, len, ctx);

   /*
      Check that the maximum degrees for each variable are correct. This can
      only happen if there was an overflow in the tightly packed exponents,
      which can only happen if the division was not exact. If no overflow
      is detected here, none happened and the division was provably exact.
   */
   if (len != 0)
   {
      mpoly_max_fields_ui_sp(max_fields1, poly1->exps, poly1->length,
                                                      poly1->bits, ctx->minfo);

      for (i = 0; i < ctx->minfo->nfields; i++)
      {
         /* we incremented max_fields2[i] by 1 previously, so add 1 here */
         if (max_fields2[i] != max_fields1[i] + max_fields3[i] + 1)
         {
            fmpz_mpoly_zero(poly1, ctx);
            len = 0;
            break;
         }
      }
   }

   /* len will be nonzero if quotient was exact */
   res = len != 0;

   if (free2)
      flint_free(exp2);

   if (free3)
      flint_free(exp3);

cleanup:

   TMP_END;

   return res;
}
