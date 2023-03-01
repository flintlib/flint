/*
    Copyright (C) 2016-2017 William Hart
    Copyright (C) 2017-2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"
#include "thread_pool.h"

int fmpz_mpoly_divides(
    fmpz_mpoly_t Q,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    thread_pool_handle * handles;
    slong num_handles;
    int divides;
    slong thread_limit;

    thread_limit = A->length/1024;

    if (B->length < 2 || A->length < 2)
    {
        if (B->length == 0)
        {
            flint_throw(FLINT_DIVZERO, "Divide by zero in fmpz_mpoly_divides");
        }

        if (A->length == 0)
        {
            fmpz_mpoly_zero(Q, ctx);
            return 1;
        }

        return fmpz_mpoly_divides_monagan_pearce(Q, A, B, ctx);
    }

    num_handles = flint_request_threads(&handles, thread_limit);

    divides = (num_handles > 0)
            ? _fmpz_mpoly_divides_heap_threaded_pool(Q, A, B, ctx,
                                                          handles, num_handles)
            : fmpz_mpoly_divides_monagan_pearce(Q, A, B, ctx);

    flint_give_back_threads(handles, num_handles);

    return divides;
}

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

/*
    a thread safe mpoly supports three mutating operations
    - init from an array of terms
    - append an array of terms
    - clear out contents to a normal mpoly
*/
typedef struct _fmpz_mpoly_ts_struct
{
    fmpz * volatile coeffs; /* this is coeff_array[idx] */
    ulong * volatile exps;       /* this is exp_array[idx] */
    volatile slong length;
    slong alloc;
    flint_bitcnt_t bits;
    flint_bitcnt_t idx;    
    ulong * exp_array[FLINT_BITS];
    fmpz * coeff_array[FLINT_BITS];
} fmpz_mpoly_ts_struct;

typedef fmpz_mpoly_ts_struct fmpz_mpoly_ts_t[1];

/* Bcoeff is changed */
void fmpz_mpoly_ts_init(fmpz_mpoly_ts_t A,
                              fmpz * Bcoeff, ulong * Bexp, slong Blen,
                                                    flint_bitcnt_t bits, slong N)
{
    slong i;
    flint_bitcnt_t idx = FLINT_BIT_COUNT(Blen);
    idx = (idx <= 8) ? 0 : idx - 8;
    for (i = 0; i < FLINT_BITS; i++)
    {
        A->exp_array[i] = NULL;
        A->coeff_array[i] = NULL;
    }
    A->bits = bits;
    A->idx = idx;
    A->alloc = WORD(256) << idx;
    A->exps = A->exp_array[idx]
            = (ulong *) flint_malloc(N*A->alloc*sizeof(ulong));
    A->coeffs = A->coeff_array[idx]
              = (fmpz *) flint_calloc(A->alloc, sizeof(fmpz));
    A->length = Blen;
    for (i = 0; i < Blen; i++)
    {
        fmpz_swap(A->coeffs + i, Bcoeff + i);
        mpoly_monomial_set(A->exps + N*i, Bexp + N*i, N);
    }
}

void fmpz_mpoly_ts_print(const fmpz_mpoly_ts_t B, const char ** x,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_t A;
    A->length = B->length;
    A->alloc = B->alloc;
    A->coeffs = B->coeffs;
    A->exps = B->exps;
    A->bits = B->bits;
    fmpz_mpoly_print_pretty(A, x, ctx);

    fmpz_mpoly_assert_canonical(A, ctx);
}

void fmpz_mpoly_ts_clear(fmpz_mpoly_ts_t A)
{
    slong i;

    for (i = 0; i < A->length; i++)
    {
        fmpz_clear(A->coeffs + i);
    }

    for (i = 0; i < FLINT_BITS; i++)
    {
        if (A->exp_array[i] != NULL)
        {
            FLINT_ASSERT(A->coeff_array[i] != NULL);
            flint_free(A->coeff_array[i]);    
            flint_free(A->exp_array[i]);
        }
    }
}

void fmpz_mpoly_ts_clear_poly(fmpz_mpoly_t Q, fmpz_mpoly_ts_t A)
{
    if (Q->alloc != 0)
    {
        slong i;

        FLINT_ASSERT(Q->exps != NULL);
        FLINT_ASSERT(Q->coeffs != NULL);

        for (i = 0; i < Q->alloc; i++)
        {
            _fmpz_demote(Q->coeffs + i);
        }
        flint_free(Q->exps);
        flint_free(Q->coeffs);
    }

    Q->exps = A->exps;
    Q->coeffs = A->coeffs;
    Q->bits = A->bits;
    Q->alloc = A->alloc;
    Q->length = A->length;

    A->length = 0;
    A->coeff_array[A->idx] = NULL;
    A->exp_array[A->idx] = NULL;
    fmpz_mpoly_ts_clear(A);
}


/* put B on the end of A - Bcoeff is changed*/
void fmpz_mpoly_ts_append(fmpz_mpoly_ts_t A,
                             fmpz * Bcoeff, ulong * Bexps, slong Blen, slong N)
{
/* TODO: this needs barriers on non-x86 */

    slong i;
    ulong * oldexps = A->exps;
    fmpz * oldcoeffs = A->coeffs;
    slong oldlength = A->length;
    slong newlength = A->length + Blen;

    if (newlength <= A->alloc)
    {
        /* write new terms first */
        for (i = 0; i < Blen; i++)
        {
            fmpz_swap(oldcoeffs + oldlength + i, Bcoeff + i);
            mpoly_monomial_set(oldexps + N*(oldlength + i), Bexps + N*i, N);
        }
    }
    else
    {
        slong newalloc;
        ulong * newexps;
        fmpz * newcoeffs;
        flint_bitcnt_t newidx;
        newidx = FLINT_BIT_COUNT(newlength - 1);
        newidx = (newidx > 8) ? newidx - 8 : 0;
        FLINT_ASSERT(newidx > A->idx);

        newalloc = UWORD(256) << newidx;
        FLINT_ASSERT(newlength <= newalloc);
        newexps = A->exp_array[newidx]
                = (ulong *) flint_malloc(N*newalloc*sizeof(ulong));
        newcoeffs = A->coeff_array[newidx]
                  = (fmpz *) flint_calloc(newalloc, sizeof(fmpz));

        for (i = 0; i < oldlength; i++)
        {
            newcoeffs[i] = oldcoeffs[i]; /* just copy the bits */
            mpoly_monomial_set(newexps + N*i, oldexps + N*i, N);
        }
        for (i = 0; i < Blen; i++)
        {
            fmpz_swap(newcoeffs + oldlength + i, Bcoeff + i);
            mpoly_monomial_set(newexps + N*(oldlength + i), Bexps + N*i, N);
        }

        A->alloc = newalloc;
        A->exps = newexps;
        A->coeffs = newcoeffs;
        A->idx = newidx;

        /* do not free oldcoeff/exps as other threads may be using them */
    }

    /* update length at the very end */
    A->length = newlength;
}


/*
    a chunk holds an exponent range on the dividend
*/
typedef struct _divides_heap_chunk_struct
{
    fmpz_mpoly_t polyC;
    struct _divides_heap_chunk_struct * next;
    ulong * emin;
    ulong * emax;
    slong startidx;
    slong endidx;
    int upperclosed;
    volatile int lock;
    volatile int producer;
    volatile slong ma;
    volatile slong mq;
    int Cinited;
} divides_heap_chunk_struct;

typedef divides_heap_chunk_struct divides_heap_chunk_t[1];

/*
    the base struct includes a linked list of chunks
*/
typedef struct
{
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
#endif
    divides_heap_chunk_struct * head;
    divides_heap_chunk_struct * tail;
    divides_heap_chunk_struct * volatile cur;
    fmpz_mpoly_t polyA;
    fmpz_mpoly_t polyB;
    fmpz_mpoly_ts_t polyQ;
    const fmpz_mpoly_ctx_struct * ctx;
    slong length;
    slong N;
    flint_bitcnt_t bits;
    slong polyBcoeff_bits;
    ulong * cmpmask;
    int failed;
} divides_heap_base_struct;

typedef divides_heap_base_struct divides_heap_base_t[1];

/*
    the worker struct has a big chunk of memory in the stripe_t
    and two polys for work space
*/
typedef struct _worker_arg_struct
{
    divides_heap_base_struct * H;
    fmpz_mpoly_stripe_t S;
    fmpz_mpoly_t polyT1;
    fmpz_mpoly_t polyT2;
} worker_arg_struct;

typedef worker_arg_struct worker_arg_t[1];


static void divides_heap_base_init(divides_heap_base_t H)
{
    H->head = NULL;
    H->tail = NULL;
    H->cur = NULL;
    H->ctx = NULL;
    H->length = 0;
    H->N = 0;
    H->bits = 0;
    H->cmpmask = NULL;
}

static void divides_heap_chunk_clear(divides_heap_chunk_t L, divides_heap_base_t H)
{
    if (L->Cinited)
    {
        fmpz_mpoly_clear(L->polyC, H->ctx);
    }
}
static int divides_heap_base_clear(fmpz_mpoly_t Q, divides_heap_base_t H)
{
    divides_heap_chunk_struct * L = H->head;
    while (L != NULL)
    {
        divides_heap_chunk_struct * nextL = L->next;
        divides_heap_chunk_clear(L, H);
        flint_free(L);
        L = nextL;
    }
    H->head = NULL;
    H->tail = NULL;
    H->cur = NULL;
    H->ctx = NULL;
    H->length = 0;
    H->N = 0;
    H->bits = 0;
    H->cmpmask = NULL;

    if (H->failed)
    {
        fmpz_mpoly_zero(Q, H->ctx);
        fmpz_mpoly_ts_clear(H->polyQ);
        return 0;
    }
    else
    {
        fmpz_mpoly_ts_clear_poly(Q, H->polyQ);
        return 1;
    }
}

static void divides_heap_base_add_chunk(divides_heap_base_t H, divides_heap_chunk_t L)
{
    L->next = NULL;

    if (H->tail == NULL)
    {
        FLINT_ASSERT(H->head == NULL);
        H->tail = L;
        H->head = L;
    }
    else
    {
        divides_heap_chunk_struct * tail = H->tail;
        FLINT_ASSERT(tail->next == NULL);
        tail->next = L;
        H->tail = L;
    }
    H->length++;
}


/*
    A = D - (a stripe of B * C)
    S->startidx and S->endidx are assumed to be correct
        that is, we expect and successive calls to keep
            B decreasing
            C the same

    saveD: 0 means we can modify coeffs of input D
           1 means we must not modify coeffs of input D
*/
slong _fmpz_mpoly_mulsub_stripe1(fmpz ** A_coeff, ulong ** A_exp, slong * A_alloc,
                         const fmpz * Dcoeff, const ulong * Dexp, slong Dlen, int saveD,
                         const fmpz * Bcoeff, const ulong * Bexp, slong Blen,
                         const fmpz * Ccoeff, const ulong * Cexp, slong Clen,
                                                   const fmpz_mpoly_stripe_t S)
{
    int upperclosed;
    slong startidx, endidx;
    ulong prev_startidx;
    ulong maskhi = S->cmpmask[0];
    ulong emax = S->emax[0];
    ulong emin = S->emin[0];
    slong i, j;
    slong next_loc = Blen + 4;   /* something bigger than heap can ever be */
    slong heap_len = 1; /* heap zero index unused */
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    slong Di;
    slong Alen;
    slong Aalloc = *A_alloc;
    fmpz * Acoeff = *A_coeff;
    ulong * Aexp = *A_exp;
    ulong exp;
    slong * ends;
    ulong texp;
    slong * hind;
    int small;
    ulong acc_sm[3], pp0, pp1;

    FLINT_ASSERT(S->N == 1);

    i = 0;
    hind = (slong *)(S->big_mem + i);
    i += Blen*sizeof(slong);
    ends = (slong *)(S->big_mem + i);
    i += Blen*sizeof(slong);
    store = store_base = (slong *) (S->big_mem + i);
    i += 2*Blen*sizeof(slong);
    heap = (mpoly_heap1_s *)(S->big_mem + i);
    i += (Blen + 1)*sizeof(mpoly_heap1_s);
    chain = (mpoly_heap_t *)(S->big_mem + i);
    i += Blen*sizeof(mpoly_heap_t);
    FLINT_ASSERT(i <= S->big_mem_alloc);

    startidx = *S->startidx;
    endidx = *S->endidx;
    upperclosed = S->upperclosed;
    emax = S->emax[0];
    emin = S->emin[0];

    /* put all the starting nodes on the heap */
    prev_startidx = -UWORD(1);
    for (i = 0; i < Blen; i++)
    {
        if (startidx < Clen)
        {
            texp = Bexp[i] + Cexp[startidx];
            FLINT_ASSERT(mpoly_monomial_cmp1(emax, texp, maskhi)
                                                               > -upperclosed);
        }
        while (startidx > 0)
        {
            texp = Bexp[i] + Cexp[startidx - 1];
            if (mpoly_monomial_cmp1(emax, texp, maskhi) <= -upperclosed)
            {
                break;
            }
            startidx--;
        }

        if (endidx < Clen)
        {
            texp = Bexp[i] + Cexp[endidx];
            FLINT_ASSERT(mpoly_monomial_cmp1(emin, texp, maskhi) > 0);
        }
        while (endidx > 0)
        {
            texp = Bexp[i] + Cexp[endidx - 1];
            if (mpoly_monomial_cmp1(emin, texp, maskhi) <= 0)
            {
                break;
            }
            endidx--;
        }

        ends[i] = endidx;

        hind[i] = 2*startidx + 1;

        if (  (startidx < endidx)
           && (((ulong)startidx) < prev_startidx) 
           )
        {
            x = chain + i;
            x->i = i;
            x->j = startidx;
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;
            _mpoly_heap_insert1(heap, Bexp[x->i] + Cexp[x->j], x,
                                      &next_loc, &heap_len, maskhi);
        }

        prev_startidx = startidx;
    }

    /* set the indices for the next time mul is called */
    *S->startidx = startidx;
    *S->endidx = endidx;

    /* whether input coeffs are small, thus output coeffs fit in three words */
    FLINT_ASSERT(ends[0] >= startidx);
    FLINT_ASSERT(S->coeff_bits == FLINT_ABS(_fmpz_vec_max_bits(Ccoeff, Clen)));
    small = S->coeff_bits <= SMALL_FMPZ_BITCOUNT_MAX
            && _fmpz_mpoly_fits_small(Bcoeff, Blen)
            && FLINT_ABS(_fmpz_vec_max_bits(Dcoeff, Dlen)) < 3*FLINT_BITS - 3;

    Alen = 0;
    Di = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;

        while (Di < Dlen && mpoly_monomial_gt1(Dexp[Di], exp, maskhi))
        {
            _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, 1);
            Aexp[Alen] = Dexp[Di];
            if (saveD)
                fmpz_set(Acoeff + Alen, Dcoeff + Di);
            else
                fmpz_swap(Acoeff + Alen, (fmpz *)(Dcoeff + Di));
            Alen++;
            Di++;
        }

        _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, 1);

        Aexp[Alen] = exp;

        if (small)
        {
            acc_sm[0] = acc_sm[1] = acc_sm[2] = 0;
            if (Di < Dlen && Dexp[Di] == exp)
            {
                _fmpz_mpoly_add_uiuiui_fmpz(acc_sm, Dcoeff + Di);
                Di++;
            }

            do
            {
                x = _mpoly_heap_pop1(heap, &heap_len, maskhi);

                hind[x->i] |= WORD(1);
                *store++ = x->i;
                *store++ = x->j;
                FLINT_ASSERT(startidx <= x->j);
                FLINT_ASSERT(x->j < ends[0]);
                FLINT_ASSERT(!COEFF_IS_MPZ(Bcoeff[x->i]));
                FLINT_ASSERT(!COEFF_IS_MPZ(Ccoeff[x->j]));
                smul_ppmm(pp1, pp0, Bcoeff[x->i], Ccoeff[x->j]);
                sub_dddmmmsss(acc_sm[2], acc_sm[1], acc_sm[0],
                              acc_sm[2], acc_sm[1], acc_sm[0],
                              FLINT_SIGN_EXT(pp1), pp1, pp0);

                while ((x = x->next) != NULL)
                {
                    hind[x->i] |= WORD(1);
                    *store++ = x->i;
                    *store++ = x->j;
                    FLINT_ASSERT(startidx <= x->j);
                    FLINT_ASSERT(x->j < ends[0]);
                    FLINT_ASSERT(!COEFF_IS_MPZ(Bcoeff[x->i]));
                    FLINT_ASSERT(!COEFF_IS_MPZ(Ccoeff[x->j]));
                    smul_ppmm(pp1, pp0, Bcoeff[x->i], Ccoeff[x->j]);
                    sub_dddmmmsss(acc_sm[2], acc_sm[1], acc_sm[0],
                                  acc_sm[2], acc_sm[1], acc_sm[0],
                                  FLINT_SIGN_EXT(pp1), pp1, pp0);
                }
            } while (heap_len > 1 && heap[1].exp == exp);

            fmpz_set_signed_uiuiui(Acoeff + Alen,
                                             acc_sm[2], acc_sm[1], acc_sm[0]);
        }
        else
        {
            if (Di < Dlen && Dexp[Di] == exp)
            {
                fmpz_set(Acoeff + Alen, Dcoeff + Di);
                Di++;
            }
            else
            {
                fmpz_zero(Acoeff + Alen);
            }

            do
            {
                x = _mpoly_heap_pop1(heap, &heap_len, maskhi);

                hind[x->i] |= WORD(1);
                *store++ = x->i;
                *store++ = x->j;
                fmpz_submul(Acoeff + Alen, Bcoeff + x->i, Ccoeff + x->j);

                while ((x = x->next) != NULL)
                {
                    hind[x->i] |= WORD(1);
                    *store++ = x->i;
                    *store++ = x->j;
                    fmpz_submul(Acoeff + Alen, Bcoeff + x->i, Ccoeff + x->j);
                }
            } while (heap_len > 1 && heap[1].exp == exp);
        }

        Alen += !fmpz_is_zero(Acoeff + Alen);

        /* process nodes taken from the heap */
        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            /* should we go right? */
            if (  (i + 1 < Blen)
               && (j + 0 < ends[i + 1])
               && (hind[i + 1] == 2*j + 1)
               )
            {
                x = chain + i + 1;
                x->i = i + 1;
                x->j = j;
                x->next = NULL;

                hind[x->i] = 2*(x->j + 1) + 0;
                _mpoly_heap_insert1(heap, Bexp[x->i] + Cexp[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
            }

            /* should we go up? */
            if (  (j + 1 < ends[i + 0])
               && ((hind[i] & 1) == 1)
               && (  (i == 0)
                  || (hind[i - 1] >= 2*(j + 2) + 1)
                  )
               )
            {
                x = chain + i;
                x->i = i;
                x->j = j + 1;
                x->next = NULL;

                hind[x->i] = 2*(x->j + 1) + 0;
                _mpoly_heap_insert1(heap, Bexp[x->i] + Cexp[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
            }
        }
    }

    FLINT_ASSERT(Di <= Dlen);
    _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + Dlen - Di, 1);
    if (saveD)
        _fmpz_vec_set(Acoeff + Alen, Dcoeff + Di, Dlen - Di);
    else
        _fmpz_vec_swap(Acoeff + Alen, (fmpz *)(Dcoeff + Di), Dlen - Di);
    mpoly_copy_monomials(Aexp + 1*Alen, Dexp + 1*Di, Dlen - Di, 1);
    Alen += Dlen - Di;

    *A_coeff = Acoeff;
    *A_exp = Aexp;
    *A_alloc = Aalloc;

    return Alen;
}


slong _fmpz_mpoly_mulsub_stripe(fmpz ** A_coeff, ulong ** A_exp, slong * A_alloc,
                         const fmpz * Dcoeff, const ulong * Dexp, slong Dlen, int saveD,
                         const fmpz * Bcoeff, const ulong * Bexp, slong Blen,
                         const fmpz * Ccoeff, const ulong * Cexp, slong Clen,
                                                   const fmpz_mpoly_stripe_t S)
{
    int upperclosed;
    slong startidx, endidx;
    ulong prev_startidx;
    ulong * emax = S->emax;
    ulong * emin = S->emin;
    slong N = S->N;
    slong i, j;
    slong next_loc = Blen + 4;   /* something bigger than heap can ever be */
    slong heap_len = 1; /* heap zero index unused */
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    slong Di;
    slong Alen;
    slong Aalloc = *A_alloc;
    fmpz * Acoeff = *A_coeff;
    ulong * Aexp = *A_exp;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    slong * ends;
    ulong * texp;
    slong * hind;
    int small;
    ulong acc_sm[3], pp0, pp1;

    i = 0;
    hind = (slong *)(S->big_mem + i);
    i += Blen*sizeof(slong);
    ends = (slong *)(S->big_mem + i);
    i += Blen*sizeof(slong);
    store = store_base = (slong *) (S->big_mem + i);
    i += 2*Blen*sizeof(slong);
    heap = (mpoly_heap_s *)(S->big_mem + i);
    i += (Blen + 1)*sizeof(mpoly_heap_s);
    chain = (mpoly_heap_t *)(S->big_mem + i);
    i += Blen*sizeof(mpoly_heap_t);
    exps = (ulong *)(S->big_mem + i);
    i +=  Blen*N*sizeof(ulong);
    exp_list = (ulong **)(S->big_mem + i);
    i +=  Blen*sizeof(ulong *);
    texp = (ulong *)(S->big_mem + i);
    i +=  N*sizeof(ulong);
    FLINT_ASSERT(i <= S->big_mem_alloc);

    exp_next = 0;

    startidx = *S->startidx;
    endidx = *S->endidx;
    upperclosed = S->upperclosed;

    for (i = 0; i < Blen; i++)
        exp_list[i] = exps + i*N;

    /* put all the starting nodes on the heap */
    prev_startidx = -UWORD(1);
    for (i = 0; i < Blen; i++)
    {
        if (startidx < Clen)
        {
            mpoly_monomial_add_mp(texp, Bexp + N*i, Cexp + N*startidx, N);
            FLINT_ASSERT(mpoly_monomial_cmp(emax, texp, N, S->cmpmask)
                                                               > -upperclosed);
        }
        while (startidx > 0)
        {
            mpoly_monomial_add_mp(texp, Bexp + N*i, Cexp + N*(startidx - 1), N);
            if (mpoly_monomial_cmp(emax, texp, N, S->cmpmask) <= -upperclosed)
            {
                break;
            }
            startidx--;
        }

        if (endidx < Clen)
        {
            mpoly_monomial_add_mp(texp, Bexp + N*i, Cexp + N*endidx, N);
            FLINT_ASSERT(mpoly_monomial_cmp(emin, texp, N, S->cmpmask) > 0);
        }
        while (endidx > 0)
        {
            mpoly_monomial_add_mp(texp, Bexp + N*i, Cexp + N*(endidx - 1), N);
            if (mpoly_monomial_cmp(emin, texp, N, S->cmpmask) <= 0)
            {
                break;
            }
            endidx--;
        }

        ends[i] = endidx;

        hind[i] = 2*startidx + 1;

        if (  (startidx < endidx)
           && (((ulong)startidx) < prev_startidx) 
           )
        {
            x = chain + i;
            x->i = i;
            x->j = startidx;
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;

            mpoly_monomial_add_mp(exp_list[exp_next], Bexp + N*x->i,
                                                      Cexp + N*x->j, N);

            exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                          &next_loc, &heap_len, N, S->cmpmask);
        }

        prev_startidx = startidx;
    }

    *S->startidx = startidx;
    *S->endidx = endidx;

    /* whether input coeffs are small, thus output coeffs fit in three words */
    FLINT_ASSERT(ends[0] >= startidx);
    FLINT_ASSERT(S->coeff_bits == FLINT_ABS(_fmpz_vec_max_bits(Ccoeff, Clen)));
    small = S->coeff_bits <= SMALL_FMPZ_BITCOUNT_MAX
            && _fmpz_mpoly_fits_small(Bcoeff, Blen)
            && FLINT_ABS(_fmpz_vec_max_bits(Dcoeff, Dlen)) < 3*FLINT_BITS - 3;

    Alen = 0;
    Di = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;

        while (Di < Dlen && mpoly_monomial_gt(Dexp + N*Di, exp, N, S->cmpmask))
        {
            _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, N);
            mpoly_monomial_set(Aexp + N*Alen, Dexp + N*Di, N);
            if (saveD)
                fmpz_set(Acoeff + Alen, Dcoeff + Di);
            else
                fmpz_swap(Acoeff + Alen, (fmpz *)(Dcoeff + Di));
            Alen++;
            Di++;
        }

        _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, N);

        mpoly_monomial_set(Aexp + N*Alen, exp, N);

        if (small)
        {
            acc_sm[0] = acc_sm[1] = acc_sm[2] = 0;
            if (Di < Dlen && mpoly_monomial_equal(Dexp + N*Di, exp, N))
            {
                _fmpz_mpoly_add_uiuiui_fmpz(acc_sm, Dcoeff + Di);
                Di++;
            }

            do
            {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, S->cmpmask);

                hind[x->i] |= WORD(1);
                *store++ = x->i;
                *store++ = x->j;
                FLINT_ASSERT(startidx <= x->j);
                FLINT_ASSERT(x->j < ends[0]);
                FLINT_ASSERT(!COEFF_IS_MPZ(Bcoeff[x->i]));
                FLINT_ASSERT(!COEFF_IS_MPZ(Ccoeff[x->j]));
                smul_ppmm(pp1, pp0, Bcoeff[x->i], Ccoeff[x->j]);
                sub_dddmmmsss(acc_sm[2], acc_sm[1], acc_sm[0],
                              acc_sm[2], acc_sm[1], acc_sm[0],
                              FLINT_SIGN_EXT(pp1), pp1, pp0);

                while ((x = x->next) != NULL)
                {
                    hind[x->i] |= WORD(1);
                    *store++ = x->i;
                    *store++ = x->j;
                    FLINT_ASSERT(startidx <= x->j);
                    FLINT_ASSERT(x->j < ends[0]);
                    FLINT_ASSERT(!COEFF_IS_MPZ(Bcoeff[x->i]));
                    FLINT_ASSERT(!COEFF_IS_MPZ(Ccoeff[x->j]));
                    smul_ppmm(pp1, pp0, Bcoeff[x->i], Ccoeff[x->j]);
                    sub_dddmmmsss(acc_sm[2], acc_sm[1], acc_sm[0],
                                  acc_sm[2], acc_sm[1], acc_sm[0],
                                  FLINT_SIGN_EXT(pp1), pp1, pp0);
                }
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

            fmpz_set_signed_uiuiui(Acoeff + Alen,
                                             acc_sm[2], acc_sm[1], acc_sm[0]);
        }
        else
        {
            if (Di < Dlen && mpoly_monomial_equal(Dexp + N*Di, exp, N))
            {
                fmpz_set(Acoeff + Alen, Dcoeff + Di);
                Di++;
            }
            else
            {
                fmpz_zero(Acoeff + Alen);
            }

            do
            {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, S->cmpmask);

                hind[x->i] |= WORD(1);
                *store++ = x->i;
                *store++ = x->j;
                FLINT_ASSERT(startidx <= x->j);
                FLINT_ASSERT(x->j < ends[0]);
                fmpz_submul(Acoeff + Alen, Bcoeff + x->i, Ccoeff + x->j);

                while ((x = x->next) != NULL)
                {
                    hind[x->i] |= WORD(1);
                    *store++ = x->i;
                    *store++ = x->j;
                    FLINT_ASSERT(startidx <= x->j);
                    FLINT_ASSERT(x->j < ends[0]);
                    fmpz_submul(Acoeff + Alen, Bcoeff + x->i, Ccoeff + x->j);
                }
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));
        }

        Alen += !fmpz_is_zero(Acoeff + Alen);

        /* process nodes taken from the heap */
        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            /* should we go right? */
            if (  (i + 1 < Blen)
               && (j + 0 < ends[i + 1])
               && (hind[i + 1] == 2*j + 1)
               )
            {
                x = chain + i + 1;
                x->i = i + 1;
                x->j = j;
                x->next = NULL;

                hind[x->i] = 2*(x->j + 1) + 0;

                mpoly_monomial_add_mp(exp_list[exp_next], Bexp + N*x->i,
                                                          Cexp + N*x->j, N);

                exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                          &next_loc, &heap_len, N, S->cmpmask);
            }

            /* should we go up? */
            if (  (j + 1 < ends[i + 0])
               && ((hind[i] & 1) == 1)
               && (  (i == 0)
                  || (hind[i - 1] >= 2*(j + 2) + 1)
                  )
               )
            {
                x = chain + i;
                x->i = i;
                x->j = j + 1;
                x->next = NULL;

                hind[x->i] = 2*(x->j + 1) + 0;

                mpoly_monomial_add_mp(exp_list[exp_next], Bexp + N*x->i,
                                                          Cexp + N*x->j, N);

                exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                          &next_loc, &heap_len, N, S->cmpmask);
            }
        }
    }

    FLINT_ASSERT(Di <= Dlen);
    _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + Dlen - Di, N);
    if (saveD)
        _fmpz_vec_set(Acoeff + Alen, Dcoeff + Di, Dlen - Di);
    else
        _fmpz_vec_swap(Acoeff + Alen, (fmpz *)(Dcoeff + Di), Dlen - Di);
    mpoly_copy_monomials(Aexp + N*Alen, Dexp + N*Di, Dlen - Di, N);
    Alen += Dlen - Di;

    *A_coeff = Acoeff;
    *A_exp = Aexp;
    *A_alloc = Aalloc;

    return Alen;
}

/*
    Q = stripe of A/B (assume A != 0)
    return Qlen = 0 if exact division is impossible
*/
slong _fmpz_mpoly_divides_stripe1(
                        fmpz ** Q_coeff,     ulong ** Q_exp, slong * Q_alloc,
                    const fmpz * Acoeff, const ulong * Aexp, slong Alen,
                    const fmpz * Bcoeff, const ulong * Bexp, slong Blen,
                                                   const fmpz_mpoly_stripe_t S)
{
    flint_bitcnt_t bits = S->bits;
    ulong emin = S->emin[0];
    ulong cmpmask = S->cmpmask[0];
    ulong texp;
    int lt_divides;
    slong i, j, s;
    slong next_loc, heap_len;
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    slong Qlen;
    slong Qalloc = * Q_alloc;
    fmpz * Qcoeff = * Q_coeff;
    ulong * Qexp = * Q_exp;
    ulong exp;
    ulong mask;
    slong * hind;
    fmpz_t acc_lg, r;
    ulong acc_sm[3];
    slong Acoeffbits;
    ulong lc_norm = 0, lc_abs = 0, lc_sign = 0, lc_n = 0, lc_i = 0;
    int small;

    FLINT_ASSERT(Alen > 0);
    FLINT_ASSERT(Blen > 0);
    FLINT_ASSERT(S->N == 1);

    Acoeffbits = _fmpz_vec_max_bits(Acoeff, Alen);
    FLINT_ASSERT(S->coeff_bits == FLINT_ABS(_fmpz_vec_max_bits(Bcoeff, Blen)));

    /* whether intermediate computations A - Q*B will fit in three words */
    /* allow one bit for sign, one bit for subtraction NOT QUITE */
    small = S->coeff_bits <= SMALL_FMPZ_BITCOUNT_MAX
         && FLINT_ABS(Acoeffbits) <= (S->coeff_bits
                                     + FLINT_BIT_COUNT(Blen) + SMALL_FMPZ_BITCOUNT_MAX);

    next_loc = Blen + 4;   /* something bigger than heap can ever be */

    i = 0;
    hind = (slong *)(S->big_mem + i);
    i += Blen*sizeof(slong);
    store = store_base = (slong *) (S->big_mem + i);
    i += 2*Blen*sizeof(slong);
    heap = (mpoly_heap1_s *)(S->big_mem + i);
    i += (Blen + 1)*sizeof(mpoly_heap1_s);
    chain = (mpoly_heap_t *)(S->big_mem + i);
    i += Blen*sizeof(mpoly_heap_t);
    FLINT_ASSERT(i <= S->big_mem_alloc);

    fmpz_init(acc_lg);
    fmpz_init(r);

    for (i = 0; i < Blen; i++)
        hind[i] = 1;

    mask = mpoly_overflow_mask_sp(bits);

    Qlen = WORD(0);

    /* s is the number of terms * (latest quotient) we should put into heap */
    s = Blen;

    /* insert (-1, 0, exp2[0]) into heap */
    heap_len = 2;
    x = chain + 0;
    x->i = -WORD(1);
    x->j = 0;
    x->next = NULL;
    heap[1].next = x;
    HEAP_ASSIGN(heap[1], Aexp[0], x);

    FLINT_ASSERT(mpoly_monomial_cmp1(Aexp[0], emin, cmpmask) >= 0);

    if (small)
    {
        FLINT_ASSERT(!COEFF_IS_MPZ(Bcoeff[0]));
        lc_abs = FLINT_ABS(Bcoeff[0]);
        lc_sign = FLINT_SIGN_EXT(Bcoeff[0]);
        count_leading_zeros(lc_norm, lc_abs);
        lc_n = lc_abs << lc_norm;
        invert_limb(lc_i, lc_n);
    }

    while (heap_len > 1)
    {
        exp = heap[1].exp;

        if (mpoly_monomial_overflows1(exp, mask))
            goto not_exact_division;

        FLINT_ASSERT(mpoly_monomial_cmp1(exp, emin, cmpmask) >= 0);

        _fmpz_mpoly_fit_length(&Qcoeff, &Qexp, &Qalloc, Qlen + 1, 1);

        lt_divides = mpoly_monomial_divides1(Qexp + Qlen, exp, Bexp[0], mask);

        if (small)
        {
            acc_sm[0] = acc_sm[1] = acc_sm[2] = 0;
            do
            {
                x = _mpoly_heap_pop1(heap, &heap_len, cmpmask);
                do
                {
                    *store++ = x->i;
                    *store++ = x->j;
                    if (x->i != -WORD(1))
                        hind[x->i] |= WORD(1);

                    if (x->i == -WORD(1))
                    {
                        _fmpz_mpoly_add_uiuiui_fmpz(acc_sm, Acoeff + x->j);
                    }
                    else
                    {
                        FLINT_ASSERT(!COEFF_IS_MPZ(Bcoeff[x->i]));
                        FLINT_ASSERT(!COEFF_IS_MPZ(Qcoeff[x->j]));
                        _fmpz_mpoly_submul_uiuiui_fmpz(acc_sm,
                                                   Bcoeff[x->i], Qcoeff[x->j]);
                    }
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);
        }
        else
        {
            fmpz_zero(acc_lg);
            do
            {
                x = _mpoly_heap_pop1(heap, &heap_len, cmpmask);
                do
                {
                    *store++ = x->i;
                    *store++ = x->j;
                    if (x->i != -WORD(1))
                        hind[x->i] |= WORD(1);

                    if (x->i == -WORD(1))
                    {
                        fmpz_add(acc_lg, acc_lg, Acoeff + x->j);
                    }
                    else
                    {
                        fmpz_submul(acc_lg, Bcoeff + x->i, Qcoeff + x->j);
                    }
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);
        }

        /* process nodes taken from the heap */
        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            if (i == -WORD(1))
            {
                /* take next dividend term */
                if (j + 1 < Alen)
                {
                    x = chain + 0;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;

                    FLINT_ASSERT(mpoly_monomial_cmp1(Aexp[x->j], emin, cmpmask)
                                                                         >= 0);

                    _mpoly_heap_insert1(heap, Aexp[x->j], x,
                                                &next_loc, &heap_len, cmpmask);
                }
            }
            else
            {
                /* should we go up */
                if (  (i + 1 < Blen)
                   && (hind[i + 1] == 2*j + 1)
                   )
                {
                    x = chain + i + 1;
                    x->i = i + 1;
                    x->j = j;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;

                    texp = Bexp[x->i] + Qexp[x->j];
                    if (mpoly_monomial_cmp1(texp, emin, cmpmask) >= 0)
                    {
                        _mpoly_heap_insert1(heap, texp, x,
                                                &next_loc, &heap_len, cmpmask);
                    }
                    else
                    {
                        hind[x->i] |= 1;
                    }
                }
                /* should we go up? */
                if (j + 1 == Qlen)
                {
                    s++;
                } else if (  ((hind[i] & 1) == 1)
                          && ((i == 1) || (hind[i - 1] >= 2*(j + 2) + 1))
                          )
                {
                    x = chain + i;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;

                    texp = Bexp[x->i] + Qexp[x->j];
                    if (mpoly_monomial_cmp1(texp, emin, cmpmask) >= 0)
                    {
                        _mpoly_heap_insert1(heap, texp, x,
                                                &next_loc, &heap_len, cmpmask);
                    }
                    else
                    {
                        hind[x->i] |= 1;
                    }
                }
            }
        }

        if (small)
        {
            ulong d0, d1, ds = acc_sm[2];

            /* d1:d0 = abs(acc_sm[1:0]) assuming ds is sign extension of acc_sm[1] */
            sub_ddmmss(d1, d0, acc_sm[1]^ds, acc_sm[0]^ds, ds, ds);

            if ((acc_sm[0] | acc_sm[1] | acc_sm[2]) == 0)
            {
                continue;
            }

            if (ds == FLINT_SIGN_EXT(acc_sm[1]) && d1 < lc_abs)
            {
                ulong qq, rr, nhi, nlo;
                FLINT_ASSERT(0 < lc_norm && lc_norm < FLINT_BITS);
                nhi = (d1 << lc_norm) | (d0 >> (FLINT_BITS - lc_norm));
                nlo = d0 << lc_norm;
                udiv_qrnnd_preinv(qq, rr, nhi, nlo, lc_n, lc_i);
                if (rr != UWORD(0))
                    goto not_exact_division;

                if (qq <= COEFF_MAX)
                {
                    _fmpz_demote(Qcoeff + Qlen);
                    Qcoeff[Qlen] = qq;
                    if (ds != lc_sign)
                        Qcoeff[Qlen] = -Qcoeff[Qlen];
                }
                else
                {
                    small = 0;
                    fmpz_set_ui(Qcoeff + Qlen, qq);
                    if (ds != lc_sign)
                        fmpz_neg(Qcoeff + Qlen, Qcoeff + Qlen);
                }
            }
            else
            {
                small = 0;
                fmpz_set_signed_uiuiui(acc_lg, acc_sm[2], acc_sm[1], acc_sm[0]);
                fmpz_fdiv_qr(Qcoeff + Qlen, r, acc_lg, Bcoeff + 0);
                if (!fmpz_is_zero(r))
                    goto not_exact_division;
            }
        }
        else
        {
            if (fmpz_is_zero(acc_lg))
            {
                continue;
            }

            fmpz_fdiv_qr(Qcoeff + Qlen, r, acc_lg, Bcoeff + 0);

            if (!fmpz_is_zero(r))
            {
                goto not_exact_division;
            }
        }

        if (!lt_divides)
        {
            goto not_exact_division;
        }

        if (s > 1)
        {
            i = 1;
            x = chain + i;
            x->i = i;
            x->j = Qlen;
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;

            texp = Bexp[x->i] + Qexp[x->j];
            if (mpoly_monomial_cmp1(texp, emin, cmpmask) >= 0)
            {
                _mpoly_heap_insert1(heap, texp, x,
                                         &next_loc, &heap_len, cmpmask);
            }
            else
            {
                hind[x->i] |= 1;
            }
        }
        s = 1;      
        Qlen++;
    }


cleanup:

    fmpz_clear(acc_lg);
    fmpz_clear(r);

    *Q_alloc = Qalloc;
    *Q_coeff = Qcoeff;
    *Q_exp = Qexp;

    return Qlen;

not_exact_division:
    Qlen = 0;
    goto cleanup;
}

slong _fmpz_mpoly_divides_stripe(
                        fmpz ** Q_coeff,     ulong ** Q_exp, slong * Q_alloc,
                    const fmpz * Acoeff, const ulong * Aexp, slong Alen,
                    const fmpz * Bcoeff, const ulong * Bexp, slong Blen,
                                                   const fmpz_mpoly_stripe_t S)
{
    flint_bitcnt_t bits = S->bits;
    slong N = S->N;
    int lt_divides;
    slong i, j, s;
    slong next_loc, heap_len;
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    slong Qlen;
    slong Qalloc = * Q_alloc;
    fmpz * Qcoeff = * Q_coeff;
    ulong * Qexp = * Q_exp;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    ulong mask;
    slong * hind;
    fmpz_t acc_lg, r;
    ulong acc_sm[3];
    slong Acoeffbits;
    ulong lc_norm = 0, lc_abs = 0, lc_sign = 0, lc_n = 0, lc_i = 0;
    int small;

    FLINT_ASSERT(Alen > 0);
    FLINT_ASSERT(Blen > 0);

    Acoeffbits = _fmpz_vec_max_bits(Acoeff, Alen);
    FLINT_ASSERT(S->coeff_bits == FLINT_ABS(_fmpz_vec_max_bits(Bcoeff, Blen)));

    /* whether intermediate computations A - Q*B will fit in three words */
    /* allow one bit for sign, one bit for subtraction NOT QUITE */
    small = S->coeff_bits <= SMALL_FMPZ_BITCOUNT_MAX
         && FLINT_ABS(Acoeffbits) <= (S->coeff_bits
                                     + FLINT_BIT_COUNT(Blen) + SMALL_FMPZ_BITCOUNT_MAX);

    FLINT_ASSERT(Alen > 0);
    FLINT_ASSERT(Blen > 0);

    next_loc = Blen + 4;   /* something bigger than heap can ever be */

    i = 0;
    hind = (slong *) (S->big_mem + i);
    i += Blen*sizeof(slong);
    store = store_base = (slong *) (S->big_mem + i);
    i += 2*Blen*sizeof(slong);
    heap = (mpoly_heap_s *)(S->big_mem + i);
    i += (Blen + 1)*sizeof(mpoly_heap_s);
    chain = (mpoly_heap_t *)(S->big_mem + i);
    i += Blen*sizeof(mpoly_heap_t);
    exps = (ulong *)(S->big_mem + i);
    i +=  Blen*N*sizeof(ulong);
    exp_list = (ulong **)(S->big_mem + i);
    i +=  Blen*sizeof(ulong *);
    exp = (ulong *)(S->big_mem + i);
    i +=  N*sizeof(ulong);
    FLINT_ASSERT(i <= S->big_mem_alloc);

    fmpz_init(acc_lg);
    fmpz_init(r);

    exp_next = 0;
    for (i = 0; i < Blen; i++)
        exp_list[i] = exps + i*N;

    for (i = 0; i < Blen; i++)
        hind[i] = 1;

    mask = bits <= FLINT_BITS ? mpoly_overflow_mask_sp(bits) : 0;

    Qlen = WORD(0);

    /* s is the number of terms * (latest quotient) we should put into heap */
    s = Blen;

    /* insert (-1, 0, exp2[0]) into heap */
    heap_len = 2;
    x = chain + 0;
    x->i = -WORD(1);
    x->j = 0;
    x->next = NULL;
    heap[1].next = x;
    heap[1].exp = exp_list[exp_next++];

    FLINT_ASSERT(mpoly_monomial_cmp(Aexp + N*0, S->emin, N, S->cmpmask) >= 0);

    mpoly_monomial_set(heap[1].exp, Aexp + N*0, N);

    if (small)
    {
        FLINT_ASSERT(!COEFF_IS_MPZ(Bcoeff[0]));
        lc_abs = FLINT_ABS(Bcoeff[0]);
        lc_sign = FLINT_SIGN_EXT(Bcoeff[0]);
        count_leading_zeros(lc_norm, lc_abs);
        lc_n = lc_abs << lc_norm;
        invert_limb(lc_i, lc_n);
    }

    while (heap_len > 1)
    {
        _fmpz_mpoly_fit_length(&Qcoeff, &Qexp, &Qalloc, Qlen + 1, N);

        mpoly_monomial_set(exp, heap[1].exp, N);

        if (bits <= FLINT_BITS)
        {
            if (mpoly_monomial_overflows(exp, N, mask))
                goto not_exact_division;
            lt_divides = mpoly_monomial_divides(Qexp + N*Qlen, exp,
                                                         Bexp + N*0, N, mask);
        }
        else
        {
            if (mpoly_monomial_overflows_mp(exp, N, bits))
                goto not_exact_division;
            lt_divides = mpoly_monomial_divides_mp(Qexp + N*Qlen, exp,
                                                         Bexp + N*0, N, bits);
        }

        FLINT_ASSERT(mpoly_monomial_cmp(exp, S->emin, N, S->cmpmask) >= 0);

        if (small)
        {
            acc_sm[0] = acc_sm[1] = acc_sm[2] = 0;
            do {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, S->cmpmask);
                do {
                    *store++ = x->i;
                    *store++ = x->j;

                    if (x->i == -WORD(1))
                    {
                        _fmpz_mpoly_add_uiuiui_fmpz(acc_sm, Acoeff + x->j);
                    }
                    else
                    {
                        hind[x->i] |= WORD(1);
                        FLINT_ASSERT(!COEFF_IS_MPZ(Bcoeff[x->i]));
                        FLINT_ASSERT(!COEFF_IS_MPZ(Qcoeff[x->j]));
                        _fmpz_mpoly_submul_uiuiui_fmpz(acc_sm,
                                                   Bcoeff[x->i], Qcoeff[x->j]);
                    }
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));
        }
        else
        {
            fmpz_zero(acc_lg);
            do {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, S->cmpmask);
                do {
                    *store++ = x->i;
                    *store++ = x->j;

                    if (x->i == -WORD(1))
                    {
                        fmpz_add(acc_lg, acc_lg, Acoeff + x->j);
                    }
                    else
                    {
                        hind[x->i] |= WORD(1);
                        fmpz_submul(acc_lg, Bcoeff + x->i, Qcoeff + x->j);
                    }
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));
        }

        /* process nodes taken from the heap */
        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            if (i == -WORD(1))
            {
                /* take next dividend term */
                if (j + 1 < Alen)
                {
                    x = chain + 0;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    mpoly_monomial_set(exp_list[exp_next], Aexp + x->j*N, N);

                    FLINT_ASSERT(mpoly_monomial_cmp(exp_list[exp_next],
                                                 S->emin, N, S->cmpmask) >= 0);

                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                          &next_loc, &heap_len, N, S->cmpmask);
                }
            }
            else
            {
                /* should we go up */
                if (  (i + 1 < Blen)
                   && (hind[i + 1] == 2*j + 1)
                   )
                {
                    x = chain + i + 1;
                    x->i = i + 1;
                    x->j = j;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;

                    mpoly_monomial_add_mp(exp_list[exp_next], Bexp + N*x->i,
                                                              Qexp + N*x->j, N);

                    if (mpoly_monomial_cmp(exp_list[exp_next], S->emin, N,
                                                              S->cmpmask) >= 0)
                    {
                        exp_next += _mpoly_heap_insert(heap, exp_list[exp_next],
                                       x, &next_loc, &heap_len, N, S->cmpmask);
                    }
                    else
                    {
                        hind[x->i] |= 1;
                    }
                }
                /* should we go up? */
                if (j + 1 == Qlen)
                {
                    s++;
                } else if (  ((hind[i] & 1) == 1)
                          && ((i == 1) || (hind[i - 1] >= 2*(j + 2) + 1))
                          )
                {
                    x = chain + i;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;

                    mpoly_monomial_add_mp(exp_list[exp_next], Bexp + N*x->i,
                                                              Qexp + N*x->j, N);

                    if (mpoly_monomial_cmp(exp_list[exp_next], S->emin, N,
                                                              S->cmpmask) >= 0)
                    {
                        exp_next += _mpoly_heap_insert(heap, exp_list[exp_next],
                                       x, &next_loc, &heap_len, N, S->cmpmask);
                    }
                    else
                    {
                        hind[x->i] |= 1;
                    }
                }
            }
        }

        if (small)
        {
            ulong d0, d1, ds = acc_sm[2];

            /* d1:d0 = abs(acc_sm[1:0]) assuming ds is sign extension of acc_sm[1] */
            sub_ddmmss(d1, d0, acc_sm[1]^ds, acc_sm[0]^ds, ds, ds);

            if ((acc_sm[0] | acc_sm[1] | acc_sm[2]) == 0)
            {
                continue;
            }

            if (ds == FLINT_SIGN_EXT(acc_sm[1]) && d1 < lc_abs)
            {
                ulong qq, rr, nhi, nlo;
                FLINT_ASSERT(0 < lc_norm && lc_norm < FLINT_BITS);
                nhi = (d1 << lc_norm) | (d0 >> (FLINT_BITS - lc_norm));
                nlo = d0 << lc_norm;
                udiv_qrnnd_preinv(qq, rr, nhi, nlo, lc_n, lc_i);
                if (rr != UWORD(0))
                    goto not_exact_division;

                if (qq <= COEFF_MAX)
                {
                    _fmpz_demote(Qcoeff + Qlen);
                    Qcoeff[Qlen] = qq;
                    if (ds != lc_sign)
                        Qcoeff[Qlen] = -Qcoeff[Qlen];
                }
                else
                {
                    small = 0;
                    fmpz_set_ui(Qcoeff + Qlen, qq);
                    if (ds != lc_sign)
                        fmpz_neg(Qcoeff + Qlen, Qcoeff + Qlen);
                }
            }
            else
            {
                small = 0;
                fmpz_set_signed_uiuiui(acc_lg, acc_sm[2], acc_sm[1], acc_sm[0]);
                fmpz_fdiv_qr(Qcoeff + Qlen, r, acc_lg, Bcoeff + 0);
                if (!fmpz_is_zero(r))
                    goto not_exact_division;
            }
        }
        else
        {
            if (fmpz_is_zero(acc_lg))
            {
                continue;
            }

            fmpz_fdiv_qr(Qcoeff + Qlen, r, acc_lg, Bcoeff + 0);

            if (!fmpz_is_zero(r))
            {
                goto not_exact_division;
            }
        }

        if (!lt_divides)
        {
            goto not_exact_division;
        }

        if (s > 1)
        {
            i = 1;
            x = chain + i;
            x->i = i;
            x->j = Qlen;
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;

            mpoly_monomial_add_mp(exp_list[exp_next], Bexp + N*x->i,
                                                             Qexp + N*x->j, N);

            if (mpoly_monomial_cmp(exp_list[exp_next], S->emin, N, S->cmpmask)
                                                                          >= 0)
            {
                exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                          &next_loc, &heap_len, N, S->cmpmask);
            }
            else
            {
                hind[x->i] |= 1;                
            }
        }
        s = 1;      
        Qlen++;
    }


cleanup:

    fmpz_clear(acc_lg);
    fmpz_clear(r);

    *Q_alloc = Qalloc;
    *Q_coeff = Qcoeff;
    *Q_exp = Qexp;

    return Qlen;

not_exact_division:
    Qlen = 0;
    goto cleanup;
}


static slong chunk_find_exp(ulong * exp, slong a, const divides_heap_base_t H)
{
    slong N = H->N;
    slong b = H->polyA->length;
    const ulong * Aexp = H->polyA->exps;

try_again:
    FLINT_ASSERT(b >= a);

    FLINT_ASSERT(a > 0);
    FLINT_ASSERT(mpoly_monomial_cmp(Aexp + N*(a - 1), exp, N, H->cmpmask) >= 0);
    FLINT_ASSERT(b >= H->polyA->length
                  ||  mpoly_monomial_cmp(Aexp + N*b, exp, N, H->cmpmask) < 0);

    if (b - a < 5)
    {
        slong i = a;
        while (i < b
                && mpoly_monomial_cmp(Aexp + N*i, exp, N, H->cmpmask) >= 0)
        {
            i++;
        }
        return i;
    }
    else
    {
        slong c = a + (b - a)/2;
        if (mpoly_monomial_cmp(Aexp + N*c, exp, N, H->cmpmask) < 0)
        {
            b = c;
        }
        else
        {
            a = c;
        }
        goto try_again;
    }
}

static void stripe_fit_length(fmpz_mpoly_stripe_struct * S, slong new_len)
{
    slong N = S->N;
    slong new_alloc;
    new_alloc = 0;
    if (N == 1)
    {
        new_alloc += new_len*sizeof(slong);
        new_alloc += new_len*sizeof(slong);
        new_alloc += 2*new_len*sizeof(slong);
        new_alloc += (new_len + 1)*sizeof(mpoly_heap1_s);
        new_alloc += new_len*sizeof(mpoly_heap_t);
    }
    else
    {
        new_alloc += new_len*sizeof(slong);
        new_alloc += new_len*sizeof(slong);
        new_alloc += 2*new_len*sizeof(slong);
        new_alloc += (new_len + 1)*sizeof(mpoly_heap_s);
        new_alloc += new_len*sizeof(mpoly_heap_t);
        new_alloc += new_len*N*sizeof(ulong);
        new_alloc += new_len*sizeof(ulong *);
        new_alloc += N*sizeof(ulong);
    }

    if (S->big_mem_alloc >= new_alloc)
    {
        return;
    }

    new_alloc = FLINT_MAX(new_alloc, S->big_mem_alloc + S->big_mem_alloc/4);
    S->big_mem_alloc = new_alloc;

    if (S->big_mem != NULL)
    {
        S->big_mem = (char *) flint_realloc(S->big_mem, new_alloc);
    }
    else
    {
        S->big_mem = (char *) flint_malloc(new_alloc);
    }

}


static void chunk_mulsub(worker_arg_t W, divides_heap_chunk_t L, slong q_prev_length)
{
    divides_heap_base_struct * H = W->H;
    slong N = H->N;
    fmpz_mpoly_struct * C = L->polyC;
    const fmpz_mpoly_struct * B = H->polyB;
    const fmpz_mpoly_struct * A = H->polyA;
    fmpz_mpoly_ts_struct * Q = H->polyQ;
    fmpz_mpoly_struct * T1 = W->polyT1;
    fmpz_mpoly_stripe_struct * S = W->S;

    S->startidx = &L->startidx;
    S->endidx = &L->endidx;
    S->emin = L->emin;
    S->emax = L->emax;
    S->upperclosed = L->upperclosed;
    FLINT_ASSERT(S->N == N);
    stripe_fit_length(S, q_prev_length - L->mq);

    if (L->Cinited)
    {
        if (N == 1)
        {
            T1->length = _fmpz_mpoly_mulsub_stripe1(
                    &T1->coeffs, &T1->exps, &T1->alloc,
                    C->coeffs, C->exps, C->length, 1,
                    Q->coeffs + L->mq, Q->exps + N*L->mq, q_prev_length - L->mq,
                    B->coeffs, B->exps, B->length, S);
        }
        else
        {
            T1->length = _fmpz_mpoly_mulsub_stripe(
                    &T1->coeffs, &T1->exps, &T1->alloc,
                    C->coeffs, C->exps, C->length, 1,
                    Q->coeffs + L->mq, Q->exps + N*L->mq, q_prev_length - L->mq,
                    B->coeffs, B->exps, B->length, S);
        }
        fmpz_mpoly_swap(C, T1, H->ctx);
    }
    else
    {
        slong startidx, stopidx;
        if (L->upperclosed)
        {
            startidx = 0;
            stopidx = chunk_find_exp(L->emin, 1, H);
        }
        else
        {
            startidx = chunk_find_exp(L->emax, 1, H);
            stopidx = chunk_find_exp(L->emin, startidx, H);
        }

        L->Cinited = 1;
        fmpz_mpoly_init2(C, 16 + stopidx - startidx, H->ctx); /*any is OK*/
        fmpz_mpoly_fit_bits(C, H->bits, H->ctx);
        C->bits = H->bits;

        if (N == 1)
        {
            C->length = _fmpz_mpoly_mulsub_stripe1(
                    &C->coeffs, &C->exps, &C->alloc,
                    A->coeffs + startidx, A->exps + N*startidx, stopidx - startidx, 1,
                    Q->coeffs + L->mq, Q->exps + N*L->mq, q_prev_length - L->mq,
                    B->coeffs, B->exps, B->length, S);
        }
        else
        {
            C->length = _fmpz_mpoly_mulsub_stripe(
                    &C->coeffs, &C->exps, &C->alloc,
                    A->coeffs + startidx, A->exps + N*startidx, stopidx - startidx, 1,
                    Q->coeffs + L->mq, Q->exps + N*L->mq, q_prev_length - L->mq,
                    B->coeffs, B->exps, B->length, S);
        }
    }

    L->mq = q_prev_length;
}

static void trychunk(worker_arg_t W, divides_heap_chunk_t L)
{
    divides_heap_base_struct * H = W->H;
    slong i;
    slong N = H->N;
    fmpz_mpoly_struct * C = L->polyC;
    slong q_prev_length;
    ulong mask;
    const fmpz_mpoly_struct * B = H->polyB;
    const fmpz_mpoly_struct * A = H->polyA;
    fmpz_mpoly_ts_struct * Q = H->polyQ;
    fmpz_mpoly_struct * T2 = W->polyT2;

    mask = 0;
    for (i = 0; i < FLINT_BITS/H->bits; i++)
        mask = (mask << H->bits) + (UWORD(1) << (H->bits - 1));

    /* return if this section has already finished processing */
    if (L->mq < 0)
    {
        return;
    }

    /* process more quotient terms if available */
    q_prev_length = Q->length;
    if (q_prev_length > L->mq)
    {
        if (L->producer == 0 && q_prev_length - L->mq < 20)
            return;

        chunk_mulsub(W, L, q_prev_length);
    }

    if (L->producer == 1)
    {
        divides_heap_chunk_struct * next;
        fmpz * Rcoeff;
        ulong * Rexp;
        slong Rlen;

        /* process the remaining quotient terms */
        q_prev_length = Q->length;
        if (q_prev_length > L->mq)
        {
            chunk_mulsub(W, L, q_prev_length);
        }

        /* find location of remaining terms */
        if (L->Cinited)
        {
            Rlen = C->length;
            Rexp = C->exps;
            Rcoeff = C->coeffs;
        }
        else
        {
            slong startidx, stopidx;
            if (L->upperclosed)
            {
                startidx = 0;
                stopidx = chunk_find_exp(L->emin, 1, H);
            }
            else
            {
                startidx = chunk_find_exp(L->emax, 1, H);
                stopidx = chunk_find_exp(L->emin, startidx, H);
            }
            Rlen = stopidx - startidx;
            Rcoeff = A->coeffs + startidx;
            Rexp = A->exps + N*startidx;
        }

        /* if we have remaining terms, add to quotient  */
        if (Rlen > 0)
        {
            fmpz_mpoly_stripe_struct * S = W->S;
            S->startidx = &L->startidx;
            S->endidx = &L->endidx;
            S->emin = L->emin;
            S->emax = L->emax;
            S->upperclosed = L->upperclosed;
            if (N == 1)
            {
                T2->length = _fmpz_mpoly_divides_stripe1(
                                    &T2->coeffs, &T2->exps, &T2->alloc,
                                       Rcoeff, Rexp, Rlen,
                                       B->coeffs, B->exps, B->length,  S);
            }
            else
            {
                T2->length = _fmpz_mpoly_divides_stripe(
                                    &T2->coeffs, &T2->exps, &T2->alloc,
                                       Rcoeff, Rexp, Rlen,
                                       B->coeffs, B->exps, B->length,  S);
            }
            if (T2->length == 0)
            {
                H->failed = 1;
                return;
            }
            else
            {
                fmpz_mpoly_ts_append(H->polyQ, T2->coeffs, T2->exps,
                                                                T2->length, N);
            }
        }

        next = L->next;
        H->length--;
        H->cur = next;

        if (next != NULL)
        {
            next->producer = 1;
        }

        L->producer = 0;
        L->mq = -1;
    }

    return;
}


static void worker_loop(void * varg)
{
    worker_arg_struct * W = (worker_arg_struct *) varg;
    divides_heap_base_struct * H = W->H;
    fmpz_mpoly_stripe_struct * S = W->S;
    const fmpz_mpoly_struct * B = H->polyB;
    fmpz_mpoly_struct * T1 = W->polyT1;
    fmpz_mpoly_struct * T2 = W->polyT2;
    slong N = H->N;
    slong Blen = B->length;

    /* initialize stripe working memory */
    S->N = N;
    S->bits = H->bits;
    S->coeff_bits = FLINT_ABS(H->polyBcoeff_bits);
    S->cmpmask = H->cmpmask;
    S->big_mem_alloc = 0;
    S->big_mem = NULL;

    stripe_fit_length(S, Blen);

    fmpz_mpoly_init2(T1, 16, H->ctx);
    fmpz_mpoly_fit_bits(T1, H->bits, H->ctx);
    T1->bits = H->bits;
    fmpz_mpoly_init2(T2, 16, H->ctx);
    fmpz_mpoly_fit_bits(T2, H->bits, H->ctx);
    T2->bits = H->bits;

    while (!H->failed)
    {
        divides_heap_chunk_struct * L;
        L = H->cur;

        if (L == NULL)
        {
            break;
        }
        while (L != NULL)
        {
#if FLINT_USES_PTHREAD
            pthread_mutex_lock(&H->mutex);
#endif
	    if (L->lock != -1)
            {
                L->lock = -1;
#if FLINT_USES_PTHREAD
                pthread_mutex_unlock(&H->mutex);
#endif
		trychunk(W, L);
#if FLINT_USES_PTHREAD
                pthread_mutex_lock(&H->mutex);
#endif
		L->lock = 0;
#if FLINT_USES_PTHREAD
                pthread_mutex_unlock(&H->mutex);
#endif
		break;
            }
            else
            {
#if FLINT_USES_PTHREAD
                pthread_mutex_unlock(&H->mutex);
#endif
	    }

            L = L->next;
        }
    }

    fmpz_mpoly_clear(T1, H->ctx);
    fmpz_mpoly_clear(T2, H->ctx);
    flint_free(S->big_mem);

    return;
}


/* return 1 if quotient is exact */
int _fmpz_mpoly_divides_heap_threaded_pool(
    fmpz_mpoly_t Q,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles)
{
    ulong mask;
    int divides;
    fmpz_mpoly_ctx_t zctx;
    fmpz_mpoly_t S;
    slong i, k, N;
    flint_bitcnt_t exp_bits;
    ulong * cmpmask;
    ulong * Aexp, * Bexp;
    int freeAexp, freeBexp;
    worker_arg_struct * worker_args;
    fmpz_t qcoeff, r;
    ulong * texps, * qexps;
    divides_heap_base_t H;
    TMP_INIT;

#if !FLINT_KNOW_STRONG_ORDER
    return fmpz_mpoly_divides_monagan_pearce(Q, A, B, ctx);
#endif

    if (B->length < 2 || A->length < 2)
    {
        return fmpz_mpoly_divides_monagan_pearce(Q, A, B, ctx);
    }

    TMP_START;

    fmpz_init(qcoeff);
    fmpz_init(r);

    exp_bits = MPOLY_MIN_BITS;
    exp_bits = FLINT_MAX(exp_bits, A->bits);
    exp_bits = FLINT_MAX(exp_bits, B->bits);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);

    N = mpoly_words_per_exp(exp_bits, ctx->minfo);
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, exp_bits, ctx->minfo);

    /* ensure input exponents packed to same size as output exponents */
    Aexp = A->exps;
    freeAexp = 0;
    if (exp_bits > A->bits)
    {
        freeAexp = 1;
        Aexp = (ulong *) flint_malloc(N*A->length*sizeof(ulong));
        mpoly_repack_monomials(Aexp, exp_bits, A->exps, A->bits,
                                                        A->length, ctx->minfo);
    }

    Bexp = B->exps;
    freeBexp = 0;
    if (exp_bits > B->bits)
    {
        freeBexp = 1;
        Bexp = (ulong *) flint_malloc(N*B->length*sizeof(ulong));
        mpoly_repack_monomials(Bexp, exp_bits, B->exps, B->bits,
                                                    B->length, ctx->minfo);
    }

    fmpz_mpoly_ctx_init(zctx, ctx->minfo->nvars, ctx->minfo->ord);
    fmpz_mpoly_init(S, zctx);

    fmpz_mod(r, A->coeffs + 0, B->coeffs + 0);
    if (!fmpz_is_zero(r))
    {
        divides = 0;
        fmpz_mpoly_zero(Q, ctx);
        goto cleanup1;
    }

    if (mpoly_divides_select_exps(S, zctx, num_handles,
                                   Aexp, A->length, Bexp, B->length, exp_bits))
    {
        divides = 0;
        fmpz_mpoly_zero(Q, ctx);
        goto cleanup1;
    }

    /*
        At this point A and B both have at least two terms
            and the leading coefficients and monomials divide
            and the exponent selection did not give an easy exit
    */

    divides_heap_base_init(H);

    H->polyA->coeffs = A->coeffs;
    H->polyA->exps = Aexp;
    H->polyA->bits = exp_bits;
    H->polyA->length = A->length;
    H->polyA->alloc = A->alloc;

    H->polyB->coeffs = B->coeffs;
    H->polyB->exps = Bexp;
    H->polyB->bits = exp_bits;
    H->polyB->length = B->length;
    H->polyB->alloc = B->alloc;

    H->polyBcoeff_bits = _fmpz_vec_max_bits(B->coeffs, B->length);

    H->ctx = ctx;
    H->bits = exp_bits;
    H->N = N;
    H->cmpmask = cmpmask;
    H->failed = 0;

    for (i = 0; i + 1 < S->length; i++)
    {
        divides_heap_chunk_struct * L;
        L = (divides_heap_chunk_struct *) flint_malloc(
                                            sizeof(divides_heap_chunk_struct));
        L->ma = 0;
        L->mq = 0;
        L->emax = S->exps + N*i;
        L->emin = S->exps + N*(i + 1);
        L->upperclosed = 0;
        L->startidx = B->length;
        L->endidx = B->length;
        L->producer = 0;
        L->Cinited = 0;
        L->lock = -2;
        divides_heap_base_add_chunk(H, L);
    }

    H->head->upperclosed = 1;
    H->head->producer = 1;
    H->cur = H->head;

    /* generate at least the first quotient term */

    texps = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    qexps = (ulong *) TMP_ALLOC(N*sizeof(ulong));

    mpoly_monomial_sub_mp(qexps + N*0, Aexp + N*0, Bexp + N*0, N);
    fmpz_divexact(qcoeff, A->coeffs + 0, B->coeffs + 0); /* already checked */

    fmpz_mpoly_ts_init(H->polyQ, qcoeff, qexps, 1, H->bits, H->N);

    mpoly_monomial_add_mp(texps, qexps + N*0, Bexp + N*1, N);

    mask = 0;
    for (i = 0; i < FLINT_BITS/exp_bits; i++)
        mask = (mask << exp_bits) + (UWORD(1) << (exp_bits - 1));

    k = 1;
    while (k < A->length && mpoly_monomial_gt(Aexp + N*k, texps, N, cmpmask))
    {
        int lt_divides;
        if (exp_bits <= FLINT_BITS)
            lt_divides = mpoly_monomial_divides(qexps, Aexp + N*k,
                                                      Bexp + N*0, N, mask);
        else
            lt_divides = mpoly_monomial_divides_mp(qexps, Aexp + N*k,
                                                      Bexp + N*0, N, exp_bits);
        if (!lt_divides)
        {
            H->failed = 1;
            break;
        }

        fmpz_fdiv_qr(qcoeff, r, A->coeffs + k, B->coeffs + 0);
        if (!fmpz_is_zero(r))
        {
            H->failed = 1;
            break;
        }

        fmpz_mpoly_ts_append(H->polyQ, qcoeff, qexps, 1, H->N);
        k++;
    }

    /* start the workers */

#if FLINT_USES_PTHREAD
    pthread_mutex_init(&H->mutex, NULL);
#endif

    worker_args = (worker_arg_struct *) flint_malloc((num_handles + 1)
                                                        *sizeof(worker_arg_t));

    for (i = 0; i < num_handles; i++)
    {
        (worker_args + i)->H = H;
        thread_pool_wake(global_thread_pool, handles[i], 0,
                                                 worker_loop, worker_args + i);
    }
    (worker_args + num_handles)->H = H;
    worker_loop(worker_args + num_handles);
    for (i = 0; i < num_handles; i++)
    {
        thread_pool_wait(global_thread_pool, handles[i]);
    }

    flint_free(worker_args);

#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&H->mutex);
#endif

    divides = divides_heap_base_clear(Q, H);

cleanup1:
    fmpz_mpoly_clear(S, zctx);
    fmpz_mpoly_ctx_clear(zctx);

    if (freeAexp)
        flint_free(Aexp);

    if (freeBexp)
        flint_free(Bexp);

    fmpz_clear(qcoeff);
    fmpz_clear(r);

    TMP_END;

    return divides;
}


int fmpz_mpoly_divides_heap_threaded(
    fmpz_mpoly_t Q,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    thread_pool_handle * handles;
    slong num_handles;
    int divides;
    slong thread_limit = A->length/32;

    if (B->length < 2 || A->length < 2)
    {
        if (B->length == 0)
        {
            flint_throw(FLINT_DIVZERO,
                         "Divide by zero in fmpz_mpoly_divides_heap_threaded");
        }

        if (A->length == 0)
        {
            fmpz_mpoly_zero(Q, ctx);
            return 1;
        }
        return fmpz_mpoly_divides_monagan_pearce(Q, A, B, ctx);
    }

    num_handles = flint_request_threads(&handles, thread_limit);

    divides = _fmpz_mpoly_divides_heap_threaded_pool(Q, A, B, ctx,
                                                         handles, num_handles);

    flint_give_back_threads(handles, num_handles);

    return divides;
}

/*
   Set poly1 to poly2/poly3 if the division is exact, and return the length
   of the quotient. Otherwise return 0. This version of the function assumes
   the exponent vectors all fit in a single word. The exponent vectors are
   assumed to have fields with the given number of bits. Assumes input polys
   are nonzero. Implements "Polynomial division using dynamic arrays, heaps
   and packed exponents" by Michael Monagan and Roman Pearce [1], except that
   we divide from right to left and use a heap with smallest exponent at head.
   [1] http://www.cecm.sfu.ca/~rpearcea/sdmp/sdmp_paper.pdf
*/
slong _fmpz_mpoly_divides_monagan_pearce1(fmpz ** poly1, ulong ** exp1,
         slong * alloc, const fmpz * poly2, const ulong * exp2, slong len2,
                const fmpz * poly3, const ulong * exp3, slong len3, slong bits,
                                                                  ulong maskhi)
{
    slong i, j, k, s;
    slong next_loc, heap_len = 2;
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    fmpz * p1 = *poly1;
    ulong * e1 = *exp1;
    slong * hind;
    ulong mask, exp, maxexp = exp2[len2 - 1];
    fmpz_t r, acc_lg;
    ulong acc_sm[3];
    int lt_divides, small;
    slong bits2, bits3;
    ulong lc_norm = 0, lc_abs = 0, lc_sign = 0, lc_n = 0, lc_i = 0;
    TMP_INIT;

    TMP_START;

    fmpz_init(acc_lg);
    fmpz_init(r);

    /* whether intermediate computations q - a*b will fit in three words */
    bits2 = _fmpz_vec_max_bits(poly2, len2);
    bits3 = _fmpz_vec_max_bits(poly3, len3);
    /* allow one bit for sign, one bit for subtraction */
    small = FLINT_ABS(bits2) <= (FLINT_ABS(bits3) + FLINT_BIT_COUNT(len3) + SMALL_FMPZ_BITCOUNT_MAX)
         && FLINT_ABS(bits3) <= SMALL_FMPZ_BITCOUNT_MAX;

    /* alloc array of heap nodes which can be chained together */
    next_loc = len3 + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap1_s *) TMP_ALLOC((len3 + 1)*sizeof(mpoly_heap1_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(len3*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*len3*sizeof(slong));

    /* space for flagged heap indices */
    hind = (slong *) TMP_ALLOC(len3*sizeof(slong));
    for (i = 0; i < len3; i++)
        hind[i] = 1;

    mask = mpoly_overflow_mask_sp(bits);

    /* output poly index starts at -1, will be immediately updated to 0 */
    k = -WORD(1);

    /* s is the number of terms * (latest quotient) we should put into heap */
    s = len3;

    /* insert (-1, 0, exp2[0]) into heap */
    x = chain + 0;
    x->i = -WORD(1);
    x->j = 0;
    x->next = NULL;
    HEAP_ASSIGN(heap[1], exp2[0], x);

    if (small)
    {
        lc_abs = FLINT_ABS(poly3[0]);
        lc_sign = FLINT_SIGN_EXT(poly3[0]);
        count_leading_zeros(lc_norm, lc_abs);
        lc_n = lc_abs << lc_norm;
        invert_limb(lc_i, lc_n);
    }

    while (heap_len > 1)
    {
        exp = heap[1].exp;

        if (mpoly_monomial_overflows1(exp, mask))
            goto not_exact_division;

        k++;
        _fmpz_mpoly_fit_length(&p1, &e1, alloc, k + 1, 1);

        lt_divides = mpoly_monomial_divides1(e1 + k, exp, exp3[0], mask);

        if (small)
        {
            acc_sm[0] = acc_sm[1] = acc_sm[2] = 0;
            do
            {
                x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
                do
                {
                    *store++ = x->i;
                    *store++ = x->j;
                    if (x->i != -WORD(1))
                        hind[x->i] |= WORD(1);

                    if (x->i == -WORD(1))
                        _fmpz_mpoly_add_uiuiui_fmpz(acc_sm, poly2 + x->j);
                    else
                        _fmpz_mpoly_submul_uiuiui_fmpz(acc_sm, poly3[x->i], p1[x->j]);
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);
        } else
        {
            fmpz_zero(acc_lg);  
            do
            {
                x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
                do
                {
                    *store++ = x->i;
                    *store++ = x->j;
                    if (x->i != -WORD(1))
                        hind[x->i] |= WORD(1);

                    if (x->i == -WORD(1))
                        fmpz_add(acc_lg, acc_lg, poly2 + x->j);
                    else
                        fmpz_submul(acc_lg, poly3 + x->i, p1 + x->j);
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);
        }

        /* process nodes taken from the heap */
        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            if (i == -WORD(1))
            {
                /* take next dividend term */
                if (j + 1 < len2)
                {
                    x = chain + 0;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    _mpoly_heap_insert1(heap, exp2[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
                }

            } else
            {
                /* should we go right? */
                if (  (i + 1 < len3)
                   && (hind[i + 1] == 2*j + 1)
                   )
                {
                    x = chain + i + 1;
                    x->i = i + 1;
                    x->j = j;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;
                    _mpoly_heap_insert1(heap, exp3[x->i] + e1[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
                }
                /* should we go up? */
                if (j + 1 == k)
                {
                    s++;
                } else if (  ((hind[i] & 1) == 1)
                          && ((i == 1) || (hind[i - 1] >= 2*(j + 2) + 1))
                          )
                {
                    x = chain + i;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;
                    _mpoly_heap_insert1(heap, exp3[x->i] + e1[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
                }
            }
        }

        /* try to divide accumulated term by leading term */
        if (small)
        {
            ulong d0, d1, ds = acc_sm[2];

            /* d1:d0 = abs(acc_sm[1:0]) assuming ds is sign extension of acc_sm[1] */
            sub_ddmmss(d1, d0, acc_sm[1]^ds, acc_sm[0]^ds, ds, ds);
            
            if ((acc_sm[0] | acc_sm[1] | acc_sm[2]) == 0)
            {
                k--;
                continue;
            }

            if (ds == FLINT_SIGN_EXT(acc_sm[1]) && d1 < lc_abs)
            {
                ulong qq, rr, nhi, nlo;
                FLINT_ASSERT(0 < lc_norm && lc_norm < FLINT_BITS);
                nhi = (d1 << lc_norm) | (d0 >> (FLINT_BITS - lc_norm));
                nlo = d0 << lc_norm;
                udiv_qrnnd_preinv(qq, rr, nhi, nlo, lc_n, lc_i);
                if (rr != WORD(0))
                    goto not_exact_division;

                if ((qq & (WORD(3) << (SMALL_FMPZ_BITCOUNT_MAX))) == WORD(0))
                {
                    _fmpz_demote(p1 + k);
                    p1[k] = (qq^ds^lc_sign) - (ds^lc_sign);
                } else
                {
                    small = 0;
                    fmpz_set_ui(p1 + k, qq);
                    if (ds != lc_sign)
                        fmpz_neg(p1 + k, p1 + k);
                }
            } else
            {
                small = 0;
                fmpz_set_signed_uiuiui(acc_lg, acc_sm[2], acc_sm[1], acc_sm[0]);
                fmpz_fdiv_qr(p1 + k, r, acc_lg, poly3 + 0);
                if (!fmpz_is_zero(r))
                    goto not_exact_division;
            }

        } else
        {
            if (fmpz_is_zero(acc_lg))
            {
                k--;
                continue;
            }
            fmpz_fdiv_qr(p1 + k, r, acc_lg, poly3 + 0);
            if (!fmpz_is_zero(r))
                goto not_exact_division;
        }

        if (!lt_divides || (exp^maskhi) < (maxexp^maskhi))
            goto not_exact_division;

        /* put newly generated quotient term back into the heap if necessary */
        if (s > 1)
        {
            i = 1;
            x = chain + i;
            x->i = i;
            x->j = k;
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;
            _mpoly_heap_insert1(heap, exp3[x->i] + e1[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
        }
        s = 1;
    }

    k++;

cleanup:

    fmpz_clear(acc_lg);
    fmpz_clear(r);

    (*poly1) = p1;
    (*exp1) = e1;

    TMP_END;

    return k;

not_exact_division:
    for (i = 0; i <= k; i++)
        _fmpz_demote(p1 + i);
    k = 0;
    goto cleanup;
}


slong _fmpz_mpoly_divides_monagan_pearce(fmpz ** poly1, ulong ** exp1,
         slong * alloc, const fmpz * poly2, const ulong * exp2, slong len2,
       const fmpz * poly3, const ulong * exp3, slong len3, flint_bitcnt_t bits, slong N,
                                                         const ulong * cmpmask)
{
    slong i, j, k, s;
    slong next_loc;
    slong heap_len = 2; /* heap zero index unused */
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    fmpz * p1 = *poly1;
    ulong * e1 = *exp1;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    fmpz_t r, acc_lg;
    ulong acc_sm[3];
    ulong mask;
    slong * hind;
    int lt_divides, small;
    slong bits2, bits3;
    ulong lc_norm = 0, lc_abs = 0, lc_sign = 0, lc_n = 0, lc_i = 0;
    TMP_INIT;

    /* if exponent vectors are all one word, call specialised version */
    if (N == 1)
        return _fmpz_mpoly_divides_monagan_pearce1(poly1, exp1, alloc,
                       poly2, exp2, len2, poly3, exp3, len3, bits, cmpmask[0]);

    TMP_START;

    fmpz_init(acc_lg);
    fmpz_init(r);

    /* whether intermediate computations q - a*b will fit in three words */
    bits2 = _fmpz_vec_max_bits(poly2, len2);
    bits3 = _fmpz_vec_max_bits(poly3, len3);
    /* allow one bit for sign, one bit for subtraction */
    small = FLINT_ABS(bits2) <= (FLINT_ABS(bits3) + FLINT_BIT_COUNT(len3) + SMALL_FMPZ_BITCOUNT_MAX)
          && FLINT_ABS(bits3) <= SMALL_FMPZ_BITCOUNT_MAX;

    /* alloc array of heap nodes which can be chained together */
    next_loc = len3 + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap_s *) TMP_ALLOC((len3 + 1)*sizeof(mpoly_heap_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(len3*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*len3*sizeof(slong));

    /* array of exponent vectors, each of "N" words */
    exps = (ulong *) TMP_ALLOC(len3*N*sizeof(ulong));
    /* list of pointers to available exponent vectors */
    exp_list = (ulong **) TMP_ALLOC(len3*sizeof(ulong *));
    /* space to save copy of current exponent vector */
    exp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    /* set up list of available exponent vectors */
    exp_next = 0;
    for (i = 0; i < len3; i++)
        exp_list[i] = exps + i*N;

    /* space for flagged heap indices */
    hind = (slong *) TMP_ALLOC(len3*sizeof(slong));
    for (i = 0; i < len3; i++)
        hind[i] = 1;

    mask = bits <= FLINT_BITS ? mpoly_overflow_mask_sp(bits) : 0;

    /* output poly index starts at -1, will be immediately updated to 0 */
    k = -WORD(1);

    /* s is the number of terms * (latest quotient) we should put into heap */
    s = len3;

    /* insert (-1, 0, exp2[0]) into heap */
    x = chain + 0;
    x->i = -WORD(1);
    x->j = 0;
    x->next = NULL;
    heap[1].next = x;
    heap[1].exp = exp_list[exp_next++];
    mpoly_monomial_set(heap[1].exp, exp2, N);

    if (small)
    {
        lc_abs = FLINT_ABS(poly3[0]);
        lc_sign = FLINT_SIGN_EXT(poly3[0]);
        count_leading_zeros(lc_norm, lc_abs);
        lc_n = lc_abs << lc_norm;
        invert_limb(lc_i, lc_n);
    }

    while (heap_len > 1)
    {
        mpoly_monomial_set(exp, heap[1].exp, N);

        if (bits <= FLINT_BITS)
        {
            if (mpoly_monomial_overflows(exp, N, mask))
                goto not_exact_division;
        } else
        {
            if (mpoly_monomial_overflows_mp(exp, N, bits))
                goto not_exact_division;
        }
      
        k++;
        _fmpz_mpoly_fit_length(&p1, &e1, alloc, k + 1, N);

        if (bits <= FLINT_BITS)
            lt_divides = mpoly_monomial_divides(e1 + k*N, exp, exp3, N, mask);
        else
            lt_divides = mpoly_monomial_divides_mp(e1 + k*N, exp, exp3, N, bits);

        if (small) 
        {
            acc_sm[0] = acc_sm[1] = acc_sm[2] = 0;
            do
            {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
                do
                {
                    *store++ = x->i;
                    *store++ = x->j;
                    if (x->i != -WORD(1))
                        hind[x->i] |= WORD(1);

                    if (x->i == -WORD(1))
                        _fmpz_mpoly_add_uiuiui_fmpz(acc_sm, poly2 + x->j);
                    else
                        _fmpz_mpoly_submul_uiuiui_fmpz(acc_sm, poly3[x->i], p1[x->j]);
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));
        } else
        {
            fmpz_zero(acc_lg);
            do
            {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
                do
                {
                    *store++ = x->i;
                    *store++ = x->j;
                    if (x->i != -WORD(1))
                        hind[x->i] |= WORD(1);

                    if (x->i == -WORD(1))
                        fmpz_add(acc_lg, acc_lg, poly2 + x->j);
                    else
                        fmpz_submul(acc_lg, poly3 + x->i, p1 + x->j);
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));
        }

        /* process nodes taken from the heap */
        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            if (i == -WORD(1))
            {
                /* take next dividend term */
                if (j + 1 < len2)
                {
                    x = chain + 0;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    mpoly_monomial_set(exp_list[exp_next], exp2 + x->j*N, N);
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
            } else
            {
                /* should we go right? */
                if (  (i + 1 < len3)
                   && (hind[i + 1] == 2*j + 1)
                   )
                {
                    x = chain + i + 1;
                    x->i = i + 1;
                    x->j = j;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;

                    if (bits <= FLINT_BITS)
                        mpoly_monomial_add(exp_list[exp_next], exp3 + x->i*N,
                                                               e1 + x->j*N, N);
                    else
                        mpoly_monomial_add_mp(exp_list[exp_next], exp3 + x->i*N,
                                                               e1 + x->j*N, N);                        

                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
                /* should we go up? */
                if (j + 1 == k)
                {
                    s++;
                } else if (  ((hind[i] & 1) == 1)
                          && ((i == 1) || (hind[i - 1] >= 2*(j + 2) + 1))
                          )
                {
                    x = chain + i;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;

                    if (bits <= FLINT_BITS)
                        mpoly_monomial_add(exp_list[exp_next], exp3 + x->i*N,
                                                               e1 + x->j*N, N);
                    else
                        mpoly_monomial_add_mp(exp_list[exp_next], exp3 + x->i*N,
                                                               e1 + x->j*N, N);

                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
            }
        }

        /* try to divide accumulated term by leading term */
        if (small)
        {
            ulong d0, d1, ds = acc_sm[2];

            /* d1:d0 = abs(acc_sm[1:0]) assuming ds is sign extension of acc_sm[1] */
            sub_ddmmss(d1, d0, acc_sm[1]^ds, acc_sm[0]^ds, ds, ds);
            
            if ((acc_sm[0] | acc_sm[1] | acc_sm[2]) == 0)
            {
                k--;
                continue;
            }

            if (ds == FLINT_SIGN_EXT(acc_sm[1]) && d1 < lc_abs)
            {
                ulong qq, rr, nhi, nlo;
                FLINT_ASSERT(0 < lc_norm && lc_norm < FLINT_BITS);
                nhi = (d1 << lc_norm) | (d0 >> (FLINT_BITS - lc_norm));
                nlo = d0 << lc_norm;
                udiv_qrnnd_preinv(qq, rr, nhi, nlo, lc_n, lc_i);
                if (rr != WORD(0))
                    goto not_exact_division;

                if ((qq & (WORD(3) << (SMALL_FMPZ_BITCOUNT_MAX))) == WORD(0))
                {
                    _fmpz_demote(p1 + k);
                    p1[k] = (qq^(ds^lc_sign)) - (ds^lc_sign);
                }
                else
                {
                    small = 0;
                    fmpz_set_ui(p1 + k, qq);
                    if (ds != lc_sign)
                        fmpz_neg(p1 + k, p1 + k);
                }
            } else
            {
                small = 0;
                fmpz_set_signed_uiuiui(acc_lg, acc_sm[2], acc_sm[1], acc_sm[0]);
                fmpz_fdiv_qr(p1 + k, r, acc_lg, poly3 + 0);
                if (!fmpz_is_zero(r))
                    goto not_exact_division;
            }

        } else
        {
            if (fmpz_is_zero(acc_lg))
            {
                k--;
                continue;
            }
            fmpz_fdiv_qr(p1 + k, r, acc_lg, poly3 + 0);
            if (!fmpz_is_zero(r))
                goto not_exact_division;
        }

        if (!lt_divides ||
                mpoly_monomial_gt(exp2 + (len2 - 1)*N, exp, N, cmpmask))
            goto not_exact_division;

        if (s > 1)
        {
            i = 1;
            x = chain + i;
            x->i = i;
            x->j = k;
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;

            if (bits <= FLINT_BITS)
                mpoly_monomial_add(exp_list[exp_next], exp3 + x->i*N,
                                                               e1 + x->j*N, N);
            else
                mpoly_monomial_add_mp(exp_list[exp_next], exp3 + x->i*N,
                                                               e1 + x->j*N, N);

            exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
        }
        s = 1;      
    }

    k++;

cleanup:

    fmpz_clear(acc_lg);
    fmpz_clear(r);

    (*poly1) = p1;
    (*exp1) = e1;

    TMP_END;

    return k;

not_exact_division:
    for (i = 0; i <= k; i++)
        _fmpz_demote(p1 + i);
    k = 0;
    goto cleanup;
}

/* return 1 if quotient is exact */
int fmpz_mpoly_divides_monagan_pearce(fmpz_mpoly_t poly1,
                  const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i, N, len = 0;
    flint_bitcnt_t exp_bits;
    fmpz * max_fields2, * max_fields3;
    ulong * cmpmask;
    ulong * exp2 = poly2->exps, * exp3 = poly3->exps, * expq;
    int easy_exit, free2 = 0, free3 = 0;
    ulong mask = 0;
    TMP_INIT;

   /* check divisor is nonzero */
   if (poly3->length == 0)
      flint_throw(FLINT_DIVZERO, "Divide by zero in fmpz_mpoly_divides_monagan_pearce");

   /* dividend zero, write out quotient */
   if (poly2->length == 0)
   {
      fmpz_mpoly_zero(poly1, ctx);

      return 1;
   }

   TMP_START;

    max_fields2 = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    max_fields3 = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_init(max_fields2 + i);
        fmpz_init(max_fields3 + i);
    }

    mpoly_max_fields_fmpz(max_fields2, poly2->exps, poly2->length,
                                                      poly2->bits, ctx->minfo);
    mpoly_max_fields_fmpz(max_fields3, poly3->exps, poly3->length,
                                                      poly3->bits, ctx->minfo);

    easy_exit = 0;
    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        /*
            cannot be exact division if any max field from poly2
            is less than corresponding max field from poly3
        */
        if (fmpz_cmp(max_fields2 + i, max_fields3 + i) < 0)
            easy_exit = 1;
    }

    exp_bits = _fmpz_vec_max_bits(max_fields2, ctx->minfo->nfields);
    exp_bits = FLINT_MAX(MPOLY_MIN_BITS, exp_bits + 1);
    exp_bits = FLINT_MAX(exp_bits, poly2->bits);
    exp_bits = FLINT_MAX(exp_bits, poly3->bits);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);

    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_clear(max_fields2 + i);
        fmpz_clear(max_fields3 + i);
    }

    if (easy_exit)
    {
        len = 0;
        goto cleanup;
    }

    N = mpoly_words_per_exp(exp_bits, ctx->minfo);
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, exp_bits, ctx->minfo);

   /* temporary space to check leading monomials divide */
   expq = (ulong *) TMP_ALLOC(N*sizeof(ulong));

   /* quick check for easy case of inexact division of leading monomials */
   if (poly2->bits == poly3->bits && N == 1 && 
       poly2->exps[0] < poly3->exps[0])
   {
      goto cleanup;
   }

   /* ensure input exponents packed to same size as output exponents */
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


    /* check leading monomial divides exactly */
    if (exp_bits <= FLINT_BITS)
    {
        /* mask with high bit of each exponent vector field set */
        for (i = 0; i < FLINT_BITS/exp_bits; i++)
            mask = (mask << exp_bits) + (UWORD(1) << (exp_bits - 1));

        if (!mpoly_monomial_divides(expq, exp2, exp3, N, mask))
        {
            len = 0;
            goto cleanup;
        }
    } else
    {
        if (!mpoly_monomial_divides_mp(expq, exp2, exp3, N, exp_bits))
        {
            len = 0;
            goto cleanup;
        }
    }

   /* deal with aliasing and divide polynomials */
   if (poly1 == poly2 || poly1 == poly3)
   {
      fmpz_mpoly_t temp;

      fmpz_mpoly_init2(temp, poly2->length/poly3->length + 1, ctx);
      fmpz_mpoly_fit_bits(temp, exp_bits, ctx);
      temp->bits = exp_bits;

      len = _fmpz_mpoly_divides_monagan_pearce(&temp->coeffs, &temp->exps,
                            &temp->alloc, poly2->coeffs, exp2, poly2->length,
                              poly3->coeffs, exp3, poly3->length, exp_bits, N,
                                                              cmpmask);

      fmpz_mpoly_swap(temp, poly1, ctx);

      fmpz_mpoly_clear(temp, ctx);
   } else
   {
      fmpz_mpoly_fit_length(poly1, poly2->length/poly3->length + 1, ctx);
      fmpz_mpoly_fit_bits(poly1, exp_bits, ctx);
      poly1->bits = exp_bits;

      len = _fmpz_mpoly_divides_monagan_pearce(&poly1->coeffs, &poly1->exps,
                            &poly1->alloc, poly2->coeffs, exp2, poly2->length,
                              poly3->coeffs, exp3, poly3->length, exp_bits, N,
                                                               cmpmask);
   }

cleanup:

   _fmpz_mpoly_set_length(poly1, len, ctx);

   if (free2)
      flint_free(exp2);

   if (free3)
      flint_free(exp3);

   TMP_END;

   /* division is exact if len is nonzero */
   return (len != 0);
}
