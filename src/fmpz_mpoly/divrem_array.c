/*
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
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
   Set polyq to the quotient and polyr to the remainder of poly2 divided
   by poly3, and return the length of the quotient. The polynomials have
   their exponents tightly packed, with mixed bases equal to the largest
   exponent for each variable, e.g. the input polys have exponents of the
   form a_0 + a_1*b1 + a_2*b_2*b_2 + .... where b_0, b_1, b_2, etc, are
   the bases, which are equal to the largest possible exponent for
   each of the respective variables in the exponent. The dividend poly3
   is assumed to be nonzero. There are assumed to be "num" variables and
   the bases b_i are passed in the array "mults". The function reallocates
   its output. The quotient and remainder are written out in reverse order.
   The quotient and remainder poly are not assumed to be zero on input.
   The quotient and remainder terms are appended to the existing terms in
   those polys. 
*/
slong _fmpz_mpoly_divrem_array_tight(slong * lenr,
 fmpz ** polyq, ulong ** expq, slong * allocq, slong len0,
       fmpz ** polyr, ulong ** expr, slong * allocr, slong len1,
                  const fmpz * poly2, const ulong * exp2, slong len2,
                        const fmpz * poly3, const ulong * exp3, slong len3,
                                                      slong * mults, slong num)
{
   slong i, j, q, r, prod, bits1, bits2, bits3, k = len0, l = len1;
   slong max3 = (slong) exp3[0]; /* largest exponent in poly3 */
   slong * prods;
   fmpz c3 = poly3[0];
   /* abs val of leading coeff of poly3 */
   ulong u3 = ((ulong) FLINT_ABS(c3)) >> 1;
   fmpz * p1 = *polyq, * p2 = *polyr;
   ulong * e1 = *expq, * e2 = *expr;
   int small;
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

   prod = prods[num];

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
      ulong * t2 = (ulong *) TMP_ALLOC(2*prod*sizeof(ulong));

      for (i = 0; i < 2*prod; i++)
         t2[i] = 0;

      /* poly2 to array format */
      _fmpz_mpoly_to_ulong_array2(t2, poly2, exp2, len2);

      /* for each term of poly2 array relevant to quotient */
      for (i = prod - 1; i >= max3; i--)
      {
         ulong * ptr = t2 + 2*i;
         ulong p[2];

         /* if coeff is nonzero */
         if (ptr[0] != 0 || ptr[1] != 0)
         {
            if (0 > (slong) ptr[1])
               mpn_neg(p, ptr, 2);
            else
               flint_mpn_copyi(p, ptr, 2);

            /* not exact quotient monomial, thus remainder monomial */
            if (!mpoly_monomial_divides_tight(i, max3, prods, num))
            {
               /* realloc remainder poly */
               _fmpz_mpoly_fit_length(&p2, &e2, allocr, l + 1, 1);

               /* set remainder coeff... */
               fmpz_set_signed_uiui(p2 + l, ptr[1], ptr[0]);

               /* ...and exponent */
               e2[l++] = i;
            } else /* monomials can be divided exactly */
            {
               /* check quotient won't overflow a word */
               if (u3 <= p[1] || (u3 == 0 && 0 > (slong) p[0])) /* quotient too large */
               {
                  for (j = len0; j < k; j++)
                     _fmpz_demote(p1 + j);
                  for (j = len1; j < l; j++)
                     _fmpz_demote(p2 + j);
                  k = len0;
                  l = len1;

                  goto big;
               }

               /* quotient and remainder of coeffs */
               sdiv_qrnnd(q, r, ptr[1], ptr[0], c3);

               /* check coefficient is small, else restart with multiprec code */
               if (COEFF_IS_MPZ(FLINT_ABS(q))) /* quotient too large */
               {
                  for (j = len0; j < k; j++)
                     _fmpz_demote(p1 + j);
                  for (j = len1; j < l; j++)
                     _fmpz_demote(p2 + j);
                  k = len0;
                  l = len1;

                  goto big;
               }

               /* check coeff quotient was exact */
               if (r != 0) /* not an exact division */
               {
                  /* realloc remainder poly */
                  _fmpz_mpoly_fit_length(&p2, &e2, allocr, l + 1, 1);

                  /* set remainder coeff... */
                  fmpz_set_si(p2 + l, (slong) r);  

                  /* ... and exponent */
                  e2[l++] = i;
               }

               if (q != 0)
               {
                  /* submul a - q*b */
                  _fmpz_mpoly_submul_array1_slong2_1(t2, q, i - max3,
                                                            poly3, exp3, len3);

                  /* realloc quotient poly */
                  _fmpz_mpoly_fit_length(&p1, &e1, allocq, k + 1, 1);
                  /* set quotient coeff and exponent */
                  fmpz_set_si(p1 + k, q);
                  e1[k++] = i - max3;
               }
            }         
         }
      }

      /* all remaining terms are remainder terms */
      for ( ; i >= 0; i--)
      {
         ulong * ptr = t2 + 2*i;

         /* if coeff nonzero */
         if (ptr[0] != 0 || ptr[1] != 0)  /* not an exact division */
         {
            /* realloc remainder poly */
            _fmpz_mpoly_fit_length(&p2, &e2, allocr, l + 1, 1);

            /* set remainder coeff... */
            fmpz_set_signed_uiui(p2 + l, ptr[1], ptr[0]);

            /* and exponent */
            e2[l++] = i;
         }
      }
   }

   /* not done, coeffs small and intermediate computations fit three words */
   if (k == len0 && l == len1 && small)
   {
      ulong * t2 = (ulong *) TMP_ALLOC(3*prod*sizeof(ulong));

      for (i = 0; i < 3*prod; i++)
         t2[i] = 0;

      /* poly2 to array format */
      _fmpz_mpoly_to_ulong_array(t2, poly2, exp2, len2);

      /* for each term of poly2 array relevant to exact quotient */
      for (i = prod - 1; i >= max3; i--)
      {
         ulong * ptr = t2 + 3*i;
         ulong p[3];

         /* if coeff is nonzero */
         if (ptr[0] != 0 || ptr[1] != 0 || ptr[2] != 0)
         {
            if (0 > (slong) ptr[2])
               mpn_neg(p, ptr, 3);
            else
               flint_mpn_copyi(p, ptr, 3);

            /* not exact quotient monomial, thus remainder monomial */
            if (!mpoly_monomial_divides_tight(i, max3, prods, num))
            {
               /* realloc remainder poly */
               _fmpz_mpoly_fit_length(&p2, &e2, allocr, l + 1, 1);

               /* set remainder coeff... */
               fmpz_set_signed_uiuiui(p2 + l, ptr[2], ptr[1], ptr[0]);

               /* ... and exponent */
               e2[l++] = i;
            } else /* monomials can be divided exact */
            {
               /* check quotient won't overflow a word */
               if (p[2] > 0 || u3 <= p[1] || (u3 == 0 && 0 > (slong) p[0])) /* quotient too large */
               {
                  for (j = len0; j < k; j++)
                     _fmpz_demote(p1 + j);
                  for (j = len1; j < l; j++)
                     _fmpz_demote(p2 + j);
                  k = len0;
                  l = len1;

                  goto big;
               }

               /* quotient and remainder of coeffs */
               sdiv_qrnnd(q, r, ptr[1], ptr[0], c3);

               /* check coefficient is small, else restart with multiprec code */
               if (COEFF_IS_MPZ(FLINT_ABS(q))) /* quotient too large */
               {
                  for (j = len0; j < k; j++)
                     _fmpz_demote(p1 + j);
                  for (j = len1; j < l; j++)
                     _fmpz_demote(p2 + j);
                  k = len0;
                  l = len1;

                  goto big;
               }

               /* check if coeff quotient was exact */
               if (r != 0) /* else remainder term */ 
               {
                  /* reallocate remainder poly */
                  _fmpz_mpoly_fit_length(&p2, &e2, allocr, l + 1, 1);

                  /* set remainder coeff... */
                  fmpz_set_si(p2 + l, (slong) r);  

                  /* and exponent */
                  e2[l++] = i;
               }

               /* if nonzero quotient */
               if (q != 0)
               {
                  /* submul a - q*b */
                  _fmpz_mpoly_submul_array1_slong_1(t2, q, i - max3,
                                                            poly3, exp3, len3);

                  /* realloc quotient poly */
                  _fmpz_mpoly_fit_length(&p1, &e1, allocq, k + 1, 1);
                  
                  /* set quotient coeff and exponent */
                  fmpz_set_si(p1 + k, q);

                  e1[k++] = i - max3;
               }
            }
         }
      }

      /* all remaining terms are remainder terms */
      for ( ; i >= 0; i--)
      {
         ulong * ptr = t2 + 3*i;

         /* if coeff nonzero */
         if (ptr[0] != 0 || ptr[1] != 0 || ptr[2] != 0) 
         {
            /* not an exact division */

            /* realloc remainder poly */
            _fmpz_mpoly_fit_length(&p2, &e2, allocr, l + 1, 1);

            /* set remainder coeff... */
            fmpz_set_signed_uiuiui(p2 + l, ptr[2], ptr[1], ptr[0]);

            /* ...and exponent */
            e2[l++] = i;
         }
      }
   }

big:

   /* if still not done, use multiprecision coeffs instead */
   if (k == len0 && l == len1)
   {
      fmpz * t2 = (fmpz *) TMP_ALLOC(prod*sizeof(fmpz));
      fmpz_t fq, fr;

      fmpz_init(fq);
      fmpz_init(fr);

      for (i = 0; i < prod; i++)
         fmpz_init(t2 + i);

      /* poly2 to array format */
      _fmpz_mpoly_to_fmpz_array(t2, poly2, exp2, len2);
      
      /* for each term of poly2 array relevant to exact quotient */
      for (i = prod - 1; i >= max3; i--)
      {
         /* if coeff is nonzero */
         if (!fmpz_is_zero(t2 + i))
         {
            /* not exact quotient monomial, thus remainder monomial */
            if (!mpoly_monomial_divides_tight(i, max3, prods, num))
            {
               /* realloc remainder poly */
               _fmpz_mpoly_fit_length(&p2, &e2, allocr, l + 1, 1);

               /* set remainder coeff... */
               fmpz_set(p2 + l, t2 + i);
               
               /* ... and exponent */
               e2[l++] = i;
            } else /* monomials can be divided exactly */
            {
               /* quotient and remainder of coeffs */
               fmpz_fdiv_qr(fq, fr, t2 + i, poly3);

               /* check if coeff quotient was exact */
               if (!fmpz_is_zero(fr)) /* else remainder term */
               {
                  /* realloc remainder poly */
                  _fmpz_mpoly_fit_length(&p2, &e2, allocr, l + 1, 1);

                  /* set remainder coeff... */
                  fmpz_set(p2 + l, fr);  

                  /* and exponent */
                  e2[l++] = i;
               }

               /* if nonzero quotient */
               if (!fmpz_is_zero(fq))
               {
                  /* submul a - q*b */
                  _fmpz_mpoly_submul_array1_fmpz_1(t2, fq, i - max3,
                                                            poly3, exp3, len3);

                  /* realloc quotient poly */
                  _fmpz_mpoly_fit_length(&p1, &e1, allocq, k + 1, 1);
            
                  /* set quotient coeff and exponent */
                  fmpz_set(p1 + k, fq);
                  e1[k++] = i - max3;
               }
            }
         }
      }

      /* all remaining terms are remainder terms */
      for ( ; i >= 0; i--)
      {
         /* if coeff nonzero */
         if (!fmpz_is_zero(t2 + i))
         {
            /* remainder */

            /* realloc remainder poly */
            _fmpz_mpoly_fit_length(&p2, &e2, allocr, l + 1, 1);

            /* set remainder coeff... */
            fmpz_set(p2 + l, t2 + i);

            /* ... and exponent */
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

   /* set remainder poly length */
   (*lenr) = l - len1;

   TMP_END;

   /* return quotient poly length */
   return k - len0;
}

/*
   Use dense array division to set polyq, polyr to poly2/poly3 in num + 1
   variables, given a list of multipliers to tightly pack exponents and a
   number of bits for the fields of the exponents of the result, assuming
   no aliasing. classical exact division in main variable, array
   multiplication (submul) for multivariate coefficients in remaining num
   variables. The function reallocates its output and returns the length
   of the quotient poly. It is assumed that poly2 is not zero. The
   quotient and remainder are written in reverse order.
*/
slong _fmpz_mpoly_divrem_array_chunked(slong * lenr,
            fmpz ** polyq, ulong ** expq, slong * allocq,
                 fmpz ** polyr, ulong ** expr, slong * allocr, 
                const fmpz * poly2, const ulong * exp2, slong len2, 
        const fmpz * poly3, const ulong * exp3, slong len3, slong * mults, 
                                                         slong num, slong bits)
{
   slong i, j, k, l = 0, m, prod, len = 0, l1, l2, l3;
   slong bits1, bits2, bits3 = 0, tlen, talloc;
   slong shift = bits*(num);
   slong * i1, * i2, * i3, * n1, * n2, * n3, * prods;
   slong * b1, * b3, * maxb1, * maxb3, * max_exp1, * max_exp3;
   ulong * e2, * e3, * texp, * p2;
   fmpz * temp;
   int small;
   TMP_INIT;

   TMP_START;
   
   /*
      compute products 1, b0, b0*b1, b0*b1*b2 ...
      from list of bases b0, b1, b2, ...
   */
   prods = (slong *) TMP_ALLOC((num + 1)*sizeof(slong));

   prods[0] = 1;
   for (i = 0; i < num; i++)
      prods[i + 1] = prods[i]*mults[i];
   prod = prods[num];

   /* lengths of poly2, poly3 and polyq in chunks */
   l2 = 1 + (slong) (exp2[0] >> shift);
   l3 = 1 + (slong) (exp3[0] >> shift);

   l1 = FLINT_MAX(l2 - l3 + 1, 0);

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
   max_exp1 = (slong *) TMP_ALLOC(l1*num*sizeof(slong));
   max_exp3 = (slong *) TMP_ALLOC(l3*num*sizeof(slong));

   mpoly_main_variable_terms1(i2, n2, exp2, l2, len2, num + 1, num + 1, bits);
   mpoly_main_variable_terms1(i3, n3, exp3, l3, len3, num + 1, num + 1, bits);

   /* work out max bits for each coeff and optimal bits */

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

   /* work out maximum exponents for each chunk */
   for (i = 0; i < l3; i++)
      mpoly_max_degrees_tight(max_exp3 + i*num, e3 + i3[i], n3[i], prods, num);
   
   /* bound poly2 coeffs and check input/output coeffs likely small */
   bits2 = _fmpz_vec_max_bits(poly2, len2);
   /* we assume a bound of SMALL_FMPZ_BITCOUNT_MAX for coefficients of the quotient */
   bits1 = FLINT_ABS(bits3) + FLINT_BITS + FLINT_BIT_COUNT(len3) - 2;

   small = FLINT_ABS(bits2) <= bits1 && FLINT_ABS(bits3) <= SMALL_FMPZ_BITCOUNT_MAX;

   /* alloc space for copy of coeff/chunk of poly2 */

   temp = (fmpz *) flint_calloc(n2[0] + 1, sizeof(fmpz));
   texp = (ulong *) flint_malloc((n2[0] + 1)*sizeof(ulong));
   talloc = n2[0] + 1; /* plus one so doubling always increases size */

   /* enough space for three words per coeff, even if only one or two needed */
   p2 = (ulong *) TMP_ALLOC(3*prod*sizeof(ulong));

   /* coefficients likely to be small */
   if (small)
   {
      /* for each chunk of poly2 */
      for (i = 0; i < l2; i++)
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

               if (k < l3 && k >= 0)
               {
                  for (m = 0; m < num; m++)
                  {
                     if (max_exp1[j*num + m] + max_exp3[k*num + m] >= mults[m])
                     {
                        for (j = 0; j < len; j++)
                           _fmpz_demote((*polyq) + j);
                        for (j = 0; j < l; j++)
                           _fmpz_demote((*polyr) + j);
                        len = 0;
                        l = 0;
                        goto cleanup3;
                     }
                  }

                  bits1 = FLINT_MAX(bits1, FLINT_MIN(b1[j] + maxb3[k], maxb1[j] + b3[k]));
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

               if (k < l3 && k >= 0)
               {
                  _fmpz_mpoly_submul_array1_slong1(p2, (*polyq) + i1[j],
                     (*expq) + i1[j], n1[j], poly3 + i3[k], e3 + i3[k], n3[k]);
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

               if (k < l3 && k >= 0)
               {
                  _fmpz_mpoly_submul_array1_slong2(p2, (*polyq) + i1[j],
                     (*expq) + i1[j], n1[j], poly3 + i3[k], e3 + i3[k], n3[k]);
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

               if (k < l3 && k >= 0)
                  _fmpz_mpoly_submul_array1_slong(p2, (*polyq) + i1[j],
                     (*expq) + i1[j], n1[j], poly3 + i3[k], e3 + i3[k], n3[k]);
            }

            /* convert chunk from array format */
            tlen = _fmpz_mpoly_from_ulong_array(&temp, &texp, &talloc, 
                                                      p2, mults, num, bits, 0);
         }

         if (tlen != 0) /* nonzero coeff/chunk */
         {
            if (i < l1) /* potentially a quotient with remainder */
            {
               /* tightly pack chunk exponents */
               mpoly_pack_monomials_tight(texp, texp, tlen, mults, num, bits);

               /* set starting index for quotient chunk we are about to compute */
               i1[i] = len;

               /* compute quotient chunk and set length of quotient chunk */
	       n1[i] = _fmpz_mpoly_divrem_array_tight(lenr, polyq,
                 expq, allocq, len, polyr, expr, allocr, l, temp, texp,
                     tlen, poly3 + i3[0], e3 + i3[0], n3[0], mults, num);
																   
	       /* unpack remainder exponents */
	       mpoly_unpack_monomials_tight(*expr + l,
			                            *expr + l, *lenr, mults, num, bits);
														  
	       /* insert main variable */
	       for (j = 0; j < *lenr; j++)
	          (*expr)[l + j] += ((l2 - i - 1) << shift);
				  
	       /* work out maximum exponents for chunk */
               mpoly_max_degrees_tight(max_exp1 + i*num,
	                                   (*expq) + i1[i], n1[i], prods, num);
						 
	       /* check there were no exponent overflows */
	       for (m = 0; m < num; m++)
	       {
	          if (max_exp3[m] + max_exp1[i*num + m] >= mults[m])
                  {
                     for (j = 0; j < len; j++)
                        _fmpz_demote((*polyq) + j);
                     for (j = 0; j < l; j++)
                        _fmpz_demote((*polyr) + j);
                     len = 0;
		     l = 0;

	             goto cleanup3;
                  }				  
	       }

               /* check the quotient didn't have large coefficients */
               if (FLINT_ABS(_fmpz_vec_max_bits((*polyq) + len,
                                                      n1[i])) > SMALL_FMPZ_BITCOUNT_MAX)
               {
                  for (j = 0; j < len; j++)
                     _fmpz_demote((*polyq) + j);
                  for (j = 0; j < l; j++)
                     _fmpz_demote((*polyr) + j);
                  len = 0;
                  l = 0;

                  goto big;
               }

               /* abs bound and sum of abs vals of coeffs of quotient chunk */
               _fmpz_vec_sum_max_bits(&b1[i], &maxb1[i], (*polyq)+i1[i], n1[i]);

               /* update length of output quotient and remainder polys */
               len += n1[i];
               l += *lenr;
            } else /* remainder terms only */
            {
               /* realloc remainder poly */
               _fmpz_mpoly_fit_length(polyr, expr, allocr, l + tlen, 1);

               /* for each term in remainder chunk */
               for (j = 0; j < tlen; j++)
               {
                  /* set remainder coeff and exponent */
                  fmpz_set(*polyr + l + j, temp + j);
                  (*expr)[l + j] = (texp[j]) + ((l2 - i - 1) << shift);
               }

               l += tlen;
            }
         } else if (i < l1) /* zero chunk, no quotient or remainder */
         {
            /* set index and length of quotient chunk */
            i1[i] = len;
            n1[i] = 0;
            b1[i] = 0;
            maxb1[i] = 0;

            /* write out maximum exponents in chunk */
            mpoly_max_degrees_tight(max_exp1 + i*num,
	                                   (*expq) + i1[i], n1[i], prods, num);
         }
      }
   }

big:

   /* if not done, use multiprecision coeffs instead */
   if (len == 0 && l == 0)
   {
      fmpz * p2 = (fmpz *) TMP_ALLOC(prod*sizeof(fmpz));

      for (j = 0; j < prod; j++)
            fmpz_init(p2 + j);
      
      /* for each chunk of poly2 */
      for (i = 0; i < l2; i++)
      {
         for (j = 0; j < prod; j++)
            fmpz_zero(p2 + j);

         /* convert relevant coeff/chunk of poly2 to array format */
         _fmpz_mpoly_to_fmpz_array(p2, poly2 + i2[i], e2 + i2[i], n2[i]);

         /* submuls */

         for (j = 0; j < i && j < l1; j++)
         {
            k = i - j;

            if (k < l3 && k >= 0)
            {
	       for (m = 0; m < num; m++)
	       {
	          if (max_exp1[j*num + m] + max_exp3[k*num + m] >= mults[m])
	          {
                     for (j = 0; j < len; j++)
                        _fmpz_demote((*polyq) + j);
                     for (j = 0; j < l; j++)
                        _fmpz_demote((*polyr) + j);
		     len = 0;
	             l = 0;
						
                for (j = 0; j < prod; j++)
                   fmpz_clear(p2 + j);
		     goto cleanup3;
	          }
	       }

               _fmpz_mpoly_submul_array1_fmpz(p2, (*polyq) + i1[j],
                    (*expq) + i1[j], n1[j], poly3 + i3[k], e3 + i3[k], n3[k]);
	    }
         }

         /* convert chunk from array format */
         tlen = _fmpz_mpoly_from_fmpz_array(&temp, &texp, &talloc, 
                                                      p2, mults, num, bits, 0);

         if (tlen != 0) /* nonzero coeff/chunk */
         {
            if (i < l1) /* potentially a quotient with remainder */
            {
               /* tightly pack chunk exponents */
               mpoly_pack_monomials_tight(texp, texp, tlen, mults, num, bits);

               /* set start index of quotient chunk we are about to compute */
               i1[i] = len;
            
               /* compute quotient chunk and set length of quotient chunk */
               n1[i] = _fmpz_mpoly_divrem_array_tight(lenr, polyq,
                      expq, allocq, len, polyr, expr, allocr, l, temp, texp, 
                           tlen, poly3 + i3[0], e3 + i3[0], n3[0], mults, num);

	       /* unpack remainder exponents */
	       mpoly_unpack_monomials_tight(*expr + l,
			                *expr + l, *lenr, mults, num, bits);
														  
	       /* insert main variable */
	       for (j = 0; j < *lenr; j++)
			             (*expr)[l + j] += ((l2 - i - 1) << shift);

	       /* work out maximum exponents for chunk */
               mpoly_max_degrees_tight(max_exp1 + i*num,
	                                   (*expq) + i1[i], n1[i], prods, num);

	       /* check there were no exponent overflows */
	       for (m = 0; m < num; m++)
	       {
	          if (max_exp3[m] + max_exp1[i*num + m] >= mults[m])
                  {
                     for (j = 0; j < len; j++)
                        _fmpz_demote((*polyq) + j);
                     for (j = 0; j < l; j++)
                        _fmpz_demote((*polyr) + j);
                     len = 0;
	             l = 0;

                for (j = 0; j < prod; j++)
                   fmpz_clear(p2 + j);
	             goto cleanup3;
                  }				  
	       }

               /* abs bound and sum of abs vals of coeffs of quotient chunk */
               _fmpz_vec_sum_max_bits(&b1[i], &maxb1[i], (*polyq)+i1[i], n1[i]);

               /* update length of output quotient and remainder polys */
               len += n1[i];
               l += *lenr;
            } else /* remainder terms only */
            {
               /* realloc remainder poly */
               _fmpz_mpoly_fit_length(polyr, expr, allocr, l + tlen, 1);

               /* for each term in chunk */
               for (j = 0; j < tlen; j++)
               {
                  /* set remainder coeff and exponent */
                  fmpz_set(*polyr + l + j, temp + j);
                  (*expr)[l + j] = (texp[j]) + ((l2 - i - 1) << shift);
               }

               /* update length of output remainder poly */
               l += tlen;
            }
         } else if (i < l1) /* zero chunk, no quotient or remainder */
         {
            /* set index and length of quotient chunk */
            i1[i] = len;
            n1[i] = 0;
            b1[i] = 0;
            maxb1[i] = 0;

            /* write out maximum exponents in chunk */
            mpoly_max_degrees_tight(max_exp1 + i*num,
	                                   (*expq) + i1[i], n1[i], prods, num);
         }
      }

      for (j = 0; j < prod; j++)
            fmpz_clear(p2 + j);
   }

   /* if there were quotient terms */
   if (len != 0)
   {
      /* unpack monomials of quotient */
      mpoly_unpack_monomials_tight((*expq), (*expq), len, mults, num, bits);

      /* put main variable back in quotient */
      for (i = 0; i < l1; i++)
      {
         for (j = 0; j < n1[i]; j++)
         {
            (*expq)[i1[i] + j] += ((l1 - i - 1) << shift);
         }
      }
   }

cleanup3:

   for (j = 0; j < talloc; j++)
      fmpz_clear(temp + j);

   flint_free(temp);
   flint_free(texp);

   TMP_END;

   /* set remainder length */
   *lenr = l;

   /* return quotient length */
   return len;
}

/*
   Use dense array division to set polyq, polyr to poly2/poly3 in num variables,
   given a list of multipliers to tightly pack exponents and a number of bits
   for the fields of the exponents of the result, assuming no aliasing. The
   array "mults" is a list of bases to be used in encoding the array indices
   from the exponents. The function reallocates its output and returns the
   length of the quotient. It is assumed that poly2 is not zero. The quotient
   and remainder are written in reverse order.
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
   slong * prods, * max_exp1, * max_exp3;
   TMP_INIT;

   TMP_START;

   /*
      compute products 1, b0, b0*b1, b0*b1*b2 ...
      from list of bases b0, b1, b2, ...
   */
   prods = (slong *) TMP_ALLOC((num + 1)*sizeof(slong));
  
   prods[0] = 1;
   for (i = 0; i < num; i++)
      prods[i + 1] = prods[i]*mults[i];
   prod = prods[num];

   /* if array size will be too large, chunk the polynomials */
   if (prod > MAX_ARRAY_SIZE)
   {
      TMP_END;
	  
	  return _fmpz_mpoly_divrem_array_chunked(lenr, polyq, expq, allocq,
                                     polyr, expr, allocr, poly2, exp2, len2,
                                      poly3, exp3, len3, mults, num - 1, bits);
   }
									  
   e2 = (ulong *) TMP_ALLOC(len2*sizeof(ulong));
   e3 = (ulong *) TMP_ALLOC(len3*sizeof(ulong));
   max_exp1 = (slong *) TMP_ALLOC(num*sizeof(slong));
   max_exp3 = (slong *) TMP_ALLOC(num*sizeof(slong));

   /* pack input exponents tightly with mixed bases specified by "mults" */

   mpoly_pack_monomials_tight(e2, exp2, len2, mults, num, bits);
   mpoly_pack_monomials_tight(e3, exp3, len3, mults, num, bits);

   /* do divrem on tightly packed polys */
   len = _fmpz_mpoly_divrem_array_tight(lenr, polyq, expq, allocq, 0,
                                  polyr, expr, allocr, 0, poly2, e2, len2,
                                                  poly3, e3, len3, mults, num);

   /* check for overflows */
   mpoly_max_degrees_tight(max_exp3, e3, len3, prods, num);
   mpoly_max_degrees_tight(max_exp1, *expq, len, prods, num);
   
   for (i = 0; i < num; i++)
   {
      if (max_exp3[i] + max_exp1[i] >= mults[i])
      {
	 len = 0;
         *lenr = 0;
		 
         break;
      }
   }
   /* unpack output quotient and remainder exponents */
   mpoly_unpack_monomials_tight((*expq), (*expq), len, mults, num, bits);
   mpoly_unpack_monomials_tight((*expr), (*expr), *lenr, mults, num, bits);

   TMP_END;

   return len;
}

/*
   Return 1 if q, r can be set to quotient and remainder of poly2 by poly3,
   else return 0 if array division not able to be performed.
*/
int fmpz_mpoly_divrem_array(fmpz_mpoly_t q, fmpz_mpoly_t r,
                    const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3, 
                                                    const fmpz_mpoly_ctx_t ctx)
{
   slong i, exp_bits, N, lenq = 0, lenr = 0, array_size;
   ulong * max_fields, * max_fields2, * max_fields3;
   ulong * exp2 = poly2->exps, * exp3 = poly3->exps;
   int free2 = 0, free3 = 0;
   fmpz_mpoly_t temp1, temp2;
   fmpz_mpoly_struct * tq, * tr;
   int res = 0;

   TMP_INIT;

   /* check divisor is not zero */
   if (poly3->length == 0)
      flint_throw(FLINT_DIVZERO, "Divide by zero in fmpz_mpoly_divrem_array");

   /* dividend is zero */
   if (poly2->length == 0)
   {
      fmpz_mpoly_zero(q, ctx);
      fmpz_mpoly_zero(r, ctx);

      return 1;
   }

   TMP_START;


   /* compute maximum exponents for each variable */
   max_fields = (ulong *) TMP_ALLOC(ctx->minfo->nfields*sizeof(ulong));
   max_fields2 = (ulong *) TMP_ALLOC(ctx->minfo->nfields*sizeof(ulong));
   max_fields3 = (ulong *) TMP_ALLOC(ctx->minfo->nfields*sizeof(ulong));
   mpoly_max_fields_ui_sp(max_fields2, poly2->exps, poly2->length,
                                                      poly2->bits, ctx->minfo);
   mpoly_max_fields_ui_sp(max_fields3, poly3->exps, poly3->length,
                                                      poly3->bits, ctx->minfo);
   for (i = 0; i < ctx->minfo->nfields; i++)
      max_fields[i] = FLINT_MAX(max_fields2[i], max_fields3[i]);

   /* compute number of bits required for output exponents */
   exp_bits = FLINT_MAX(poly2->bits, poly3->bits);
   N = mpoly_words_per_exp(exp_bits, ctx->minfo);

   /* array division expects each exponent vector in one word */
   /* current code is wrong for reversed orderings */
   if (N != 1 || mpoly_ordering_isrev(ctx->minfo))
      goto cleanup;

   /* compute bounds on output exps, used as mixed bases for packing exps */
   array_size = 1;
   for (i = 0; i < ctx->minfo->nfields - 1; i++)
   {
      max_fields2[i] = max_fields[i] + 1;
      array_size *= max_fields2[i];
   }  
   max_fields2[ctx->minfo->nfields - 1] = max_fields[ctx->minfo->nfields - 1] + 1;
   
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

   if (exp2[0] < exp3[0])
   {
      fmpz_mpoly_set(r, poly2, ctx);
      fmpz_mpoly_zero(q, ctx);
	  
      res = 1;
	  
      goto cleanup2;
   }

   /* handle aliasing and do array division */

   if (q == poly2 || q == poly3)
   {
      fmpz_mpoly_init2(temp1, poly2->length/poly3->length + 1, ctx);
      fmpz_mpoly_fit_bits(temp1, exp_bits, ctx);
      temp1->bits = exp_bits;

      tq = temp1;
   } else
   {
      fmpz_mpoly_fit_length(q, poly2->length/poly3->length + 1, ctx);
      fmpz_mpoly_fit_bits(q, exp_bits, ctx);
      q->bits = exp_bits;

      tq = q;
   }

   if (r == poly2 || r == poly3)
   {
      fmpz_mpoly_init2(temp2, poly3->length, ctx);
      fmpz_mpoly_fit_bits(temp2, exp_bits, ctx);
      temp2->bits = exp_bits;

      tr = temp2;
   } else
   {
      fmpz_mpoly_fit_length(r, poly3->length, ctx);
      fmpz_mpoly_fit_bits(r, exp_bits, ctx);
      r->bits = exp_bits;

      tr = r;
   }

   lenq = _fmpz_mpoly_divrem_array(&lenr, &tq->coeffs, &tq->exps,
        &tq->alloc, &tr->coeffs, &tr->exps, &tr->alloc, poly2->coeffs,
                  exp2, poly2->length, poly3->coeffs, exp3, poly3->length,
                         (slong *) max_fields2, ctx->minfo->nfields, exp_bits);

    res = (lenq != 0 || lenr != 0);

    if (res)
    {
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
    }
    else
    {
        if (q == poly2 || q == poly3)
        {
            fmpz_mpoly_clear(temp1, ctx);
        }

        if (r == poly2 || r == poly3)
        {
            fmpz_mpoly_clear(temp2, ctx);
        }

        for (i = q->length; i < q->alloc; i++)
        {
            _fmpz_demote(q->coeffs + i);
        }
        for (i = r->length; i < r->alloc; i++)
        {
            _fmpz_demote(r->coeffs + i);
        }
    }

    _fmpz_mpoly_set_length(q, lenq, ctx);
    _fmpz_mpoly_set_length(r, lenr, ctx);


cleanup2:

   if (free2)
      flint_free(exp2);

   if (free3)
      flint_free(exp3);

cleanup:

   TMP_END;

   return res;
}
