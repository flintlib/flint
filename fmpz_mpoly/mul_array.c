/*
    Copyright (C) 2016 William Hart
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

/*
    NOTE: this file is dirty - it assumes that a zero fmpz is zero
*/

/* improve locality */
#define BLOCK 128
#define MAX_ARRAY_SIZE (WORD(300000))
#define MAX_LEX_SIZE (WORD(300))

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




/****************************************************
    LEX
****************************************************/

#define LEX_UNPACK_MACRO(fxn_name, coeff_decl, nonzero_test, swapper)          \
slong fxn_name(fmpz_mpoly_t P, slong Plen, coeff_decl,                         \
           const ulong * mults, slong num, slong array_size, slong top)        \
{                                                                              \
    slong off, j;                                                              \
    slong topmult = num == 0 ? 1 : mults[num - 1];                             \
    slong lastd   = topmult - 1;                                               \
    slong reset   = array_size/topmult;                                        \
    slong counter = reset;                                                     \
    ulong startexp = (top << (P->bits*num)) + (lastd << (P->bits*(num-1)));    \
    for (off = array_size - 1; off >= 0; off--)                                \
    {                                                                          \
        if (nonzero_test)                                                      \
        {                                                                      \
            slong d = off;                                                     \
            ulong exp = startexp;                                              \
            for (j = 0; j + 1 < num; j++) {                                    \
                exp += (d % mults[j]) << (P->bits*j);                          \
                d = d / mults[j];                                              \
            }                                                                  \
            _fmpz_mpoly_fit_length(&P->coeffs, &P->exps, &P->alloc, Plen + 1, 1); \
            P->exps[Plen] = exp;                                               \
            swapper                                                            \
            Plen++;                                                            \
        }                                                                      \
        counter--;                                                             \
        if (counter <= 0) {                                                    \
            counter = reset;                                                   \
            lastd--;                                                           \
            startexp -= UWORD(1) << (P->bits*(num-1));                         \
        }                                                                      \
    }                                                                          \
    return Plen;                                                               \
}

/*
    These four functions will replace
        _fmpz_mpoly_from_ulong_array, ..., _fmpz_mpoly_from_fmpz_array
    defined above.

    WARNING: In order for the coeff_array to be used interchangeably between
    the small versions and fmpz version, it is necessary that zero is
    represented as a real zero in the fmpz type.
*/

LEX_UNPACK_MACRO(
    fmpz_mpoly_append_array_sm1_LEX, ulong * coeff_array
,
    coeff_array[off] != WORD(0)
,
    fmpz_set_si(P->coeffs + Plen, coeff_array[off]);
    coeff_array[off] = 0;
)

LEX_UNPACK_MACRO(
    fmpz_mpoly_append_array_sm2_LEX, ulong * coeff_array
,
    (coeff_array[2*off+0] || coeff_array[2*off+1]) != WORD(0)
,
    fmpz_set_signed_uiui(P->coeffs + Plen, coeff_array[2*off+1], coeff_array[2*off+0]);
    coeff_array[2*off+0] = coeff_array[2*off+1] = 0;
)

LEX_UNPACK_MACRO(
    fmpz_mpoly_append_array_sm3_LEX, ulong * coeff_array
,
    (coeff_array[3*off+0] || coeff_array[3*off+1] || coeff_array[3*off+2]) != WORD(0)
,
    fmpz_set_signed_uiuiui(P->coeffs + Plen, coeff_array[3*off+2], coeff_array[3*off+1], coeff_array[3*off+0]);
    coeff_array[3*off+0] = coeff_array[3*off+1] = coeff_array[3*off+2] = 0;
)

LEX_UNPACK_MACRO(
    fmpz_mpoly_append_array_fmpz_LEX, fmpz * coeff_array
,
    !fmpz_is_zero(coeff_array + off)
,
    fmpz_swap(P->coeffs + Plen, coeff_array + off);
    fmpz_zero(coeff_array + off);
)



void _fmpz_mpoly_mul_array_chunked_LEX(
    fmpz_mpoly_t P,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const ulong * mults,
    const fmpz_mpoly_ctx_t ctx)
{
    slong num = ctx->minfo->nfields - 1;
    slong Pi, i, j, Plen, Pl, Al, Bl, array_size;
    slong Abits, * Asum, * Amax, Bbits, * Bsum, * Bmax;
    slong * Amain, * Bmain;
    ulong * Apexp, * Bpexp;
    int small;
    TMP_INIT;

    array_size = 1;
    for (i = 0; i < num; i++) {
        array_size *= mults[i];
    }

    /* compute lengths of poly2 and poly3 in chunks */
    Al = 1 + (slong) (A->exps[0] >> (A->bits*num));
    Bl = 1 + (slong) (B->exps[0] >> (B->bits*num));

    TMP_START;

    /* compute indices and lengths of coefficients of polys in main variable */
    Amain = (slong *) TMP_ALLOC(((Al + 1) + Al + Al + (Bl + 1) + Bl + Bl)*sizeof(slong));
    Asum  = Amain + Al + 1;
    Amax  = Asum + Al;
    Bmain = Amax + Al;
    Bsum  = Bmain + Bl + 1;
    Bmax  = Bsum + Bl;
    Apexp = (ulong *) flint_malloc(A->length*sizeof(ulong));
    Bpexp = (ulong *) flint_malloc(B->length*sizeof(ulong));
    mpoly_main_variable_split_LEX(Amain, Apexp, A->exps, Al, A->length, mults, num, A->bits);
    mpoly_main_variable_split_LEX(Bmain, Bpexp, B->exps, Bl, B->length, mults, num, B->bits);

    /* work out bit counts for each chunk */
    Abits = 0;
    for (i = 0; i < Al; i++)
    {
        _fmpz_vec_sum_max_bits(&Asum[i], &Amax[i], A->coeffs + Amain[i], Amain[i + 1] - Amain[i]);
        Abits = FLINT_MAX(Abits, Amax[i]);
    }

    Bbits = 0;
    for (j = 0; j < Bl; j++)
    {
        _fmpz_vec_sum_max_bits(&Bsum[j], &Bmax[j], B->coeffs + Bmain[j], Bmain[j + 1] - Bmain[j]);
        Bbits = FLINT_MAX(Bbits, Bmax[j]);
    }

    /* whether the output coefficients are "small" */
    small = Abits <= (SMALL_FMPZ_BITCOUNT_MAX) && Bbits <= (SMALL_FMPZ_BITCOUNT_MAX);

    Pl = Al + Bl - 1;
    Plen = 0;

    if (small)
    {
        ulong * coeff_array = (ulong *) TMP_ALLOC(3*array_size*sizeof(ulong));
        for (j = 0; j < 3*array_size; j++)
            coeff_array[j] = 0;

        /* for each output chunk */
        for (Pi = 0; Pi < Pl; Pi++)
        {
            /* compute bound on coeffs of output chunk */
            slong number = 0;
            slong Pbits = 0;
            for (i = 0, j = Pi; i < Al && j >= 0; i++, j--)
            {
                if (j < Bl)
                {
                    Pbits = FLINT_MAX(Pbits,
                              FLINT_MIN(Asum[i] + Bmax[j], Amax[i] + Bsum[j]));
                    number++;
                }
            }
            Pbits += FLINT_BIT_COUNT(number) + 1; /* includes one bit for sign */

            if (Pbits <= FLINT_BITS)
            {
                /* fits into one word */
                for (i = 0, j = Pi; i < Al && j >= 0; i++, j--)
                {
                    if (j < Bl)
                    {
                        _fmpz_mpoly_addmul_array1_slong1(coeff_array, 
                            (slong *) A->coeffs + Amain[i],
                                Apexp + Amain[i], Amain[i + 1] - Amain[i],
                            (slong *) B->coeffs + Bmain[j],
                                Bpexp + Bmain[j], Bmain[j + 1] - Bmain[j]);
                    }
                }

                Plen = fmpz_mpoly_append_array_sm1_LEX(P, Plen, coeff_array,
                                          mults, num, array_size, Pl - Pi - 1);

            } else if (Pbits <= 2*FLINT_BITS)
            {
                /* fits into two words */
                for (i = 0, j = Pi; i < Al && j >= 0; i++, j--)
                {
                    if (j < Bl)
                    {
                        _fmpz_mpoly_addmul_array1_slong2(coeff_array, 
                            (slong *) A->coeffs + Amain[i],
                                    Apexp + Amain[i], Amain[i + 1] - Amain[i],
                            (slong *) B->coeffs + Bmain[j],
                                    Bpexp + Bmain[j], Bmain[j + 1] - Bmain[j]);
                    }
                }

                Plen = fmpz_mpoly_append_array_sm2_LEX(P, Plen, coeff_array,
                                          mults, num, array_size, Pl - Pi - 1);

            } else 
            {
                /* fits into three word */
                for (i = 0, j = Pi; i < Al && j >= 0; i++, j--)
                {
                    if (j < Bl)
                    {
                        _fmpz_mpoly_addmul_array1_slong(coeff_array, 
                            (slong *) A->coeffs + Amain[i],
                                    Apexp + Amain[i], Amain[i + 1] - Amain[i],
                            (slong *) B->coeffs + Bmain[j],
                                    Bpexp + Bmain[j], Bmain[j + 1] - Bmain[j]);
                    }
                }

                Plen = fmpz_mpoly_append_array_sm3_LEX(P, Plen, coeff_array,
                                          mults, num, array_size, Pl - Pi - 1);
            }
        }
    } else
    {

        fmpz * coeff_array = (fmpz *) TMP_ALLOC(array_size*sizeof(fmpz));
        for (j = 0; j < array_size; j++)
            fmpz_init(coeff_array + j);

        /* for each output chunk */
        for (Pi = 0; Pi < Pl; Pi++)
        {
            /* addmuls for each cross product of chunks */
            for (i = 0, j = Pi; i < Al && j >= 0; i++, j--)
            {
                if (j < Bl)
                {
                    _fmpz_mpoly_addmul_array1_fmpz(coeff_array, 
                        A->coeffs + Amain[i],
                            Apexp + Amain[i], Amain[i + 1] - Amain[i],
                        B->coeffs + Bmain[j],
                            Bpexp + Bmain[j], Bmain[j + 1] - Bmain[j]);
                }
            }

            Plen = fmpz_mpoly_append_array_fmpz_LEX(P, Plen, coeff_array,
                                          mults, num, array_size, Pl - Pi - 1);
        }
    }

    _fmpz_mpoly_set_length(P, Plen, ctx);

    flint_free(Apexp);
    flint_free(Bpexp);
    TMP_END;
}



int _fmpz_mpoly_mul_array_LEX(
    fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    fmpz * maxBfields,
    const fmpz_mpoly_t C,
    fmpz * maxCfields,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, exp_bits, array_size;
    ulong max, * mults;
    int success;
    TMP_INIT;

    FLINT_ASSERT(ctx->minfo->nvars > 0);
    FLINT_ASSERT(B->length != 0);
    FLINT_ASSERT(C->length != 0);

    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    FLINT_ASSERT(1 == mpoly_words_per_exp(B->bits, ctx->minfo));
    FLINT_ASSERT(1 == mpoly_words_per_exp(C->bits, ctx->minfo));

    TMP_START;

    /* compute maximum exponents for each variable */
    mults = (ulong *) TMP_ALLOC(ctx->minfo->nfields*sizeof(ulong));

    /* the field of index n-1 is the one that wil be pulled out */
    i = ctx->minfo->nfields - 1;
    FLINT_ASSERT(fmpz_fits_si(maxBfields + i));
    FLINT_ASSERT(fmpz_fits_si(maxCfields + i));
    mults[i] = 1 + fmpz_get_ui(maxBfields + i) + fmpz_get_ui(maxCfields + i);
    max = mults[i];
    if (((slong) mults[i]) <= 0 || mults[i] > MAX_LEX_SIZE)
    {
        success = 0;
        goto cleanup;
    }

    /* the fields of index n-2...0, contribute to the array size */
    array_size = WORD(1);
    for (i--; i >= 0; i--)
    {
        ulong hi;
        FLINT_ASSERT(fmpz_fits_si(maxBfields + i));
        FLINT_ASSERT(fmpz_fits_si(maxCfields + i));
        mults[i] = 1 + fmpz_get_ui(maxBfields + i) + fmpz_get_ui(maxCfields + i);
        max |= mults[i];
        umul_ppmm(hi, array_size, array_size, mults[i]);
        if (hi != 0 || (slong) mults[i] <= 0
                    || array_size <= 0
                    || array_size > MAX_ARRAY_SIZE)
        {
            success = 0;
            goto cleanup;
        }
    }

    exp_bits = FLINT_MAX(MPOLY_MIN_BITS, FLINT_BIT_COUNT(max) + 1);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);

    /* array multiplication assumes result fits into 1 word */
    if (1 != mpoly_words_per_exp(exp_bits, ctx->minfo))
    {
        success = 0;
        goto cleanup;
    }

    /* handle aliasing and do array multiplication */
    if (A == B || A == C)
    {
        fmpz_mpoly_t T;
        fmpz_mpoly_init3(T, B->length + C->length, exp_bits, ctx);
        _fmpz_mpoly_mul_array_chunked_LEX(T, C, B, mults, ctx);
        fmpz_mpoly_swap(T, A, ctx);
        fmpz_mpoly_clear(T, ctx);
    }
    else
    {
        fmpz_mpoly_fit_length_reset_bits(A, B->length + C->length, exp_bits, ctx);
        _fmpz_mpoly_mul_array_chunked_LEX(A, C, B, mults, ctx);
    }
    success = 1;

cleanup:

    TMP_END;

    return success;
}




/****************************************************
    DEGLEX and DEGREVLEX
****************************************************/

#define DEGLEX_UNPACK_MACRO(fxn_name, coeff_decl, nonzero_test, swapper)       \
slong fxn_name(fmpz_mpoly_t P, slong Plen, coeff_decl,                         \
                                     slong top, slong nvars, slong degb)       \
{                                                                              \
    slong i;                                                                   \
    ulong exp, lomask = (UWORD(1) << (P->bits - 1)) - 1;                       \
    slong off, array_size;                                                     \
    slong * curexp, * degpow;                                                  \
    ulong * oneexp;                                                            \
    int carry;                                                                 \
    TMP_INIT;                                                                  \
                                                                               \
    TMP_START;                                                                 \
    curexp = (slong *) TMP_ALLOC(nvars*sizeof(slong));                         \
    degpow = (slong *) TMP_ALLOC(nvars*sizeof(slong));                         \
    oneexp = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));                         \
    array_size = 1;                                                            \
    curexp[0] = 0;                                                             \
    oneexp[0] = 0;                                                             \
    degpow[0] = 1;                                                             \
    for (i = 0; i < nvars-1; i++)                                              \
    {                                                                          \
        curexp[i] = 0;                                                         \
        degpow[i] = array_size;                                                \
        oneexp[i] = (UWORD(1) << (P->bits*(i+1))) - UWORD(1);                  \
        array_size *= degb;                                                    \
    }                                                                          \
    off = 0;                                                                   \
    if (nvars > 1)                                                             \
    {                                                                          \
        curexp[nvars - 2] = top;                                               \
        off = top * degpow[nvars - 2];                                         \
    }                                                                          \
    exp = (top << (P->bits*nvars)) + (top << (P->bits*(nvars-1)));             \
                                                                               \
    carry = 1;                                                                 \
    do {                                                                       \
        if (nonzero_test)                                                      \
        {                                                                      \
            _fmpz_mpoly_fit_length(&P->coeffs, &P->exps, &P->alloc, Plen + 1, 1); \
            P->exps[Plen] = exp;                                               \
            swapper                                                            \
            Plen++;                                                            \
        }                                                                      \
                                                                               \
        exp -= oneexp[0];                                                      \
        off -= 1;                                                              \
        curexp[0] -= 1;                                                        \
        if (curexp[0] >= 0)                                                    \
        {                                                                      \
            carry = 0;                                                         \
        } else                                                                 \
        {                                                                      \
            exp -= curexp[0]*oneexp[0];                                        \
            off -= curexp[0];                                                  \
            curexp[0] = 0;                                                     \
            carry = 1;                                                         \
                                                                               \
            for (i = 1; i < nvars - 1; i++)                                    \
            {                                                                  \
                exp -= oneexp[i];                                              \
                off -= degpow[i];                                              \
                curexp[i] -= 1;                                                \
                if (curexp[i] < 0)                                             \
                {                                                              \
                    exp -= curexp[i]*oneexp[i];                                \
                    off -= curexp[i]*degpow[i];                                \
                    curexp[i] = 0;                                             \
                    carry = 1;                                                 \
                } else                                                         \
                {                                                              \
                    ulong t = exp & lomask;                                    \
                    off += t*degpow[i - 1];                                    \
                    curexp[i - 1] = t;                                         \
                    exp += t*oneexp[i - 1];                                    \
                    carry = 0;                                                 \
                    break;                                                     \
                }                                                              \
            }                                                                  \
        }                                                                      \
    } while (!carry);                                                          \
                                                                               \
    TMP_END;                                                                   \
                                                                               \
    return Plen;                                                               \
}

DEGLEX_UNPACK_MACRO(
    fmpz_mpoly_append_array_sm1_DEGLEX, ulong * coeff_array
,
    coeff_array[off] != WORD(0)
,
    fmpz_set_si(P->coeffs + Plen, coeff_array[off]);
    coeff_array[off] = 0;
)

DEGLEX_UNPACK_MACRO(
    fmpz_mpoly_append_array_sm2_DEGLEX, ulong * coeff_array
,
    (coeff_array[2*off + 0] || coeff_array[2*off + 1]) != WORD(0)
,
    fmpz_set_signed_uiui(P->coeffs + Plen, coeff_array[2*off + 1], coeff_array[2*off + 0]);
    coeff_array[2*off + 0] = coeff_array[2*off + 1] = 0;
)

DEGLEX_UNPACK_MACRO(
    fmpz_mpoly_append_array_sm3_DEGLEX, ulong * coeff_array
,
    (coeff_array[3*off + 0] || coeff_array[3*off + 1] || coeff_array[3*off + 2]) != WORD(0)
,
    fmpz_set_signed_uiuiui(P->coeffs + Plen, coeff_array[3*off + 2], coeff_array[3*off + 1], coeff_array[3*off + 0]);
    coeff_array[3*off + 0] = coeff_array[3*off + 1] = coeff_array[3*off + 2] = 0;
)

DEGLEX_UNPACK_MACRO(
    fmpz_mpoly_append_array_fmpz_DEGLEX, fmpz * coeff_array
,
    !fmpz_is_zero(coeff_array + off)
,
    fmpz_swap(P->coeffs + Plen, coeff_array + off);
    fmpz_zero(coeff_array + off);
)





#define DEGREVLEX_UNPACK_MACRO(fxn_name, coeff_decl, nonzero_test, swapper)    \
slong fxn_name(fmpz_mpoly_t P, slong Plen, coeff_decl,                         \
                                     slong top, slong nvars, slong degb)       \
{                                                                              \
    slong i;                                                                   \
    ulong exp, mask = UWORD(1) << (P->bits - 1);                               \
    slong off, array_size;                                                     \
    slong * curexp, * degpow;                                                  \
    ulong * oneexp;                                                            \
    int carry;                                                                 \
    TMP_INIT;                                                                  \
                                                                               \
    TMP_START;                                                                 \
    curexp = (slong *) TMP_ALLOC(nvars*sizeof(slong));                         \
    degpow = (slong *) TMP_ALLOC(nvars*sizeof(slong));                         \
    oneexp = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));                         \
    array_size = 1;                                                            \
    oneexp[0] = 0;                                                             \
    for (i = 0; i < nvars-1; i++) {                                            \
        curexp[i] = 0;                                                         \
        degpow[i] = array_size;                                                \
        oneexp[i] = (UWORD(1) << (P->bits*(i+1))) - UWORD(1);                  \
        array_size *= degb;                                                    \
    }                                                                          \
                                                                               \
    off = 0;                                                                   \
    exp = (top << (P->bits*nvars)) + top;                                      \
                                                                               \
    do {                                                                       \
        if (nonzero_test)                                                      \
        {                                                                      \
            _fmpz_mpoly_fit_length(&P->coeffs, &P->exps, &P->alloc, Plen + 1, 1); \
            P->exps[Plen] = exp;                                               \
            swapper                                                            \
            Plen++;                                                            \
        }                                                                      \
                                                                               \
        exp += oneexp[0];                                                      \
        off += 1;                                                              \
        curexp[0] += 1;                                                        \
        if ((exp & mask) == 0)                                                 \
        {                                                                      \
            carry = (nvars - 1 == 0);                                          \
        } else                                                                 \
        {                                                                      \
            carry = 1;                                                         \
            exp -= curexp[0]*oneexp[0];                                        \
            off -= curexp[0];                                                  \
            curexp[0] = 0;                                                     \
            for (i = 1; i < nvars - 1; i++)                                    \
            {                                                                  \
                exp += oneexp[i];                                              \
                off += degpow[i];                                              \
                curexp[i] += 1;                                                \
                if ((exp & mask) == 0)                                         \
                {                                                              \
                    carry = 0;                                                 \
                    break;                                                     \
                } else {                                                       \
                    carry = 1;                                                 \
                    exp -= curexp[i]*oneexp[i];                                \
                    off -= curexp[i]*degpow[i];                                \
                    curexp[i] = 0;                                             \
                }                                                              \
            }                                                                  \
        }                                                                      \
    } while (!carry);                                                          \
                                                                               \
    TMP_END;                                                                   \
                                                                               \
    return Plen;                                                               \
}

DEGREVLEX_UNPACK_MACRO(
    fmpz_mpoly_append_array_sm1_DEGREVLEX, ulong * coeff_array
,
    coeff_array[off] != WORD(0)
,
    fmpz_set_si(P->coeffs + Plen, coeff_array[off]);
    coeff_array[off] = 0;
)

DEGREVLEX_UNPACK_MACRO(
    fmpz_mpoly_append_array_sm2_DEGREVLEX, ulong * coeff_array
,
    (coeff_array[2*off + 0] || coeff_array[2*off + 1]) != WORD(0)
,
    fmpz_set_signed_uiui(P->coeffs + Plen, coeff_array[2*off + 1], coeff_array[2*off + 0]);
    coeff_array[2*off + 0] = coeff_array[2*off + 1] = 0;
)

DEGREVLEX_UNPACK_MACRO(
    fmpz_mpoly_append_array_sm3_DEGREVLEX, ulong * coeff_array
,
    (coeff_array[3*off + 0] || coeff_array[3*off + 1] || coeff_array[3*off + 2]) != WORD(0)
,
    fmpz_set_signed_uiuiui(P->coeffs + Plen, coeff_array[3*off + 2], coeff_array[3*off + 1], coeff_array[3*off + 0]);
    coeff_array[3*off + 0] = coeff_array[3*off + 1] = coeff_array[3*off + 2] = 0;
)

DEGREVLEX_UNPACK_MACRO(
    fmpz_mpoly_append_array_fmpz_DEGREVLEX, fmpz * coeff_array
,
    !fmpz_is_zero(coeff_array + off)
,
    fmpz_swap(P->coeffs + Plen, coeff_array + off);
    fmpz_zero(coeff_array + off);
)



void _fmpz_mpoly_mul_array_chunked_DEG(
    fmpz_mpoly_t P,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    ulong degb,
    const fmpz_mpoly_ctx_t ctx)
{
    slong nvars = ctx->minfo->nvars;
    slong Pi, i, j, Plen, Pl, Al, Bl, array_size;
    slong Abits, * Asum, * Amax, Bbits, * Bsum, * Bmax;
    slong * Amain, * Bmain;
    ulong * Apexp, * Bpexp;
    slong (* upack_sm1)(fmpz_mpoly_t, slong, ulong *, slong, slong, slong); 
    slong (* upack_sm2)(fmpz_mpoly_t, slong, ulong *, slong, slong, slong); 
    slong (* upack_sm3)(fmpz_mpoly_t, slong, ulong *, slong, slong, slong); 
    slong (* upack_fmpz)(fmpz_mpoly_t, slong, fmpz *, slong, slong, slong);
    TMP_INIT;

    TMP_START;

    /* compute lengths of poly2 and poly3 in chunks */
    Al = 1 + (slong) (A->exps[0] >> (A->bits*nvars));
    Bl = 1 + (slong) (B->exps[0] >> (B->bits*nvars));

    array_size = 1;
    for (i = 0; i < nvars-1; i++)
    {
        array_size *= degb;
    }

    upack_sm1  = &fmpz_mpoly_append_array_sm1_DEGLEX;
    upack_sm2  = &fmpz_mpoly_append_array_sm2_DEGLEX;
    upack_sm3  = &fmpz_mpoly_append_array_sm3_DEGLEX;
    upack_fmpz = &fmpz_mpoly_append_array_fmpz_DEGLEX;
    if (ctx->minfo->ord == ORD_DEGREVLEX)
    {
        upack_sm1  = &fmpz_mpoly_append_array_sm1_DEGREVLEX;
        upack_sm2  = &fmpz_mpoly_append_array_sm2_DEGREVLEX;
        upack_sm3  = &fmpz_mpoly_append_array_sm3_DEGREVLEX;
        upack_fmpz = &fmpz_mpoly_append_array_fmpz_DEGREVLEX;
    }

    /* compute indices and lengths of coefficients of polys in main variable */
    Amain = (slong *) TMP_ALLOC(((Al + 1) + Al + Al + (Bl + 1) + Bl + Bl)*sizeof(slong));
    Asum  = Amain + Al + 1;
    Amax  = Asum + Al;
    Bmain = Amax + Al;
    Bsum  = Bmain + Bl + 1;
    Bmax  = Bsum + Bl;
    Apexp = (ulong *) flint_malloc(A->length*sizeof(ulong));
    Bpexp = (ulong *) flint_malloc(B->length*sizeof(ulong));
    mpoly_main_variable_split_DEG(Amain, Apexp, A->exps, Al, A->length,
                                                         degb, nvars, A->bits);
    mpoly_main_variable_split_DEG(Bmain, Bpexp, B->exps, Bl, B->length,
                                                         degb, nvars, B->bits);

    /* work out bit counts for each chunk */
    Abits = 0;
    for (i = 0; i < Al; i++)
    {
        _fmpz_vec_sum_max_bits(&Asum[i], &Amax[i],
                                A->coeffs + Amain[i], Amain[i + 1] - Amain[i]);
        Abits = FLINT_MAX(Abits, Amax[i]);
    }

    Bbits = 0;
    for (j = 0; j < Bl; j++)
    {
        _fmpz_vec_sum_max_bits(&Bsum[j], &Bmax[j],
                                B->coeffs + Bmain[j], Bmain[j + 1] - Bmain[j]);
        Bbits = FLINT_MAX(Bbits, Bmax[j]);
    }

    Pl = Al + Bl - 1;
    FLINT_ASSERT(Pl == degb);
    Plen = 0;

    if (Abits <= (SMALL_FMPZ_BITCOUNT_MAX) && Bbits <= (SMALL_FMPZ_BITCOUNT_MAX))
    {
        ulong * coeff_array = (ulong *) TMP_ALLOC(3*array_size*sizeof(ulong));
        for (j = 0; j < 3*array_size; j++)
            coeff_array[j] = 0;

        /* for each output chunk */
        for (Pi = 0; Pi < Pl; Pi++)
        {
            /* compute bound on coeffs of output chunk */
            slong number = 0;
            slong Pbits = 0;
            for (i = 0, j = Pi; i < Al && j >= 0; i++, j--)
            {
                if (j < Bl)
                {
                    Pbits = FLINT_MAX(Pbits,
                              FLINT_MIN(Asum[i] + Bmax[j], Amax[i] + Bsum[j]));
                    number++;
                }
            }
            Pbits += FLINT_BIT_COUNT(number) + 1; /* includes one bit for sign */

            if (Pbits <= FLINT_BITS)
            {
                /* fits into one word */
                for (i = 0, j = Pi; i < Al && j >= 0; i++, j--)
                {
                    if (j < Bl)
                    {
                        _fmpz_mpoly_addmul_array1_slong1(coeff_array, 
                            (slong *) A->coeffs + Amain[i],
                                Apexp + Amain[i], Amain[i + 1] - Amain[i],
                            (slong *) B->coeffs + Bmain[j],
                                Bpexp + Bmain[j], Bmain[j + 1] - Bmain[j]);
                    }
                }
                Plen = upack_sm1(P, Plen, coeff_array, Pl - Pi - 1, nvars, degb);

            } else if (Pbits <= 2*FLINT_BITS)
            {
                /* fits into two words */
                for (i = 0, j = Pi; i < Al && j >= 0; i++, j--)
                {
                    if (j < Bl)
                    {
                        _fmpz_mpoly_addmul_array1_slong2(coeff_array, 
                            (slong *) A->coeffs + Amain[i],
                                Apexp + Amain[i], Amain[i + 1] - Amain[i],
                            (slong *) B->coeffs + Bmain[j],
                                Bpexp + Bmain[j], Bmain[j + 1] - Bmain[j]);

                    }
                }
                Plen = upack_sm2(P, Plen, coeff_array, Pl - Pi - 1, nvars, degb);

            } else
            {
                /* fits into three words */
                for (i = 0, j = Pi; i < Al && j >= 0; i++, j--)
                {
                    if (j < Bl)
                    {
                        _fmpz_mpoly_addmul_array1_slong(coeff_array, 
                            (slong *) A->coeffs + Amain[i],
                                Apexp + Amain[i], Amain[i + 1] - Amain[i],
                            (slong *) B->coeffs + Bmain[j],
                                Bpexp + Bmain[j], Bmain[j + 1] - Bmain[j]);
                    }
                }
                Plen = upack_sm3(P, Plen, coeff_array, Pl - Pi - 1, nvars, degb);
            }
        }

    } else
    {
        fmpz * coeff_array = (fmpz *) TMP_ALLOC(array_size*sizeof(fmpz));
        for (j = 0; j < array_size; j++)
            fmpz_init(coeff_array + j);

        /* for each output chunk */
        for (Pi = 0; Pi < Pl; Pi++)
        {
            /* addmuls for each cross product of chunks */
            for (i = 0, j = Pi; i < Al && j >= 0; i++, j--)
            {
                if (j < Bl)
                {
                    _fmpz_mpoly_addmul_array1_fmpz(coeff_array, 
                        A->coeffs + Amain[i],
                            Apexp + Amain[i], Amain[i + 1] - Amain[i],
                        B->coeffs + Bmain[j],
                            Bpexp + Bmain[j], Bmain[j + 1] - Bmain[j]);
                }
            }
            Plen = upack_fmpz(P, Plen, coeff_array, Pl - Pi - 1, nvars, degb);
        }
    }

    _fmpz_mpoly_set_length(P, Plen, ctx);

    flint_free(Apexp);
    flint_free(Bpexp);
    TMP_END;
}



int _fmpz_mpoly_mul_array_DEG(
    fmpz_mpoly_t A,
    const fmpz_mpoly_t B, fmpz * maxBfields,
    const fmpz_mpoly_t C, fmpz * maxCfields,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, exp_bits, array_size;
    ulong deg;
    int success;

    FLINT_ASSERT(ctx->minfo->nvars > 0);
    FLINT_ASSERT(B->length != 0);
    FLINT_ASSERT(C->length != 0);

    FLINT_ASSERT(  ctx->minfo->ord == ORD_DEGREVLEX
                || ctx->minfo->ord == ORD_DEGLEX);

    FLINT_ASSERT(1 == mpoly_words_per_exp(B->bits, ctx->minfo));
    FLINT_ASSERT(1 == mpoly_words_per_exp(C->bits, ctx->minfo));

    /* the field of index n-1 is the degree and will be pulled out */
    i = ctx->minfo->nfields - 1;
    FLINT_ASSERT(fmpz_fits_si(maxBfields + i));
    FLINT_ASSERT(fmpz_fits_si(maxCfields + i));
    deg = 1 + fmpz_get_ui(maxBfields + i) + fmpz_get_ui(maxCfields + i);
    if (((slong) deg) <= 0 || deg > MAX_ARRAY_SIZE)
    {
        success = 0;
        goto cleanup;
    }

    /* the fields of index n-2...1, contribute to the array size */
    array_size = WORD(1);
    for (i--; i >= 1; i--)
    {
        ulong hi;
        umul_ppmm(hi, array_size, array_size, deg);
        if (hi != WORD(0) || array_size <= 0
                          || array_size > MAX_ARRAY_SIZE)
        {
            success = 0;
            goto cleanup;
        }
    }

    exp_bits = FLINT_MAX(MPOLY_MIN_BITS, FLINT_BIT_COUNT(deg) + 1);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);

    /* array multiplication assumes result fits into 1 word */
    if (1 != mpoly_words_per_exp(exp_bits, ctx->minfo))
    {
        success = 0;
        goto cleanup;
    }

    /* handle aliasing and do array multiplication */
    if (A == B || A == C)
    {
        fmpz_mpoly_t T;
        fmpz_mpoly_init3(T, B->length + C->length, exp_bits, ctx);
        _fmpz_mpoly_mul_array_chunked_DEG(T, C, B, deg, ctx);
        fmpz_mpoly_swap(T, A, ctx);
        fmpz_mpoly_clear(T, ctx);
    }
    else
    {
        fmpz_mpoly_fit_length_reset_bits(A, B->length + C->length, exp_bits, ctx);
        _fmpz_mpoly_mul_array_chunked_DEG(A, C, B, deg, ctx);
    }
    success = 1;

cleanup:

    return success;
}



int fmpz_mpoly_mul_array(
    fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_t C,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    int success;
    fmpz * maxBfields, * maxCfields;
    TMP_INIT;

    if (B->length == 0 || C->length == 0)
    {
        fmpz_mpoly_zero(A, ctx);
        return 1;
    }

    if (ctx->minfo->nvars < 1 ||
        1 != mpoly_words_per_exp(B->bits, ctx->minfo) ||
        1 != mpoly_words_per_exp(C->bits, ctx->minfo))
    {
        return 0;
    }

    TMP_START;

    maxBfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    maxCfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_init(maxBfields + i);
        fmpz_init(maxCfields + i);
    }
    mpoly_max_fields_fmpz(maxBfields, B->exps, B->length, B->bits, ctx->minfo);
    mpoly_max_fields_fmpz(maxCfields, C->exps, C->length, C->bits, ctx->minfo);

    switch (ctx->minfo->ord)
    {
        case ORD_LEX:
        {
            success = _fmpz_mpoly_mul_array_LEX(A, B, maxBfields,
                                                   C, maxCfields, ctx);
            break;
        }
        case ORD_DEGLEX:
        case ORD_DEGREVLEX:
        {
            success = _fmpz_mpoly_mul_array_DEG(A, B, maxBfields,
                                                   C, maxCfields, ctx);
            break;
        }
        default:
        {
            success = 0;
            break;
        }
    }

    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_clear(maxBfields + i);
        fmpz_clear(maxCfields + i);
    }

    TMP_END;
    return success;
}
