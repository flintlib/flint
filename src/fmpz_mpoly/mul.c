/*
    Copyright (C) 2016 William Hart
    Copyright (C) 2017-2019, 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "long_extras.h"
#include "fmpz_mpoly.h"

static int _try_dense(int try_array, slong * Bdegs, slong * Cdegs,
                                           slong Blen, slong Clen, slong nvars)
{
    const int max_bit_size = FLINT_MIN(FLINT_BITS/3 + 16, FLINT_BITS - 3);
    slong i, product_count, dense_size;
    ulong hi;

    FLINT_ASSERT(Blen > 0);
    FLINT_ASSERT(Clen > 0);

    dense_size = WORD(1);
    for (i = 0; i < nvars; i++)
    {
        umul_ppmm(hi, dense_size, dense_size, Bdegs[i] + Cdegs[i] + 1);

        if (hi != 0 || dense_size <= 0)
            return 0;
    }

    if (dense_size >= WORD(1) << max_bit_size)
        return 0;

    umul_ppmm(hi, product_count, Blen, Clen);

    if (hi != 0 || product_count < 0)
        return 1;

    /*
        Assume that the running time of the dense method is linear
        in "dense_size" and that the running time of the array|heap
        method is linear in "product_count".
        Assume further that the array method is 4x faster than heap.
    */

    if (try_array)
        return dense_size < product_count/128;
    else
        return dense_size < product_count/32;
}


static int _try_array_LEX(slong * Bdegs, slong * Cdegs,
                                           slong Blen, slong Clen, slong nvars)
{
    slong i, dense_size;
    ulong hi;

    FLINT_ASSERT(Blen > 0);
    FLINT_ASSERT(Clen > 0);

    /* accept array method if the array is probably at least 10% full */

    dense_size = WORD(1);
    for (i = 0; i < nvars; i++)
    {
        umul_ppmm(hi, dense_size, dense_size, Bdegs[i] + Cdegs[i] + 1);
        if (hi != 0 || dense_size <= 0)
            return 0;
    }

    return dense_size <= WORD(50000000) &&
           dense_size/Blen/Clen < WORD(10);
}


static int _try_array_DEG(slong Btotaldeg, slong Ctotaldeg,
                                           slong Blen, slong Clen, slong nvars)
{
    slong i, dense_size, total_degree;
    ulong hi;

    FLINT_ASSERT(Blen > 0);
    FLINT_ASSERT(Clen > 0);

    total_degree = Btotaldeg + Btotaldeg;
    if (total_degree <= 0)
        return 0;

    /* the relevant portion of the array has approx size d^nvars/nvars!*/
    dense_size = WORD(1);
    for (i = 0; i < nvars; i++)
    {
        umul_ppmm(hi, dense_size, dense_size, total_degree);
        if (hi != WORD(0) || dense_size < 0)
            return 0;
    }
    for (i = 0; i < nvars; i++)
    {
        dense_size /= i + 1;
    }

    return dense_size <= WORD(5000000) &&
           dense_size/Blen/Clen < WORD(10);
}

/* !!! this function DOES need to change with new orderings */
static int _try_dense_univar(
    fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_t C,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    ulong maskB = (-UWORD(1)) >> (FLINT_BITS - B->bits);
    ulong maskC = (-UWORD(1)) >> (FLINT_BITS - C->bits);
    slong NB = mpoly_words_per_exp_sp(B->bits, ctx->minfo);
    slong NC = mpoly_words_per_exp_sp(C->bits, ctx->minfo);
    slong Blen = B->length;
    slong Clen = C->length;
    slong BClen;
    const ulong * Bexps = B->exps;
    const ulong * Cexps = C->exps;
    fmpz * Acoeffs, * Bcoeffs, * Ccoeffs;
    slong Adeg, Bdeg, Cdeg;
    TMP_INIT;

    /* the variable power is stored at offset 0 */
    Bdeg = Bexps[NB*0 + 0] & maskB;
    Cdeg = Cexps[NC*0 + 0] & maskC;

    if (z_mul_checked(&BClen, Blen, Clen) || z_add_checked(&Adeg, Bdeg, Cdeg))
        return 0;

    if (Adeg > WORD_MAX/FLINT_BITS)
        return 0;

    if (Adeg > BClen)
        return 0;

    {
        slong Bcoeffbits = _fmpz_vec_max_bits(B->coeffs, Blen);
        slong Ccoeffbits = _fmpz_vec_max_bits(C->coeffs, Clen);
        slong t = FLINT_ABS(Bcoeffbits) + FLINT_ABS(Ccoeffbits);

        if (t > FLINT_BITS && Adeg > BClen/4)
            return 0;
    }

    TMP_START;

    Acoeffs = TMP_ARRAY_ALLOC(Adeg + 1 + Bdeg + 1 + Cdeg + 1, fmpz);
    Bcoeffs = Acoeffs + Adeg + 1;
    Ccoeffs = Bcoeffs + Bdeg + 1;

    /* we own the fmpz's in Acoeffs */
    for (i = 0; i < Adeg + 1; i++)
        fmpz_init(Acoeffs + i);

    if (A != B && A != C)
    {
        for (i = FLINT_MIN(A->length - 1, Adeg); i >= 0; i--)
            fmpz_swap(Acoeffs + i, A->coeffs + i);
    }

    /* Bcoeffs and Ccoeffs are shallow copies */
    for (i = 0; i < Bdeg + 1 + Cdeg + 1; i++)
        Bcoeffs[i] = 0;

    for (i = 0; i < Blen; i++)
        Bcoeffs[Bexps[NB*i + 0] & maskB] = B->coeffs[i];

    for (i = 0; i < Clen; i++)
        Ccoeffs[Cexps[NC*i + 0] & maskC] = C->coeffs[i];

    if (Bdeg >= Cdeg)
        _fmpz_poly_mul(Acoeffs, Bcoeffs, Bdeg + 1, Ccoeffs, Cdeg + 1);
    else
        _fmpz_poly_mul(Acoeffs, Ccoeffs, Cdeg + 1, Bcoeffs, Bdeg + 1);

    /* Acoeffs cleared */
    _fmpz_mpoly_set_fmpz_poly_one_var(A, FLINT_MAX(B->bits, C->bits),
                                                           Acoeffs, Adeg, ctx);
    TMP_END;
    return 1;
}

void fmpz_mpoly_mul(
    fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_t C,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    slong nvars = ctx->minfo->nvars;
    int success, try_array;
    slong * Bdegs, * Cdegs;
    fmpz * maxBfields, * maxCfields;
    thread_pool_handle * handles;
    slong num_handles;
    slong min_length, max_length;
    slong thread_limit;
    TMP_INIT;

    if (B->length < 1 || C->length < 1)
    {
        fmpz_mpoly_zero(A, ctx);
        return;
    }
    else if (B->length == 1)
    {
        fmpz_mpoly_mul_monomial(A, C, B, ctx);
        return;
    }
    else if (C->length == 1)
    {
        fmpz_mpoly_mul_monomial(A, B, C, ctx);
        return;
    }

    FLINT_ASSERT(nvars > 0);

    if (nvars == 1 && B->bits <= FLINT_BITS && C->bits <= FLINT_BITS)
    {
        if (_try_dense_univar(A, B, C, ctx))
            return;
    }

    TMP_START;

    /*
        All methods require a linear scan of the exponents.
        Do it here once and for all.
    */
    maxBfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    maxCfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_init(maxBfields + i);
        fmpz_init(maxCfields + i);
    }
    mpoly_max_fields_fmpz(maxBfields, B->exps, B->length, B->bits, ctx->minfo);
    mpoly_max_fields_fmpz(maxCfields, C->exps, C->length, C->bits, ctx->minfo);

    min_length = FLINT_MIN(B->length, C->length);
    max_length = FLINT_MAX(B->length, C->length);
    thread_limit = min_length/512;

    /*
        If one polynomial is tiny or if both polynomials are small,
        heap method with operational complexity O(B->length*C->length) is fine.
    */
    if (min_length < 20 || max_length < 50)
    {
        _fmpz_mpoly_mul_johnson_maxfields(A, B, maxBfields, C, maxCfields, ctx);
        goto cleanup;
    }

    /*
        If either polynomial has multi-word fields, only heap will do.
    */
    if (B->bits > FLINT_BITS || C->bits > FLINT_BITS)
    {
        num_handles = flint_request_threads(&handles, thread_limit);
        goto do_heap;
    }

    /*
        The multiplication is not trivial and each packed field fits
        into one word. In particular, the degrees must fit an slong.
    */
    Bdegs = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    Cdegs = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    mpoly_get_monomial_ui_unpacked_ffmpz((ulong *)Bdegs, maxBfields, ctx->minfo);
    mpoly_get_monomial_ui_unpacked_ffmpz((ulong *)Cdegs, maxCfields, ctx->minfo);

    /*
        See if array method is applicable.
        If so, it should be about 4x faster than heap.
    */
    try_array = 0;
    if (nvars > WORD(1) &&
        nvars < WORD(8) &&
        1 == mpoly_words_per_exp(B->bits, ctx->minfo) &&
        1 == mpoly_words_per_exp(C->bits, ctx->minfo))
    {
        if (ctx->minfo->ord == ORD_LEX)
        {
            try_array = _try_array_LEX(Bdegs, Cdegs, B->length, C->length, nvars);
        }
        else if (ctx->minfo->ord == ORD_DEGLEX || ctx->minfo->ord == ORD_DEGREVLEX)
        {
            slong Btdeg = fmpz_get_si(maxBfields + nvars);
            slong Ctdeg = fmpz_get_si(maxCfields + nvars);
            try_array = _try_array_DEG(Btdeg, Ctdeg, B->length, C->length, nvars);
        }
    }

    success = 0;
    if (_try_dense(try_array, Bdegs, Cdegs, B->length, C->length, nvars))
    {
        success = _fmpz_mpoly_mul_dense(A, B, maxBfields, C, maxCfields, ctx);
        if (success)
        {
            goto cleanup;
        }
    }

    num_handles = flint_request_threads(&handles, thread_limit);

    if (!try_array)
    {
        goto do_heap;
    }

    if (ctx->minfo->ord == ORD_LEX)
    {
        success = (num_handles > 0)
                ? _fmpz_mpoly_mul_array_threaded_pool_LEX(
                                    A, B, maxBfields, C, maxCfields, ctx,
                                                         handles, num_handles)
                : _fmpz_mpoly_mul_array_LEX(
                                    A, B, maxBfields, C, maxCfields, ctx);
    }
    else if (ctx->minfo->ord == ORD_DEGLEX || ctx->minfo->ord == ORD_DEGREVLEX)
    {
        success = (num_handles > 0)
                ? _fmpz_mpoly_mul_array_threaded_pool_DEG(
                                    A, B, maxBfields, C, maxCfields, ctx,
                                                         handles, num_handles)
                : _fmpz_mpoly_mul_array_DEG(
                                    A, B, maxBfields, C, maxCfields, ctx);
    }

    if (success)
    {
        goto cleanup_threads;
    }

do_heap:

    if (num_handles > 0)
    {
        _fmpz_mpoly_mul_heap_threaded_pool_maxfields(A,
                      B, maxBfields, C, maxCfields, ctx, handles, num_handles);
    }
    else
    {
        _fmpz_mpoly_mul_johnson_maxfields(A, B, maxBfields, C, maxCfields, ctx);
    }

cleanup_threads:

    flint_give_back_threads(handles, num_handles);

cleanup:

    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_clear(maxBfields + i);
        fmpz_clear(maxCfields + i);
    }

    TMP_END;
}

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

    /* the field of index n-1 is the one that will be pulled out */
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

#undef BLOCK
#undef MAX_ARRAY_SIZE
#undef MAX_LEX_SIZE

/*
    NOTE: this file is dirty - it assumes that a zero fmpz is zero
*/

/* improve locality */
#define BLOCK 128
#define MAX_ARRAY_SIZE (WORD(300000))
#define MAX_LEX_SIZE (WORD(300))


typedef struct
{
    slong idx;
    slong work;
    slong len;
    fmpz_mpoly_t poly;
}
_chunk_struct;


typedef struct
{
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
#endif
    volatile int idx;
    slong nthreads;
    slong Al, Bl, Pl;
    fmpz * Acoeffs, * Bcoeffs;
    slong * Amax, * Bmax, * Asum, * Bsum;
    slong * Amain, * Bmain;
    ulong * Apexp, * Bpexp;
    slong * perm;
    slong nvars;
    const ulong * mults;
    slong array_size;
    slong degb;
    _chunk_struct * Pchunks;
    int rev;
}
_base_struct1;

typedef _base_struct1 _base_t1[1];


typedef struct
{
    slong idx;
    slong time;
    _base_struct1 * base;
    ulong * exp;
}
_worker_arg_struct1;



/******************
    LEX
******************/

void _fmpz_mpoly_mul_array_threaded_worker_LEX(void * varg)
{
    slong i, j, Pi;
    _worker_arg_struct1 * arg = (_worker_arg_struct1 *) varg;
    _base_struct1 * base = arg->base;
    ulong * coeff_array;
    TMP_INIT;

    TMP_START;
    coeff_array = (ulong *) TMP_ALLOC(3*base->array_size*sizeof(ulong));
    for (j = 0; j < 3*base->array_size; j++)
        coeff_array[j] = 0;

#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&base->mutex);
#endif
    Pi = base->idx;
    base->idx = Pi + 1;
#if FLINT_USES_PTHREAD
    pthread_mutex_unlock(&base->mutex);
#endif

    while (Pi < base->Pl)
    {
        /* work out bit counts for this chunk */
        slong Abits = 0;
        slong Bbits = 0;
        slong Pbits = 0;
        slong number = 0;
        for (i = 0, j = base->perm[Pi]; i < base->Al && j >= 0; i++, j--)
        {
            if (j < base->Bl)
            {
                number++;
                Abits = FLINT_MAX(Abits, base->Amax[i]);
                Bbits = FLINT_MAX(Bbits, base->Bmax[j]);
                Pbits = FLINT_MAX(Pbits,
                            FLINT_MIN(base->Asum[i] + base->Bmax[j],
                                      base->Amax[i] + base->Bsum[j]));
            }
        }
        Pbits += FLINT_BIT_COUNT(number) + 1; /* includes one bit for sign */

        if (Abits <= SMALL_FMPZ_BITCOUNT_MAX && Bbits <= SMALL_FMPZ_BITCOUNT_MAX)
        {
            if (Pbits <= FLINT_BITS)
            {
                for (i = 0, j = base->perm[Pi]; i < base->Al && j >= 0; i++, j--)
                {
                    if (j < base->Bl)
                    {
                        _fmpz_mpoly_addmul_array1_slong1((ulong *)coeff_array, 
                            base->Acoeffs + base->Amain[i],
                                base->Apexp + base->Amain[i],
                                base->Amain[i + 1] - base->Amain[i],
                            base->Bcoeffs + base->Bmain[j],
                                base->Bpexp + base->Bmain[j],
                                base->Bmain[j + 1] - base->Bmain[j]);
                    }
                }
                (base->Pchunks + base->perm[Pi])->len = 
                    fmpz_mpoly_append_array_sm1_LEX(
                        (base->Pchunks + base->perm[Pi])->poly, 0,
                        (ulong *)coeff_array, base->mults, base->nvars - 1,
                        base->array_size, base->Pl - base->perm[Pi] - 1);

            } else if (Pbits <= 2*FLINT_BITS)
            {
                for (i = 0, j = base->perm[Pi]; i < base->Al && j >= 0; i++, j--)
                {
                    if (j < base->Bl)
                    {
                        _fmpz_mpoly_addmul_array1_slong2((ulong *)coeff_array, 
                            base->Acoeffs + base->Amain[i],
                                base->Apexp + base->Amain[i],
                                base->Amain[i + 1] - base->Amain[i],
                            base->Bcoeffs + base->Bmain[j],
                                base->Bpexp + base->Bmain[j],
                                base->Bmain[j + 1] - base->Bmain[j]);
                    }
                }
                (base->Pchunks + base->perm[Pi])->len = 
                    fmpz_mpoly_append_array_sm2_LEX(
                        (base->Pchunks + base->perm[Pi])->poly, 0,
                        (ulong *)coeff_array, base->mults, base->nvars - 1,
                        base->array_size, base->Pl - base->perm[Pi] - 1);
            } else
            {
                for (i = 0, j = base->perm[Pi]; i < base->Al && j >= 0; i++, j--)
                {
                    if (j < base->Bl)
                    {
                        _fmpz_mpoly_addmul_array1_slong((ulong *)coeff_array, 
                            base->Acoeffs + base->Amain[i],
                                base->Apexp + base->Amain[i],
                                base->Amain[i + 1] - base->Amain[i],
                            base->Bcoeffs + base->Bmain[j],
                                base->Bpexp + base->Bmain[j],
                                base->Bmain[j + 1] - base->Bmain[j]);
                    }
                }
                (base->Pchunks + base->perm[Pi])->len = 
                    fmpz_mpoly_append_array_sm3_LEX(
                        (base->Pchunks + base->perm[Pi])->poly, 0,
                        (ulong *)coeff_array, base->mults, base->nvars - 1,
                        base->array_size, base->Pl - base->perm[Pi] - 1);
            }
        } else
        {
            for (i = 0, j = base->perm[Pi]; i < base->Al && j >= 0; i++, j--)
            {
                if (j < base->Bl)
                {
                    _fmpz_mpoly_addmul_array1_fmpz((fmpz *)coeff_array, 
                        base->Acoeffs + base->Amain[i],
                            base->Apexp + base->Amain[i],
                            base->Amain[i + 1] - base->Amain[i],
                        base->Bcoeffs + base->Bmain[j],
                            base->Bpexp + base->Bmain[j],
                            base->Bmain[j + 1] - base->Bmain[j]);
                }
            }
            (base->Pchunks + base->perm[Pi])->len = 
                fmpz_mpoly_append_array_fmpz_LEX(
                    (base->Pchunks + base->perm[Pi])->poly, 0,
                    (fmpz *)coeff_array, base->mults, base->nvars - 1,
                    base->array_size, base->Pl - base->perm[Pi] - 1);
        }

#if FLINT_USES_PTHREAD
        pthread_mutex_lock(&base->mutex);
#endif
        Pi = base->idx;
        base->idx = Pi + 1;
#if FLINT_USES_PTHREAD
        pthread_mutex_unlock(&base->mutex);
#endif
    }

    TMP_END;
}


void _fmpz_mpoly_mul_array_chunked_threaded_LEX(
    fmpz_mpoly_t P,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const ulong * mults,
    const fmpz_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles)
{
    slong nvars = ctx->minfo->nvars;
    slong Pi, i, j, Plen, Pl, Al, Bl, array_size;
    slong * Asum, * Amax, * Bsum, * Bmax;
    slong * Amain, * Bmain;
    ulong * Apexp, * Bpexp;
    _base_t1 base;
    _worker_arg_struct1 * args;
    _chunk_struct * Pchunks;
    slong * perm;
    TMP_INIT;

    array_size = 1;
    for (i = 0; i < nvars - 1; i++) {
        array_size *= mults[i];
    }

    /* compute lengths of poly2 and poly3 in chunks */
    Al = 1 + (slong) (A->exps[0] >> (A->bits*(nvars - 1)));
    Bl = 1 + (slong) (B->exps[0] >> (B->bits*(nvars - 1)));

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
    mpoly_main_variable_split_LEX(Amain, Apexp, A->exps, Al, A->length,
                                                    mults, nvars - 1, A->bits);
    mpoly_main_variable_split_LEX(Bmain, Bpexp, B->exps, Bl, B->length,
                                                    mults, nvars - 1, B->bits);

    /* work out bit counts for each chunk */
    for (i = 0; i < Al; i++)
    {
        _fmpz_vec_sum_max_bits(&Asum[i], &Amax[i],
                                A->coeffs + Amain[i], Amain[i + 1] - Amain[i]);
    }
    for (j = 0; j < Bl; j++)
    {
        _fmpz_vec_sum_max_bits(&Bsum[j], &Bmax[j],
                                B->coeffs + Bmain[j], Bmain[j + 1] - Bmain[j]);
    }

    Pl = Al + Bl - 1;

    /* work out data for each chunk of the output */
    Pchunks = (_chunk_struct *) TMP_ALLOC(Pl*sizeof(_chunk_struct));
    perm = (slong *) TMP_ALLOC(Pl*sizeof(slong));
    for (Pi = 0; Pi < Pl; Pi++)
    {
        fmpz_mpoly_init3((Pchunks + Pi)->poly, 8, P->bits, ctx);
        (Pchunks + Pi)->work = 0;
        perm[Pi] = Pi;
        for (i = 0, j = Pi; i < Al && j >= 0; i++, j--)
        {
            if (j < Bl)
            {
                (Pchunks + Pi)->work += (Amain[i + 1] - Amain[i])
                                       *(Bmain[j + 1] - Bmain[j]);
            }
        }
    }

    for (i = 0; i < Pl; i++)
    {
        for (j = i; j > 0 && (Pchunks + perm[j - 1])->work
                                   < (Pchunks + perm[j])->work; j--)
        {
            slong t = perm[j - 1];
            perm[j - 1] = perm[j];
            perm[j] = t;
        }
    }

    base->nthreads = num_handles + 1;
    base->Al = Al;
    base->Bl = Bl;
    base->Pl = Pl;
    base->Amax = Amax;
    base->Bmax = Bmax;
    base->Asum = Asum;
    base->Bsum = Bsum;
    base->Acoeffs = A->coeffs;
    base->Amain = Amain;
    base->Apexp = Apexp;
    base->Bcoeffs = B->coeffs;
    base->Bmain = Bmain;
    base->Bpexp = Bpexp;
    base->idx = 0;
    base->perm = perm;
    base->nvars = nvars;
    base->Pchunks = Pchunks;
    base->array_size = array_size;
    base->mults = mults;

    args = (_worker_arg_struct1 *) TMP_ALLOC(base->nthreads
                                                  *sizeof(_worker_arg_struct1));

#if FLINT_USES_PTHREAD
    pthread_mutex_init(&base->mutex, NULL);
#endif
    for (i = 0; i < num_handles; i++)
    {
        args[i].idx = i;
        args[i].base = base;
        thread_pool_wake(global_thread_pool, handles[i], 0,
                          _fmpz_mpoly_mul_array_threaded_worker_LEX, &args[i]);
    }
    i = num_handles;
    args[i].idx = i;
    args[i].base = base;
    _fmpz_mpoly_mul_array_threaded_worker_LEX(&args[i]);
    for (i = 0; i < num_handles; i++)
    {
        thread_pool_wait(global_thread_pool, handles[i]);
    }
#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&base->mutex);
#endif

    /* join answers */
    Plen = 0;
    for (Pi = 0; Pi < Pl; Pi++)
    {
        _fmpz_mpoly_fit_length(&P->coeffs, &P->exps, &P->alloc,
                                                Plen + (Pchunks + Pi)->len, 1);
        for (i = 0; i < (Pchunks + Pi)->len; i++)
        {
            P->exps[Plen] = (Pchunks + Pi)->poly->exps[i];
            fmpz_swap(P->coeffs + Plen, (Pchunks + Pi)->poly->coeffs + i);
            fmpz_clear((Pchunks + Pi)->poly->coeffs + i);
            Plen++;
        }

        flint_free((Pchunks + Pi)->poly->coeffs);
        flint_free((Pchunks + Pi)->poly->exps);
    }

    _fmpz_mpoly_set_length(P, Plen, ctx);

    flint_free(Apexp);
    flint_free(Bpexp);
    TMP_END;
}


int _fmpz_mpoly_mul_array_threaded_pool_LEX(
    fmpz_mpoly_t A,
    const fmpz_mpoly_t B, fmpz * maxBfields,
    const fmpz_mpoly_t C, fmpz * maxCfields,
    const fmpz_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles)
{
    slong i, exp_bits, array_size;
    ulong max, * mults;
    int success;
    TMP_INIT;

    FLINT_ASSERT(B->length != 0);
    FLINT_ASSERT(C->length != 0);

    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    FLINT_ASSERT(1 == mpoly_words_per_exp(B->bits, ctx->minfo));
    FLINT_ASSERT(1 == mpoly_words_per_exp(C->bits, ctx->minfo));

    TMP_START;

    /* compute maximum exponents for each variable */
    mults = (ulong *) TMP_ALLOC(ctx->minfo->nfields*sizeof(ulong));

    /* the field of index n-1 is the one that will be pulled out */
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
        if (hi != WORD(0) || (array_size | (slong) mults[i]) <= 0
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
        _fmpz_mpoly_mul_array_chunked_threaded_LEX(T, C, B, mults, ctx,
                                                         handles, num_handles);
        fmpz_mpoly_swap(T, A, ctx);
        fmpz_mpoly_clear(T, ctx);
    }
    else
    {
        fmpz_mpoly_fit_length_reset_bits(A, B->length + C->length, exp_bits, ctx);
        _fmpz_mpoly_mul_array_chunked_threaded_LEX(A, C, B, mults, ctx,
                                                         handles, num_handles);
    }
    success = 1;

cleanup:

    TMP_END;

    return success;
}




/*****************************
    DEGLEX and DEGREVLEX
*****************************/


void _fmpz_mpoly_mul_array_threaded_worker_DEG(void * varg)
{
    slong i, j, Pi;
    _worker_arg_struct1 * arg = (_worker_arg_struct1 *) varg;
    _base_struct1 * base = arg->base;
    ulong * coeff_array;
    slong (* upack_sm1)(fmpz_mpoly_t, slong, ulong *, slong, slong, slong); 
    slong (* upack_sm2)(fmpz_mpoly_t, slong, ulong *, slong, slong, slong); 
    slong (* upack_sm3)(fmpz_mpoly_t, slong, ulong *, slong, slong, slong); 
    slong (* upack_fmpz)(fmpz_mpoly_t, slong, fmpz *, slong, slong, slong); 
    TMP_INIT;

    upack_sm1  = &fmpz_mpoly_append_array_sm1_DEGLEX;
    upack_sm2  = &fmpz_mpoly_append_array_sm2_DEGLEX;
    upack_sm3  = &fmpz_mpoly_append_array_sm3_DEGLEX;
    upack_fmpz = &fmpz_mpoly_append_array_fmpz_DEGLEX;
    if (base->rev)
    {
        upack_sm1  = &fmpz_mpoly_append_array_sm1_DEGREVLEX;
        upack_sm2  = &fmpz_mpoly_append_array_sm2_DEGREVLEX;
        upack_sm3  = &fmpz_mpoly_append_array_sm3_DEGREVLEX;
        upack_fmpz = &fmpz_mpoly_append_array_fmpz_DEGREVLEX;
    }

    TMP_START;
    coeff_array = (ulong *) TMP_ALLOC(3*base->array_size*sizeof(ulong));
    for (j = 0; j < 3*base->array_size; j++)
        coeff_array[j] = 0;

#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&base->mutex);
#endif
    Pi = base->idx;
    base->idx = Pi + 1;
#if FLINT_USES_PTHREAD
    pthread_mutex_unlock(&base->mutex);
#endif

    while (Pi < base->Pl)
    {
        /* work out bit counts for this chunk */
        slong Abits = 0;
        slong Bbits = 0;
        slong Pbits = 0;
        slong number = 0;
        for (i = 0, j = base->perm[Pi]; i < base->Al && j >= 0; i++, j--)
        {
            if (j < base->Bl)
            {
                number++;
                Abits = FLINT_MAX(Abits, base->Amax[i]);
                Bbits = FLINT_MAX(Bbits, base->Bmax[j]);
                Pbits = FLINT_MAX(Pbits,
                            FLINT_MIN(base->Asum[i] + base->Bmax[j],
                                      base->Amax[i] + base->Bsum[j]));
            }
        }
        Pbits += FLINT_BIT_COUNT(number) + 1; /* includes one bit for sign */

        if (Abits <= SMALL_FMPZ_BITCOUNT_MAX && Bbits <= SMALL_FMPZ_BITCOUNT_MAX)
        {
            if (Pbits <= FLINT_BITS)
            {
                for (i = 0, j = base->perm[Pi]; i < base->Al && j >= 0; i++, j--)
                {
                    if (j < base->Bl)
                    {
                        _fmpz_mpoly_addmul_array1_slong1((ulong *)coeff_array, 
                            base->Acoeffs + base->Amain[i],
                                base->Apexp + base->Amain[i],
                                base->Amain[i + 1] - base->Amain[i],
                            base->Bcoeffs + base->Bmain[j],
                                base->Bpexp + base->Bmain[j],
                                base->Bmain[j + 1] - base->Bmain[j]);
                    }
                }
                (base->Pchunks + base->perm[Pi])->len = 
                    upack_sm1((base->Pchunks + base->perm[Pi])->poly, 0,
                        (ulong *)coeff_array, base->Pl - base->perm[Pi] - 1,
                                                      base->nvars, base->degb);

            } else if (Pbits <= 2*FLINT_BITS)
            {
                for (i = 0, j = base->perm[Pi]; i < base->Al && j >= 0; i++, j--)
                {
                    if (j < base->Bl)
                    {
                        _fmpz_mpoly_addmul_array1_slong2((ulong *)coeff_array, 
                            base->Acoeffs + base->Amain[i],
                                base->Apexp + base->Amain[i],
                                base->Amain[i + 1] - base->Amain[i],
                            base->Bcoeffs + base->Bmain[j],
                                base->Bpexp + base->Bmain[j],
                                base->Bmain[j + 1] - base->Bmain[j]);
                    }
                }
                (base->Pchunks + base->perm[Pi])->len = 
                    upack_sm2((base->Pchunks + base->perm[Pi])->poly, 0,
                        (ulong *)coeff_array, base->Pl - base->perm[Pi] - 1,
                                                      base->nvars, base->degb);
            } else
            {
                for (i = 0, j = base->perm[Pi]; i < base->Al && j >= 0; i++, j--)
                {
                    if (j < base->Bl)
                    {
                        _fmpz_mpoly_addmul_array1_slong((ulong *)coeff_array, 
                            base->Acoeffs + base->Amain[i],
                                base->Apexp + base->Amain[i],
                                base->Amain[i + 1] - base->Amain[i],
                            base->Bcoeffs + base->Bmain[j],
                                base->Bpexp + base->Bmain[j],
                                base->Bmain[j + 1] - base->Bmain[j]);
                    }
                }
                (base->Pchunks + base->perm[Pi])->len = 
                    upack_sm3((base->Pchunks + base->perm[Pi])->poly, 0,
                        (ulong *)coeff_array, base->Pl - base->perm[Pi] - 1,
                                                      base->nvars, base->degb);
            }
        } else
        {
            for (i = 0, j = base->perm[Pi]; i < base->Al && j >= 0; i++, j--)
            {
                if (j < base->Bl)
                {
                    _fmpz_mpoly_addmul_array1_fmpz((fmpz *)coeff_array, 
                        base->Acoeffs + base->Amain[i],
                            base->Apexp + base->Amain[i],
                            base->Amain[i + 1] - base->Amain[i],
                        base->Bcoeffs + base->Bmain[j],
                            base->Bpexp + base->Bmain[j],
                            base->Bmain[j + 1] - base->Bmain[j]);
                }
            }
            (base->Pchunks + base->perm[Pi])->len = 
                upack_fmpz((base->Pchunks + base->perm[Pi])->poly, 0,
                    (fmpz *)coeff_array, base->Pl - base->perm[Pi] - 1,
                                                      base->nvars, base->degb);
        }

#if FLINT_USES_PTHREAD
        pthread_mutex_lock(&base->mutex);
#endif
        Pi = base->idx;
        base->idx = Pi + 1;
#if FLINT_USES_PTHREAD
        pthread_mutex_unlock(&base->mutex);
#endif
    }

    TMP_END;
}



void _fmpz_mpoly_mul_array_chunked_threaded_DEG(
    fmpz_mpoly_t P,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    ulong degb,
    const fmpz_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles)
{
    slong nvars = ctx->minfo->nvars;
    slong Pi, i, j, Plen, Pl, Al, Bl, array_size;
    slong * Asum, * Amax, * Bsum, * Bmax;
    slong * Amain, * Bmain;
    ulong * Apexp, * Bpexp;
    _base_t1 base;
    _worker_arg_struct1 * args;
    _chunk_struct * Pchunks;
    slong * perm;
    TMP_INIT;

    /* compute lengths of poly2 and poly3 in chunks */
    Al = 1 + (slong) (A->exps[0] >> (A->bits*nvars));
    Bl = 1 + (slong) (B->exps[0] >> (B->bits*nvars));

    array_size = 1;
    for (i = 0; i < nvars-1; i++) {
        array_size *= degb;
    }

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
    mpoly_main_variable_split_DEG(Amain, Apexp, A->exps, Al, A->length,
                                                         degb, nvars, A->bits);
    mpoly_main_variable_split_DEG(Bmain, Bpexp, B->exps, Bl, B->length,
                                                         degb, nvars, B->bits);

    /* work out bit counts for each chunk */
    for (i = 0; i < Al; i++)
    {
        _fmpz_vec_sum_max_bits(&Asum[i], &Amax[i],
                                A->coeffs + Amain[i], Amain[i + 1] - Amain[i]);
    }
    for (j = 0; j < Bl; j++)
    {
        _fmpz_vec_sum_max_bits(&Bsum[j], &Bmax[j],
                                B->coeffs + Bmain[j], Bmain[j + 1] - Bmain[j]);
    }

    Pl = Al + Bl - 1;
    FLINT_ASSERT(Pl == degb);

    /* work out data for each chunk of the output */
    Pchunks = (_chunk_struct *) TMP_ALLOC(Pl*sizeof(_chunk_struct));
    perm = (slong *) TMP_ALLOC(Pl*sizeof(slong));
    for (Pi = 0; Pi < Pl; Pi++)
    {
        fmpz_mpoly_init3((Pchunks + Pi)->poly, 8, P->bits, ctx);
        (Pchunks + Pi)->work = 0;
        perm[Pi] = Pi;
        for (i = 0, j = Pi; i < Al && j >= 0; i++, j--)
        {
            if (j < Bl)
            {
                (Pchunks + Pi)->work += (Amain[i + 1] - Amain[i])
                                       *(Bmain[j + 1] - Bmain[j]);
            }
        }
    }

    for (i = 0; i < Pl; i++)
    {
        for (j = i; j > 0 && (Pchunks + perm[j-1])->work
                                         < (Pchunks + perm[j])->work; j--)
        {
            slong t = perm[j - 1];
            perm[j - 1] = perm[j];
            perm[j] = t;
        }
    }

    base->nthreads = num_handles + 1;
    base->Al = Al;
    base->Bl = Bl;
    base->Pl = Pl;
    base->Amax = Amax;
    base->Bmax = Bmax;
    base->Asum = Asum;
    base->Bsum = Bsum;
    base->Acoeffs = A->coeffs;
    base->Amain = Amain;
    base->Apexp = Apexp;
    base->Bcoeffs = B->coeffs;
    base->Bmain = Bmain;
    base->Bpexp = Bpexp;
    base->idx = 0;
    base->perm = perm;
    base->nvars = nvars;
    base->Pchunks = Pchunks;
    base->array_size = array_size;
    base->degb = degb;
    base->rev = (ctx->minfo->ord == ORD_DEGREVLEX);

    args = (_worker_arg_struct1 *) TMP_ALLOC(base->nthreads
                                                  *sizeof(_worker_arg_struct1));

#if FLINT_USES_PTHREAD
    pthread_mutex_init(&base->mutex, NULL);
#endif
    for (i = 0; i < num_handles; i++)
    {
        args[i].idx = i;
        args[i].base = base;

        thread_pool_wake(global_thread_pool, handles[i], 0,
                          _fmpz_mpoly_mul_array_threaded_worker_DEG, &args[i]);
    }
    i = num_handles;
    args[i].idx = i;
    args[i].base = base;
    _fmpz_mpoly_mul_array_threaded_worker_DEG(&args[i]);
    for (i = 0; i < num_handles; i++)
    {
        thread_pool_wait(global_thread_pool, handles[i]);
    }
#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&base->mutex);
#endif

    /* join answers */
    Plen = 0;
    for (Pi = 0; Pi < Pl; Pi++)
    {
        _fmpz_mpoly_fit_length(&P->coeffs, &P->exps, &P->alloc,
                                                Plen + (Pchunks + Pi)->len, 1);
        for (i = 0; i < (Pchunks + Pi)->len; i++)
        {
            P->exps[Plen] = (Pchunks + Pi)->poly->exps[i];
            fmpz_swap(P->coeffs + Plen, (Pchunks + Pi)->poly->coeffs + i);
            fmpz_clear((Pchunks + Pi)->poly->coeffs + i);
            Plen++;
        }

        flint_free((Pchunks + Pi)->poly->coeffs);
        flint_free((Pchunks + Pi)->poly->exps);
    }

    _fmpz_mpoly_set_length(P, Plen, ctx);

    flint_free(Apexp);
    flint_free(Bpexp);
    TMP_END;
}


int _fmpz_mpoly_mul_array_threaded_pool_DEG(
    fmpz_mpoly_t A,
    const fmpz_mpoly_t B, fmpz * maxBfields,
    const fmpz_mpoly_t C, fmpz * maxCfields,
    const fmpz_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles)
{
    slong i, exp_bits, array_size;
    ulong deg;
    int success;

    FLINT_ASSERT(B->length != 0);
    FLINT_ASSERT(C->length != 0);

    FLINT_ASSERT(  ctx->minfo->ord == ORD_DEGREVLEX
                || ctx->minfo->ord == ORD_DEGLEX);

    FLINT_ASSERT(1 == mpoly_words_per_exp(B->bits, ctx->minfo));
    FLINT_ASSERT(1 == mpoly_words_per_exp(C->bits, ctx->minfo));

    /* the field of index n-1 is the one that will be pulled out */
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

    /* array multiplication assumes result fit into 1 word */
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
        _fmpz_mpoly_mul_array_chunked_threaded_DEG(T, C, B, deg, ctx,
                                                         handles, num_handles);
        fmpz_mpoly_swap(T, A, ctx);
        fmpz_mpoly_clear(T, ctx);
    }
    else
    {
        fmpz_mpoly_fit_length_reset_bits(A, B->length + C->length, exp_bits, ctx);
        _fmpz_mpoly_mul_array_chunked_threaded_DEG(A, C, B, deg, ctx,
                                                         handles, num_handles);
    }
    success = 1;

cleanup:

    return success;
}


int fmpz_mpoly_mul_array_threaded(
    fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_t C,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    int success;
    fmpz * maxBfields, * maxCfields;
    thread_pool_handle * handles;
    slong num_handles;
    slong thread_limit = FLINT_MIN(A->length, B->length)/16;
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

    num_handles = flint_request_threads(&handles, thread_limit);

    switch (ctx->minfo->ord)
    {
        case ORD_LEX:
        {
            success = _fmpz_mpoly_mul_array_threaded_pool_LEX(A,
                      B, maxBfields, C, maxCfields, ctx, handles, num_handles);
            break;
        }
        case ORD_DEGREVLEX:
        case ORD_DEGLEX:
        {
            success = _fmpz_mpoly_mul_array_threaded_pool_DEG(A,
                      B, maxBfields, C, maxCfields, ctx, handles, num_handles);
            break;
        }
        default:
        {
            success = 0;
            break;
        }
    }

    flint_give_back_threads(handles, num_handles);

    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_clear(maxBfields + i);
        fmpz_clear(maxCfields + i);
    }

    TMP_END;
    return success;
}

#undef BLOCK
#undef MAX_ARRAY_SIZE
#undef MAX_LEX_SIZE

slong fmpz_mpolyd_length(const fmpz_mpolyd_t A)
{
    slong i, j;
    slong degb_prod;

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++)
    {
        degb_prod *= A->deg_bounds[j];
    }

    for (i = degb_prod; i > 0; i--)
    {
        if (!fmpz_is_zero(A->coeffs + i - 1))
        {
            break;
        }
    }

    return i;
}


int fmpz_mpolyd_set_degbounds(fmpz_mpolyd_t A, slong * bounds)
{
    slong i;
    int success = 0;
    slong degb_prod;

    degb_prod = 1;
    for (i = 0; i < A->nvars; i++)
    {
        ulong hi;
        A->deg_bounds[i] = bounds[i];
        umul_ppmm(hi, degb_prod, degb_prod, A->deg_bounds[i]);
        if (hi != WORD(0) || degb_prod < 0)
        {
            goto done;
        }
    }

    success = 1;
    fmpz_mpolyd_fit_length(A, degb_prod);

done:
    return success;
}


/*
    assuming poly1 has valid degree bounds set, pack poly2 into it.
*/
void fmpz_mpoly_convert_to_fmpz_mpolyd_degbound(fmpz_mpolyd_t poly1,
                          const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx)
{
    slong degb_prod;
    slong i, j, N;
    ulong * exps;
    slong nvars = ctx->minfo->nvars;
    TMP_INIT;

    FLINT_ASSERT(poly1->nvars == nvars);
    FLINT_ASSERT(poly2->bits <= FLINT_BITS);

    degb_prod = WORD(1);
    for (i = 0; i < nvars; i++)
    {
        degb_prod *= poly1->deg_bounds[i];
    }

    for (i = 0; i < degb_prod; i++)
    {
        fmpz_zero(poly1->coeffs + i);
    }

    if (poly2->length == 0)
    {
        return;
    }

    TMP_START;
    exps = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));
    N = mpoly_words_per_exp(poly2->bits, ctx->minfo);
    for (i = 0; i < poly2->length; i++)
    {
        slong off;

        mpoly_get_monomial_ui(exps, poly2->exps + N*i, poly2->bits, ctx->minfo);
        off = 0;
        for (j = 0; j < nvars; j++)
        {
            off = exps[j] + poly1->deg_bounds[j]*off;
        }
        fmpz_set(poly1->coeffs + off, poly2->coeffs + i);
    }

    TMP_END;
}


/*
    Convert B to A and clear B in the process
*/
void fmpz_mpoly_consume_fmpz_mpolyd_clear(fmpz_mpoly_t A, fmpz_mpolyd_t B,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j, k, N;
    slong bits, nvars = ctx->minfo->nvars;
    slong Alen;
    ulong diff;
    ulong topmask;
    ulong * exps, * ptempexp, * plastexp;
    TMP_INIT;

    FLINT_ASSERT(nvars == B->nvars);

    TMP_START;
    exps = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));

    /* clear all irrelevant coefficients */
    for (i = B->coeff_alloc - 1; i >= B->length; i--)
        fmpz_clear(B->coeffs + i);

    /* find bits needed for the result */
    for (j = 0; j < nvars; j++)
    {
        exps[j] = B->deg_bounds[j] - 1;
    }
    bits = mpoly_exp_bits_required_ui(exps, ctx->minfo);
    bits = mpoly_fix_bits(bits, ctx->minfo);
    N = mpoly_words_per_exp(bits, ctx->minfo);

    /* we are going to push back terms manually */
    Alen = 0;
    fmpz_mpoly_zero(A, ctx);
    fmpz_mpoly_fit_length_reset_bits(A, 0, bits, ctx);

    /* find exponent vector for least significant variable */
    plastexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    for (j = 0; j < nvars; j++)
        exps[j] = (j == nvars - 1);
    mpoly_set_monomial_ui(plastexp, exps, bits, ctx->minfo);

    /* get most significant exponent in exps and its vector in ptempexp */
    ptempexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    k = i;
    for (j = nvars - 1; j >= 0; j--) 
    {
        exps[j] = k % B->deg_bounds[j];
        k = k / B->deg_bounds[j];
    }
    mpoly_set_monomial_ui(ptempexp, exps, bits, ctx->minfo);
    diff = 0;

    /* scan down through the exponents */
    topmask = 0;
    for (; i >= 0; i--)
    {
        if (!fmpz_is_zero(B->coeffs + i))
        {
            _fmpz_mpoly_fit_length(&A->coeffs, &A->exps, &A->alloc, Alen + 1, N);
            fmpz_swap(A->coeffs + Alen, B->coeffs + i);
            mpoly_monomial_msub_mp(A->exps + N*Alen, ptempexp, diff, plastexp, N);
            topmask |= (A->exps + N*Alen)[N - 1];
            Alen++;
        }
        fmpz_clear(B->coeffs + i);

        diff++;
        --exps[nvars - 1];
        if ((slong)(exps[nvars - 1]) < WORD(0))
        {
            exps[nvars - 1] = B->deg_bounds[nvars - 1] - 1;            
            for (j = nvars - 2; j >= 0; j--)
            {
                --exps[j];
                if ((slong)(exps[j]) < WORD(0))
                {
                    FLINT_ASSERT(i == 0 || j > 0);
                    exps[j] = B->deg_bounds[j] - 1;
                } else
                {
                    break;
                }
            }

            mpoly_set_monomial_ui(ptempexp, exps, bits, ctx->minfo);
            diff = 0;
        }
    }
    _fmpz_mpoly_set_length(A, Alen, ctx);

    /* sort the exponents if needed */
    if (ctx->minfo->ord != ORD_LEX)
    {
        slong msb;
        mpoly_get_cmpmask(ptempexp, N, bits, ctx->minfo);
        if (topmask != WORD(0))
        {
            count_leading_zeros(msb, topmask);
            msb = (FLINT_BITS - 1)^msb;
        } else
        {
            msb = -WORD(1);
        }
        if (N == 1) {
            if (msb >= WORD(0))
            {
                _fmpz_mpoly_radix_sort1(A, 0, A->length,
                                                   msb, ptempexp[0], topmask);
            }
        } else {
            _fmpz_mpoly_radix_sort(A, 0, A->length,
                                        (N - 1)*FLINT_BITS + msb, N, ptempexp);
        }
    }

    flint_free(B->deg_bounds);
    flint_free(B->coeffs);
    B->deg_bounds = NULL;
    B->coeffs = NULL;
    TMP_END;
}



int _fmpz_mpoly_mul_dense(fmpz_mpoly_t P,
                                 const fmpz_mpoly_t A, fmpz * maxAfields,
                                 const fmpz_mpoly_t B, fmpz * maxBfields,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    int success, P_is_stolen;
    slong i;
    slong nvars = ctx->minfo->nvars;
    fmpz_mpolyd_t Ad, Bd, Pd;
    fmpz_poly_t Au, Bu, Pu;
    slong * Abounds, * Bbounds, * Pbounds;
    TMP_INIT;

    FLINT_ASSERT(A->length != 0);
    FLINT_ASSERT(B->length != 0);
    FLINT_ASSERT(nvars > 0);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);

    TMP_START;
    /*
        for each variable v except for the outer variable,
        we need to pack to degree deg_v(A) + deg_v(A)
    */
    Abounds = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    Bbounds = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    Pbounds = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    mpoly_get_monomial_ui_unpacked_ffmpz((ulong *)Abounds, maxAfields, ctx->minfo);
    mpoly_get_monomial_ui_unpacked_ffmpz((ulong *)Bbounds, maxBfields, ctx->minfo);
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        Abounds[i] = Abounds[i] + 1;
        Bbounds[i] = Bbounds[i] + 1;
        Pbounds[i] = Abounds[i] + Bbounds[i] - 1;
        if ((Abounds[i] | Bbounds[i] | Pbounds[i]) < 0)
        {
            goto failed_stage1;
        }
        if (i > 0)
        {
            Abounds[i] = Pbounds[i];
            Bbounds[i] = Pbounds[i];
        }
    }

    fmpz_mpolyd_init(Ad, nvars);
    fmpz_mpolyd_init(Bd, nvars);

    P_is_stolen = 0;
    if (P != A && P != B && P->alloc > 0)
    {
        /* we may steal the coeffs of P in this case to init Pd */
        Pd->nvars = nvars;
        Pd->degb_alloc = nvars;
        Pd->deg_bounds = (slong *) flint_malloc(Pd->degb_alloc*sizeof(slong));
        for (i = 0; i < nvars; i++)
        {
            Pd->deg_bounds[i] = WORD(1);
        }
        Pd->coeffs = P->coeffs;
        Pd->coeff_alloc = P->alloc;

        P->coeffs = (fmpz *) flint_calloc(P->alloc, sizeof(fmpz));

        P_is_stolen = 1;
    }
    else
    {
        fmpz_mpolyd_init(Pd, ctx->minfo->nvars);
    }

    success = 1;
    success = success && fmpz_mpolyd_set_degbounds(Ad, Abounds);
    success = success && fmpz_mpolyd_set_degbounds(Bd, Bbounds);
    success = success && fmpz_mpolyd_set_degbounds(Pd, Pbounds);
    if (!success)
    {
        goto failed_stage2;
    }

    fmpz_mpoly_convert_to_fmpz_mpolyd_degbound(Ad, A, ctx);
    fmpz_mpoly_convert_to_fmpz_mpolyd_degbound(Bd, B, ctx);

    /* let Au and Bu borrow Ad and Bd */
    Au->alloc = Ad->coeff_alloc;
    Au->coeffs = Ad->coeffs;
    Au->length = fmpz_mpolyd_length(Ad);

    Bu->alloc = Bd->coeff_alloc;
    Bu->coeffs = Bd->coeffs;
    Bu->length = fmpz_mpolyd_length(Bd);

    /* manually move P to Pu */
    Pu->alloc = Pd->coeff_alloc;
    Pu->coeffs = Pd->coeffs;
    Pu->length = 0;

    fmpz_poly_mul(Pu, Au, Bu);

    /* manually move Pu to P */
    Pd->coeff_alloc = Pu->alloc;
    Pd->coeffs = Pu->coeffs;
    Pd->length = Pu->length;

    fmpz_mpolyd_clear(Bd);
    fmpz_mpolyd_clear(Ad);
    fmpz_mpoly_consume_fmpz_mpolyd_clear(P, Pd, ctx);

done:
    TMP_END;
    return success;

failed_stage2:
    fmpz_mpolyd_clear(Ad);
    fmpz_mpolyd_clear(Bd);
    if (P_is_stolen)
    {
        fmpz * t = Pd->coeffs;
        Pd->coeffs = P->coeffs;
        P->coeffs = t;
    }
    fmpz_mpolyd_clear(Pd);

failed_stage1:
    success = 0;
    goto done;
}



int fmpz_mpoly_mul_dense(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                              const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx)
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

    if (B->bits > FLINT_BITS || C->bits > FLINT_BITS ||
        ctx->minfo->nvars < 1)
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

    success = _fmpz_mpoly_mul_dense(A, B, maxBfields, C, maxCfields, ctx);

    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_clear(maxBfields + i);
        fmpz_clear(maxCfields + i);
    }

    TMP_END;
    return success;
}

slong _fmpz_mpoly_mul_heap_part1(fmpz ** A_coeff, ulong ** A_exp, slong * A_alloc,
              const fmpz * Bcoeff, const ulong * Bexp, slong Blen,
              const fmpz * Ccoeff, const ulong * Cexp, slong Clen,
         slong * start, slong * end, slong * hind, const fmpz_mpoly_stripe_t S)
{
    const int flint_small = S->flint_small;
    const ulong cmpmask = S->cmpmask[0];
    slong i, j;
    ulong exp;
    mpoly_heap_t * x;
    slong next_loc;
    slong heap_len;
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    slong Alen;
    fmpz * Acoeff = *A_coeff;
    ulong * Aexp = *A_exp;
    slong Aalloc = *A_alloc;
    ulong acc[3], p[3];
    int first_prod;

    i = 0;
    store = store_base = (slong *) (S->big_mem + i);
    i += 2*Blen*sizeof(slong);
    heap = (mpoly_heap1_s *)(S->big_mem + i);
    i += (Blen + 1)*sizeof(mpoly_heap1_s);
    chain = (mpoly_heap_t *)(S->big_mem + i);
    i += Blen*sizeof(mpoly_heap_t);
    FLINT_ASSERT(i <= S->big_mem_alloc);

    /* put all the starting nodes on the heap */
    heap_len = 1; /* heap zero index unused */
    next_loc = Blen + 4;   /* something bigger than heap can ever be */
    for (i = 0; i < Blen; i++)
    {
        hind[i] = 2*start[i] + 1;
    }
    for (i = 0; i < Blen; i++)
    {
        if (  (start[i] < end[i])
           && (  (i == 0)
              || (start[i] < start[i - 1])
              )
           )
        {
            x = chain + i;
            x->i = i;
            x->j = start[i];
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;
            _mpoly_heap_insert1(heap, Bexp[x->i] + Cexp[x->j], x,
                                                &next_loc, &heap_len, cmpmask);
        }
    }

    Alen = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;

        _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, 1);

        Aexp[Alen] = exp;

        acc[0] = acc[1] = acc[2] = 0;
        first_prod = 1;
        while (heap_len > 1 && heap[1].exp == exp)
        {
            x = _mpoly_heap_pop1(heap, &heap_len, cmpmask);

            hind[x->i] |= WORD(1);
            *store++ = x->i;
            *store++ = x->j;

            if (flint_small)
            {
                smul_ppmm(p[1], p[0], Bcoeff[x->i], Ccoeff[x->j]);
                p[2] = FLINT_SIGN_EXT(p[1]);
                add_sssaaaaaa(acc[2], acc[1], acc[0], acc[2], acc[1], acc[0],
                                                         p[2], p[1], p[0]);
                first_prod = 0;
                while ((x = x->next) != NULL)
                {          
                    smul_ppmm(p[1], p[0], Bcoeff[x->i], Ccoeff[x->j]);
                    p[2] = FLINT_SIGN_EXT(p[1]);
                    add_sssaaaaaa(acc[2], acc[1], acc[0], acc[2], acc[1], acc[0],
                                                             p[2], p[1], p[0]);
                    hind[x->i] |= WORD(1);
                    *store++ = x->i;
                    *store++ = x->j;
                }
            }
            else /* output coeffs require multiprecision */
            {
                if (first_prod)
                    fmpz_mul(Acoeff + Alen, Bcoeff + x->i, Ccoeff + x->j);
                else
                    fmpz_addmul(Acoeff + Alen, Bcoeff + x->i, Ccoeff + x->j);
                first_prod = 0; 
                while ((x = x->next) != NULL)
                {
                    fmpz_addmul(Acoeff + Alen, Bcoeff + x->i, Ccoeff + x->j);
                    hind[x->i] |= WORD(1);
                    *store++ = x->i;
                    *store++ = x->j;
                }
            }
        }

        /* for each node temporarily stored */
        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            /* should we go right? */
            if (  (i + 1 < Blen)
               && (j + 0 < end[i + 1])
               && (hind[i + 1] == 2*j + 1)
               )
            {
                x = chain + i + 1;
                x->i = i + 1;
                x->j = j;
                x->next = NULL;
                hind[x->i] = 2*(x->j + 1) + 0;
                _mpoly_heap_insert1(heap, Bexp[x->i] + Cexp[x->j], x,
                                                &next_loc, &heap_len, cmpmask);
            }

            /* should we go up? */
            if (  (j + 1 < end[i + 0])
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
                hind[x->i] = 2*(x->j+1) + 0;
                _mpoly_heap_insert1(heap, Bexp[x->i] + Cexp[x->j], x,
                                                &next_loc, &heap_len, cmpmask);
            }
        }

        /* set output poly coeff from temporary accumulation, if not multiprec */
        if (flint_small)
        {
            fmpz_set_signed_uiuiui(Acoeff + Alen, acc[2], acc[1], acc[0]);
        }

        Alen += !fmpz_is_zero(Acoeff + Alen);
    }

    *A_coeff = Acoeff;
    *A_exp = Aexp;
    *A_alloc = Aalloc;
    return Alen;
}


slong _fmpz_mpoly_mul_heap_part(fmpz ** A_coeff, ulong ** A_exp, slong * A_alloc,
                 const fmpz * Bcoeff, const ulong * Bexp, slong Blen,
                 const fmpz * Ccoeff, const ulong * Cexp, slong Clen,
         slong * start, slong * end, slong * hind, const fmpz_mpoly_stripe_t S)
{
    const int flint_small = S->flint_small;
    flint_bitcnt_t bits = S->bits;
    slong N = S->N;
    const ulong * cmpmask = S->cmpmask;
    slong i, j;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    mpoly_heap_t * x;
    slong next_loc;
    slong heap_len;
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    slong Alen;
    ulong * Aexp = *A_exp;
    slong Aalloc = *A_alloc;
    fmpz * Acoeff = *A_coeff;
    ulong acc[3], p[3];
    int first_prod;

    i = 0;
    store = store_base = (slong *) (S->big_mem + i);
    i += 2*Blen*sizeof(slong);
    exp_list = (ulong **) (S->big_mem + i);
    i += Blen*sizeof(ulong *);
    exps = (ulong *) (S->big_mem + i);
    i += Blen*N*sizeof(ulong);
    heap = (mpoly_heap_s *) (S->big_mem + i);
    i += (Blen + 1)*sizeof(mpoly_heap_s);
    chain = (mpoly_heap_t *) (S->big_mem + i);
    i += Blen*sizeof(mpoly_heap_t);
    FLINT_ASSERT(i <= S->big_mem_alloc);

    /* put all the starting nodes on the heap */
    heap_len = 1; /* heap zero index unused */
    next_loc = Blen + 4;   /* something bigger than heap can ever be */
    exp_next = 0;
    for (i = 0; i < Blen; i++)
        exp_list[i] = exps + i*N;
    for (i = 0; i < Blen; i++)
        hind[i] = 2*start[i] + 1;
    for (i = 0; i < Blen; i++)
    {
        if (  (start[i] < end[i])
           && (  (i == 0)
              || (start[i] < start[i - 1])
              )
           )
        {
            x = chain + i;
            x->i = i;
            x->j = start[i];
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;
            if (bits <= FLINT_BITS)
                mpoly_monomial_add(exp_list[exp_next], Bexp + x->i*N,
                                                       Cexp + x->j*N, N);
            else
                mpoly_monomial_add_mp(exp_list[exp_next], Bexp + x->i*N,
                                                          Cexp + x->j*N, N);

            exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
        }
    }

    Alen = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;

        _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, N);

        mpoly_monomial_set(Aexp + N*Alen, exp, N);

        acc[0] = acc[1] = acc[2] = 0;
        first_prod = 1;
        while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N))
        {
            exp_list[--exp_next] = heap[1].exp;

            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);

            hind[x->i] |= WORD(1);
            *store++ = x->i;
            *store++ = x->j;

            /* if output coeffs will fit in three words */
            if (flint_small)
            {
                smul_ppmm(p[1], p[0], Bcoeff[x->i], Ccoeff[x->j]);
                p[2] = FLINT_SIGN_EXT(p[1]);
                add_sssaaaaaa(acc[2], acc[1], acc[0], acc[2], acc[1], acc[0],
                                                         p[2], p[1], p[0]);
                first_prod = 0;
                while ((x = x->next) != NULL)
                {
                    smul_ppmm(p[1], p[0], Bcoeff[x->i], Ccoeff[x->j]);
                    p[2] = FLINT_SIGN_EXT(p[1]);
                    add_sssaaaaaa(acc[2], acc[1], acc[0], acc[2], acc[1], acc[0],
                                                             p[2], p[1], p[0]);
                    hind[x->i] |= WORD(1);
                    *store++ = x->i;
                    *store++ = x->j;
                }
            }
            else /* output coeffs require multiprecision */
            {
                if (first_prod)
                    fmpz_mul(Acoeff + Alen, Bcoeff + x->i, Ccoeff + x->j);
                else
                    fmpz_addmul(Acoeff + Alen, Bcoeff + x->i, Ccoeff + x->j);
                first_prod = 0;
                while ((x = x->next) != NULL)
                {
                    fmpz_addmul(Acoeff + Alen, Bcoeff + x->i, Ccoeff + x->j);
                    hind[x->i] |= WORD(1);
                    *store++ = x->i;
                    *store++ = x->j;
                }
            }
        }
      
        /* for each node temporarily stored */
        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            /* should we go right? */
            if (  (i + 1 < Blen)
               && (j + 0 < end[i + 1])
               && (hind[i + 1] == 2*j + 1)
               )
            {
                x = chain + i + 1;
                x->i = i + 1;
                x->j = j;
                x->next = NULL;

                hind[x->i] = 2*(x->j + 1) + 0;

                if (bits <= FLINT_BITS)
                    mpoly_monomial_add(exp_list[exp_next], Bexp + x->i*N,
                                                           Cexp + x->j*N, N);
                else
                    mpoly_monomial_add_mp(exp_list[exp_next], Bexp + x->i*N,
                                                              Cexp + x->j*N, N);

                exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
            }

            /* should we go up? */
            if (  (j + 1 < end[i + 0])
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

                if (bits <= FLINT_BITS)
                    mpoly_monomial_add(exp_list[exp_next], Bexp + x->i*N,
                                                           Cexp + x->j*N, N);
                else
                    mpoly_monomial_add_mp(exp_list[exp_next], Bexp + x->i*N,
                                                              Cexp + x->j*N, N);

                exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
            }
        }     

        /* set output poly coeff from temporary accumulation, if not multiprec */
        if (flint_small)
        {
            fmpz_set_signed_uiuiui(Acoeff + Alen, acc[2], acc[1], acc[0]);
        }

        Alen += !fmpz_is_zero(Acoeff + Alen);
    }

    *A_coeff = Acoeff;
    *A_exp = Aexp;
    *A_alloc = Aalloc;
    return Alen;
}


/*
    The workers calculate product terms from 4*n divisions, where n is the
    number of threads.
*/

typedef struct
{
    volatile int idx;
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
#endif
    slong nthreads;
    slong ndivs;
    fmpz * Acoeff;
    ulong * Aexp;
    const fmpz * Bcoeff;
    const ulong * Bexp;
    slong Blen;
    const fmpz * Ccoeff;
    const ulong * Cexp;
    slong Clen;
    slong N;
    flint_bitcnt_t bits;
    const ulong * cmpmask;
    int flint_small;
}
_base_struct2;

typedef _base_struct2 _base_t2[1];

typedef struct
{
    slong lower;
    slong upper;
    slong thread_idx;
    slong Aoffset;
    slong Alen;
    slong Aalloc;
    ulong * Aexp;
    fmpz * Acoeff;
}
_div_struct;

typedef struct
{
    fmpz_mpoly_stripe_t S;
    slong idx;
    slong time;
    _base_struct2 * base;
    _div_struct * divs;
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
    pthread_cond_t cond;
#endif
    slong * t1, * t2, * t3, * t4;
    ulong * exp;
}
_worker_arg_struct2;


/*
    The workers simply take the next available division and calculate all
    product terms in this division.
*/

#define SWAP_PTRS(xx, yy) \
   do { \
      tt = xx; \
      xx = yy; \
      yy = tt; \
   } while (0)

static void _fmpz_mpoly_mul_heap_threaded_worker(void * varg)
{
    _worker_arg_struct2 * arg = (_worker_arg_struct2 *) varg;
    fmpz_mpoly_stripe_struct * S = arg->S;
    _div_struct * divs = arg->divs;
    _base_struct2 * base = arg->base;
    slong Blen = base->Blen;
    slong N = base->N;
    slong i, j;
    ulong *exp;
    slong score;
    slong *start, *end, *t1, *t2, *t3, *t4, *tt;

    exp = (ulong *) flint_malloc(N*sizeof(ulong));
    t1 = (slong *) flint_malloc(Blen*sizeof(slong));
    t2 = (slong *) flint_malloc(Blen*sizeof(slong));
    t3 = (slong *) flint_malloc(Blen*sizeof(slong));
    t4 = (slong *) flint_malloc(Blen*sizeof(slong));

    S->N = N;
    S->bits = base->bits;
    S->cmpmask = base->cmpmask;
    S->flint_small = base->flint_small;

    S->big_mem_alloc = 0;
    if (N == 1)
    {
        S->big_mem_alloc += 2*Blen*sizeof(slong);
        S->big_mem_alloc += (Blen + 1)*sizeof(mpoly_heap1_s);
        S->big_mem_alloc += Blen*sizeof(mpoly_heap_t);
    }
    else
    {
        S->big_mem_alloc += 2*Blen*sizeof(slong);
        S->big_mem_alloc += (Blen + 1)*sizeof(mpoly_heap_s);
        S->big_mem_alloc += Blen*sizeof(mpoly_heap_t);
        S->big_mem_alloc += Blen*S->N*sizeof(ulong);
        S->big_mem_alloc += Blen*sizeof(ulong *);
    }
    S->big_mem = (char *) flint_malloc(S->big_mem_alloc);

    /* get index to start working on */
    if (arg->idx + 1 < base->nthreads)
    {
#if FLINT_USES_PTHREAD
        pthread_mutex_lock(&base->mutex);
#endif
        i = base->idx - 1;
        base->idx = i;
#if FLINT_USES_PTHREAD
        pthread_mutex_unlock(&base->mutex);
#endif
    }
    else
    {
        i = base->ndivs - 1;
    }

    while (i >= 0)
    {
        FLINT_ASSERT(divs[i].thread_idx == -WORD(1));
        divs[i].thread_idx = arg->idx;

        /* calculate start */
        if (i + 1 < base-> ndivs)
        {
            mpoly_search_monomials(
                &start, exp, &score, t1, t2, t3,
                            divs[i].lower, divs[i].lower,
                            base->Bexp, base->Blen, base->Cexp, base->Clen,
                                          base->N, base->cmpmask);
            if (start == t2)
            {
                SWAP_PTRS(t1, t2);
            }
            else if (start == t3)
            {
                SWAP_PTRS(t1, t3);
            }
        }
        else
        {
            start = t1;
            for (j = 0; j < base->Blen; j++)
                start[j] = 0;
        }

        /* calculate end */
        if (i > 0)
        {
            mpoly_search_monomials(
                &end, exp, &score, t2, t3, t4,
                            divs[i - 1].lower, divs[i - 1].lower,
                            base->Bexp, base->Blen, base->Cexp, base->Clen,
                                          base->N, base->cmpmask);
            if (end == t3)
            {
                SWAP_PTRS(t2, t3);
            }
            else if (end == t4)
            {
                SWAP_PTRS(t2, t4);
            }
        }
        else
        {
            end = t2;
            for (j = 0; j < base->Blen; j++)
                end[j] = base->Clen;
        }
        /* t3 and t4 are free for workspace at this point */

        /* join code assumes all divisions have been allocated */
        _fmpz_mpoly_fit_length(&divs[i].Acoeff, &divs[i].Aexp, &divs[i].Aalloc,
                                                                       256, N);
        /* calculate products in [start, end) */
        if (N == 1)
        {
            divs[i].Alen = _fmpz_mpoly_mul_heap_part1(
                         &divs[i].Acoeff, &divs[i].Aexp, &divs[i].Aalloc,
                              base->Bcoeff, base->Bexp, base->Blen,
                              base->Ccoeff, base->Cexp, base->Clen,
                                                          start, end, t3, S);
        }
        else
        {
            divs[i].Alen = _fmpz_mpoly_mul_heap_part(
                         &divs[i].Acoeff, &divs[i].Aexp, &divs[i].Aalloc,
                              base->Bcoeff, base->Bexp, base->Blen,
                              base->Ccoeff, base->Cexp, base->Clen,
                                                          start, end, t3, S);
        }

        /* get next index to work on */
#if FLINT_USES_PTHREAD
        pthread_mutex_lock(&base->mutex);
#endif
	i = base->idx - 1;
        base->idx = i;
#if FLINT_USES_PTHREAD
	pthread_mutex_unlock(&base->mutex);
#endif
    }

    flint_free(S->big_mem);
    flint_free(t4);
    flint_free(t3);
    flint_free(t2);
    flint_free(t1);
    flint_free(exp);
}

static void _join_worker(void * varg)
{
    _worker_arg_struct2 * arg = (_worker_arg_struct2 *) varg;
    _div_struct * divs = arg->divs;
    _base_struct2 * base = arg->base;
    slong N = base->N;
    slong i;

    for (i = base->ndivs - 2; i >= 0; i--)
    {
        FLINT_ASSERT(divs[i].thread_idx != -WORD(1));

        if (divs[i].thread_idx != arg->idx)
            continue;

        FLINT_ASSERT(divs[i].Acoeff != NULL);
        FLINT_ASSERT(divs[i].Aexp != NULL);

        memcpy(base->Acoeff + divs[i].Aoffset, divs[i].Acoeff,
                                                    divs[i].Alen*sizeof(fmpz));

        memcpy(base->Aexp + N*divs[i].Aoffset, divs[i].Aexp,
                                                 N*divs[i].Alen*sizeof(ulong));

        flint_free(divs[i].Acoeff);
        flint_free(divs[i].Aexp);
    }
}

void _fmpz_mpoly_mul_heap_threaded(
    fmpz_mpoly_t A,
    const fmpz * Bcoeff, const ulong * Bexp, slong Blen,
    const fmpz * Ccoeff, const ulong * Cexp, slong Clen,
    flint_bitcnt_t bits,
    slong N,
    const ulong * cmpmask,
    const thread_pool_handle * handles,
    slong num_handles)
{
    slong i, j;
    slong BClen, hi;
    _base_t2 base;
    _div_struct * divs;
    _worker_arg_struct2 * args;
    slong Aalloc;
    slong Alen;
    fmpz * Acoeff;
    ulong * Aexp;

    /* bail if product of lengths overflows a word */
    umul_ppmm(hi, BClen, Blen, Clen);
    if (hi != 0 || BClen < 0)
    {
        Alen = _fmpz_mpoly_mul_johnson(&A->coeffs, &A->exps, &A->alloc,
                                         Bcoeff, Bexp, Blen,
                                         Ccoeff, Cexp, Clen, bits, N, cmpmask);
        _fmpz_mpoly_set_length(A, Alen, NULL);
        return;

    }

    base->nthreads = num_handles + 1;
    base->ndivs = base->nthreads*4;  /* number of divisions */
    base->Bcoeff = Bcoeff;
    base->Bexp = Bexp;
    base->Blen = Blen;
    base->Ccoeff = Ccoeff;
    base->Cexp = Cexp;
    base->Clen = Clen;
    base->bits = bits;
    base->N = N;
    base->cmpmask = cmpmask;
    base->idx = base->ndivs - 1;    /* decremented by worker threads */
    base->flint_small =   _fmpz_mpoly_fits_small(Bcoeff, Blen)
                       && _fmpz_mpoly_fits_small(Ccoeff, Clen);

    divs = (_div_struct *) flint_malloc(base->ndivs*sizeof(_div_struct));
    args = (_worker_arg_struct2 *) flint_malloc(base->nthreads
                                                  *sizeof(_worker_arg_struct2));

    /* allocate space and set the boundary for each division */
    FLINT_ASSERT(BClen/Blen == Clen);
    for (i = base->ndivs - 1; i >= 0; i--)
    {
        double d = (double)(i + 1) / (double)(base->ndivs);

        /* divisions decrease in size so that no worker finishes too early */
        divs[i].lower = (d * d) * BClen;
        divs[i].lower = FLINT_MIN(divs[i].lower, BClen);
        divs[i].lower = FLINT_MAX(divs[i].lower, WORD(0));
        divs[i].upper = divs[i].lower;
        divs[i].Aoffset = -WORD(1);
        divs[i].thread_idx = -WORD(1);

        divs[i].Alen = 0;
        if (i == base->ndivs - 1)
        {
            /* highest division writes to original poly */
            divs[i].Aalloc = A->alloc;
            divs[i].Aexp = A->exps;
            divs[i].Acoeff = A->coeffs;
            /* must clear output coefficients before working in parallel */
            for (j = 0; j < A->length; j++)
               _fmpz_demote(A->coeffs + j);
        }
        else
        {
            /* lower divisions write to a new worker poly */
            divs[i].Aalloc = 0;
            divs[i].Aexp = NULL;
            divs[i].Acoeff = NULL;
        }
    }

    /* compute each chunk in parallel */
#if FLINT_USES_PTHREAD
    pthread_mutex_init(&base->mutex, NULL);
#endif
    for (i = 0; i < num_handles; i++)
    {
        args[i].idx = i;
        args[i].base = base;
        args[i].divs = divs;
        thread_pool_wake(global_thread_pool, handles[i], 0,
                               _fmpz_mpoly_mul_heap_threaded_worker, &args[i]);
    }
    i = num_handles;
    args[i].idx = i;
    args[i].base = base;
    args[i].divs = divs;
    _fmpz_mpoly_mul_heap_threaded_worker(&args[i]);
    for (i = 0; i < num_handles; i++)
    {
        thread_pool_wait(global_thread_pool, handles[i]);
    }

    /* calculate and allocate space for final answer */
    i = base->ndivs - 1;
    Alen = divs[i].Alen;
    Acoeff = divs[i].Acoeff;
    Aexp = divs[i].Aexp;
    Aalloc = divs[i].Aalloc;
    for (i = base->ndivs - 2; i >= 0; i--)
    {
        divs[i].Aoffset = Alen;
        Alen += divs[i].Alen;
    }
    if (Alen > Aalloc)
    {
        Acoeff = (fmpz *) flint_realloc(Acoeff, Alen*sizeof(fmpz));
        Aexp = (ulong *) flint_realloc(Aexp, Alen*N*sizeof(ulong));
        Aalloc = Alen;
    }
    base->Acoeff = Acoeff;
    base->Aexp = Aexp;

    /* join answers */
    for (i = 0; i < num_handles; i++)
    {
        thread_pool_wake(global_thread_pool, handles[i], 0, _join_worker, &args[i]);
    }
    _join_worker(&args[num_handles]);
    for (i = 0; i < num_handles; i++)
    {
        thread_pool_wait(global_thread_pool, handles[i]);
    }

#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&base->mutex);
#endif

    flint_free(args);
    flint_free(divs);

    /* we should have managed to keep coefficients past length demoted */
    FLINT_ASSERT(Alen <= Aalloc);
#if FLINT_WANT_ASSERT
    for (i = Alen; i < Aalloc; i++)
    {
        FLINT_ASSERT(!COEFF_IS_MPZ(*(Acoeff + i)));
    }
#endif

    A->coeffs = Acoeff;
    A->exps = Aexp;
    A->alloc = Aalloc;
    A->length = Alen;
}



/* maxBfields gets clobbered */
void _fmpz_mpoly_mul_heap_threaded_pool_maxfields(
    fmpz_mpoly_t A,
    const fmpz_mpoly_t B, fmpz * maxBfields,
    const fmpz_mpoly_t C, fmpz * maxCfields,
    const fmpz_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles)
{
    slong N;
    flint_bitcnt_t exp_bits;
    ulong * cmpmask;
    ulong * Bexp, * Cexp;
    int freeBexp, freeCexp;
    TMP_INIT;

    TMP_START;

    _fmpz_vec_add(maxBfields, maxBfields, maxCfields, ctx->minfo->nfields);

    exp_bits = _fmpz_vec_max_bits(maxBfields, ctx->minfo->nfields);
    exp_bits = FLINT_MAX(MPOLY_MIN_BITS, exp_bits + 1);
    exp_bits = FLINT_MAX(exp_bits, B->bits);
    exp_bits = FLINT_MAX(exp_bits, C->bits);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);

    N = mpoly_words_per_exp(exp_bits, ctx->minfo);
    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, exp_bits, ctx->minfo);

    /* ensure input exponents are packed into same sized fields as output */
    freeBexp = 0;
    Bexp = B->exps;
    if (exp_bits > B->bits)
    {
        freeBexp = 1;
        Bexp = (ulong *) flint_malloc(N*B->length*sizeof(ulong));
        mpoly_repack_monomials(Bexp, exp_bits, B->exps, B->bits,
                                                        B->length, ctx->minfo);
    }

    freeCexp = 0;
    Cexp = C->exps;
    if (exp_bits > C->bits)
    {
        freeCexp = 1;
        Cexp = (ulong *) flint_malloc(N*C->length*sizeof(ulong));
        mpoly_repack_monomials(Cexp, exp_bits, C->exps, C->bits,
                                                        C->length, ctx->minfo);
    }

    /* deal with aliasing and do multiplication */
    if (A == B || A == C)
    {
        fmpz_mpoly_t T;
        fmpz_mpoly_init3(T, 0, exp_bits, ctx);

        /* algorithm more efficient if smaller poly first */
        if (B->length >= C->length)
        {
            _fmpz_mpoly_mul_heap_threaded(T, C->coeffs, Cexp, C->length,
                                             B->coeffs, Bexp, B->length,
                                   exp_bits, N, cmpmask, handles, num_handles);
        }
        else
        {
            _fmpz_mpoly_mul_heap_threaded(T, B->coeffs, Bexp, B->length,
                                             C->coeffs, Cexp, C->length,
                                   exp_bits, N, cmpmask, handles, num_handles);
        }

        fmpz_mpoly_swap(T, A, ctx);
        fmpz_mpoly_clear(T, ctx);
    }
    else
    {
        fmpz_mpoly_fit_length_reset_bits(A, A->length, exp_bits, ctx);

        /* algorithm more efficient if smaller poly first */
        if (B->length > C->length)
        {
            _fmpz_mpoly_mul_heap_threaded(A, C->coeffs, Cexp, C->length,
                                             B->coeffs, Bexp, B->length,
                                   exp_bits, N, cmpmask, handles, num_handles);
        }
        else
        {
            _fmpz_mpoly_mul_heap_threaded(A, B->coeffs, Bexp, B->length,
                                             C->coeffs, Cexp, C->length,
                                   exp_bits, N, cmpmask, handles, num_handles);
        }
    }

    if (freeBexp)
        flint_free(Bexp);

    if (freeCexp)
        flint_free(Cexp);

    TMP_END;
}


void fmpz_mpoly_mul_heap_threaded(
    fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_t C,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    fmpz * maxBfields, * maxCfields;
    thread_pool_handle * handles;
    slong num_handles;
    slong thread_limit = FLINT_MIN(A->length, B->length)/16;
    TMP_INIT;

    if (B->length == 0 || C->length == 0)
    {
        fmpz_mpoly_zero(A, ctx);
        return;
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

    num_handles = flint_request_threads(&handles, thread_limit);

    _fmpz_mpoly_mul_heap_threaded_pool_maxfields(A, B, maxBfields, C, maxCfields,
                                                    ctx, handles, num_handles);

    flint_give_back_threads(handles, num_handles);

    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_clear(maxBfields + i);
        fmpz_clear(maxCfields + i);
    }

    TMP_END;
}

/*
   Set poly1 to poly2*poly3 using Johnson's heap method. The function
   realocates its output and returns the length of the product. This
   version of the function assumes the exponent vectors all fit in a
   single word. Assumes input polys are nonzero.
*/
slong _fmpz_mpoly_mul_johnson1(fmpz ** poly1, ulong ** exp1, slong * alloc,
              const fmpz * poly2, const ulong * exp2, slong len2,
              const fmpz * poly3, const ulong * exp3, slong len3, ulong maskhi)
{
   slong i, j, k;
   slong next_loc;
   slong Q_len = 0, heap_len = 2; /* heap zero index unused */
   mpoly_heap1_s * heap;
   mpoly_heap_t * chain;
   slong * Q;
   mpoly_heap_t * x;
   fmpz * p1 = *poly1;
   ulong * e1 = *exp1;
   slong * hind;
   ulong exp, cy;
   ulong c[3], p[2]; /* for accumulating coefficients */
   int first, small;
   TMP_INIT;

   TMP_START;

   /* whether input coeffs are small, thus output coeffs fit in three words */
   small = _fmpz_mpoly_fits_small(poly2, len2) &&
                                           _fmpz_mpoly_fits_small(poly3, len3);

   next_loc = len2 + 4;   /* something bigger than heap can ever be */
   heap = (mpoly_heap1_s *) TMP_ALLOC((len2 + 1)*sizeof(mpoly_heap1_s));
   /* alloc array of heap nodes which can be chained together */
   chain = (mpoly_heap_t *) TMP_ALLOC(len2*sizeof(mpoly_heap_t));
   /* space for temporary storage of pointers to heap nodes */
   Q = (slong *) TMP_ALLOC(2*len2*sizeof(slong));

    /* space for heap indices */
    hind = (slong *) TMP_ALLOC(len2*sizeof(slong));
    for (i = 0; i < len2; i++)
        hind[i] = 1;

   /* put (0, 0, exp2[0] + exp3[0]) on heap */
   x = chain + 0;
   x->i = 0;
   x->j = 0;
   x->next = NULL;

   HEAP_ASSIGN(heap[1], exp2[0] + exp3[0], x);
   hind[0] = 2*1 + 0;

   /* output poly index starts at -1, will be immediately updated to 0 */
   k = -WORD(1);

   /* while heap is nonempty */
   while (heap_len > 1)
   {
      /* get exponent field of heap top */
      exp = heap[1].exp;
      
      /* realloc output poly ready for next product term */
      k++;
      _fmpz_mpoly_fit_length(&p1, &e1, alloc, k + 1, 1);

      /* whether we are on first coeff product for this output exponent */
      first = 1;

      /* set temporary coeff to zero */
      c[0] = c[1] = c[2] = 0;

      /* while heap nonempty and contains chain with current output exponent */
      while (heap_len > 1 && heap[1].exp == exp)
      {
         /* pop chain from heap */
         x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
         
         /* take node out of heap and put into store */
         hind[x->i] |= WORD(1);
         Q[Q_len++] = x->i;
         Q[Q_len++] = x->j;

         /* if output coeffs will fit in three words */
         if (small)
         {
            /* compute product of input poly coeffs */
            if (first)
            {
               smul_ppmm(c[1], c[0], poly2[x->i], poly3[x->j]);
               c[2] = -(c[1] >> (FLINT_BITS - 1));

               /* set output monomial */
               e1[k] = exp;
               first = 0; 
            } else /* addmul product of input poly coeffs */
            {
               smul_ppmm(p[1], p[0], poly2[x->i], poly3[x->j]);
               add_sssaaaaaa(cy, c[1], c[0], 0, c[1], c[0], 0, p[1], p[0]);
               c[2] += (0 <= (slong) p[1]) ? cy : cy - 1;
            }

            /* for every node in this chain */
            while ((x = x->next) != NULL)
            {
               /* addmul product of input poly coeffs */
               smul_ppmm(p[1], p[0], poly2[x->i], poly3[x->j]);
               add_sssaaaaaa(cy, c[1], c[0], 0, c[1], c[0], 0, p[1], p[0]);
               c[2] += (0 <= (slong) p[1]) ? cy : cy - 1;

               /* take node out of heap and put into store */
               hind[x->i] |= WORD(1);
               Q[Q_len++] = x->i;
               Q[Q_len++] = x->j;
            }
         } else /* output coeffs require multiprecision */
         {
            if (first) /* compute product of input poly coeffs */
            {
               fmpz_mul(p1 + k, poly2 + x->i, poly3 + x->j);
               
               e1[k] = exp;
               first = 0; 
            } else
            {  /* addmul product of input poly coeffs */
               fmpz_addmul(p1 + k, poly2 + x->i, poly3 + x->j);
            }

            /* for each node in this chain */
            while ((x = x->next) != NULL)
            {
               /* addmul product of input poly coeffs */
               fmpz_addmul(p1 + k, poly2 + x->i, poly3 + x->j);

               /* take node out of heap and put into store */
               hind[x->i] |= WORD(1);
               Q[Q_len++] = x->i;
               Q[Q_len++] = x->j;
            }
         }
      }
      
      /* for each node temporarily stored */
      while (Q_len > 0)
      {
         /* take node from store */
         j = Q[--Q_len];
         i = Q[--Q_len];

         /* should we go right? */
         if (  (i + 1 < len2)
            && (hind[i + 1] == 2*j + 1)
            )
         {
            x = chain + i + 1;
            x->i = i + 1;
            x->j = j;
            x->next = NULL;

            hind[x->i] = 2*(x->j+1) + 0;
            _mpoly_heap_insert1(heap, exp2[x->i] + exp3[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
         }

         /* should we go up? */
         if (  (j + 1 < len3)
            && ((hind[i] & 1) == 1)
            && (  (i == 0)
               || (hind[i - 1] >  2*(j + 2) + 1)
               || (hind[i - 1] == 2*(j + 2) + 1) /* gcc should fuse */
               )
            )
         {
            x = chain + i;
            x->i = i;
            x->j = j + 1;
            x->next = NULL;

            hind[x->i] = 2*(x->j+1) + 0;
            _mpoly_heap_insert1(heap, exp2[x->i] + exp3[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
         }
      }

      /* set output poly coeff from temporary accumulation, if not multiprec */
      if (small)
         fmpz_set_signed_uiuiui(p1 + k, c[2], c[1], c[0]);

      if (fmpz_is_zero(p1 + k))
         k--;
   }

   k++;

   (*poly1) = p1;
   (*exp1) = e1;
   
   TMP_END;

   return k;
}

/*
   Set poly1 to poly2*poly3 using Johnson's heap method. The function
   realocates its output and returns the length of the product. This
   version of the function assumes the exponent vectors take N words.
*/
slong _fmpz_mpoly_mul_johnson(fmpz ** poly1, ulong ** exp1, slong * alloc,
                 const fmpz * poly2, const ulong * exp2, slong len2,
                 const fmpz * poly3, const ulong * exp3, slong len3,
                              flint_bitcnt_t bits, slong N, const ulong * cmpmask)
{
   slong i, j, k;
   slong next_loc;
   slong Q_len = 0, heap_len = 2; /* heap zero index unused */
   mpoly_heap_s * heap;
   mpoly_heap_t * chain;
   slong * Q;
   mpoly_heap_t * x;
   fmpz * p1 = *poly1;
   ulong * e1 = *exp1;
   ulong cy;
   ulong c[3], p[2]; /* for accumulating coefficients */
   ulong * exp, * exps;
   ulong ** exp_list;
   slong exp_next;
   slong * hind;
   int first, small;
   TMP_INIT;

   /* if exponent vectors fit in single word, call special version */
   if (N == 1)
      return _fmpz_mpoly_mul_johnson1(poly1, exp1, alloc,
                             poly2, exp2, len2, poly3, exp3, len3, cmpmask[0]);

   TMP_START;

   /* whether input coeffs are small, thus output coeffs fit in three words */
   small = _fmpz_mpoly_fits_small(poly2, len2) &&
                                           _fmpz_mpoly_fits_small(poly3, len3);

   next_loc = len2 + 4;   /* something bigger than heap can ever be */
   heap = (mpoly_heap_s *) TMP_ALLOC((len2 + 1)*sizeof(mpoly_heap_s));
   /* alloc array of heap nodes which can be chained together */
   chain = (mpoly_heap_t *) TMP_ALLOC(len2*sizeof(mpoly_heap_t));
   /* space for temporary storage of pointers to heap nodes */
   Q = (slong *) TMP_ALLOC(2*len2*sizeof(slong));
   /* allocate space for exponent vectors of N words */
   exps = (ulong *) TMP_ALLOC(len2*N*sizeof(ulong));
   /* list of pointers to allocated exponent vectors */
   exp_list = (ulong **) TMP_ALLOC(len2*sizeof(ulong *));
   for (i = 0; i < len2; i++)
      exp_list[i] = exps + i*N;

   /* space for heap indices */
   hind = (slong *) TMP_ALLOC(len2*sizeof(slong));
   for (i = 0; i < len2; i++)
       hind[i] = 1;

   /* start with no heap nodes and no exponent vectors in use */
   exp_next = 0;

   /* put (0, 0, exp2[0] + exp3[0]) on heap */
   x = chain + 0;
   x->i = 0;
   x->j = 0;
   x->next = NULL;

   heap[1].next = x;
   heap[1].exp = exp_list[exp_next++];

    if (bits <= FLINT_BITS)
        mpoly_monomial_add(heap[1].exp, exp2, exp3, N);
    else
        mpoly_monomial_add_mp(heap[1].exp, exp2, exp3, N);

    hind[0] = 2*1 + 0;

   /* output poly index starts at -1, will be immediately updated to 0 */
   k = -WORD(1);

   /* while heap is nonempty */
   while (heap_len > 1)
   {
      /* get pointer to exponent field of heap top */
      exp = heap[1].exp;

      /* realloc output poly ready for next product term */
      k++;
      _fmpz_mpoly_fit_length(&p1, &e1, alloc, k + 1, N);

      /* whether we are on first coeff product for this output exponent */
      first = 1;

      /* set temporary coeff to zero */
      c[0] = c[1] = c[2] = 0;

      /* while heap nonempty and contains chain with current output exponent */
      while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N))
      {
         /* pop chain from heap and set exponent field to be reused */
         exp_list[--exp_next] = heap[1].exp;

         x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);

         /* take node out of heap and put into store */
         hind[x->i] |= WORD(1);
         Q[Q_len++] = x->i;
         Q[Q_len++] = x->j;

         /* if output coeffs will fit in three words */
         if (small)
         {
            /* compute product of input poly coeffs */
            if (first)
            {
               smul_ppmm(c[1], c[0], poly2[x->i], poly3[x->j]);
               c[2] = -(c[1] >> (FLINT_BITS - 1));

               /* set output monomial */
               mpoly_monomial_set(e1 + k*N, exp, N);

               first = 0; 
            } else /* addmul product of input poly coeffs */
            {
               smul_ppmm(p[1], p[0], poly2[x->i], poly3[x->j]);
               add_sssaaaaaa(cy, c[1], c[0], 0, c[1], c[0], 0, p[1], p[0]);
               c[2] += (0 <= (slong) p[1]) ? cy : cy - 1;
            }
      
            /* for every node in this chain */
            while ((x = x->next) != NULL)
            {
               /* addmul product of input poly coeffs */
               smul_ppmm(p[1], p[0], poly2[x->i], poly3[x->j]);
               add_sssaaaaaa(cy, c[1], c[0], 0, c[1], c[0], 0, p[1], p[0]);
               c[2] += (0 <= (slong) p[1]) ? cy : cy - 1;

               /* take node out of heap and put into store */
               hind[x->i] |= WORD(1);
               Q[Q_len++] = x->i;
               Q[Q_len++] = x->j;
            }
         } else /* output coeffs require multiprecision */
         {
            if (first) /* compute product of input poly coeffs */
            {
               fmpz_mul(p1 + k, poly2 + x->i, poly3 + x->j);
               
               /* set output monomial */
               mpoly_monomial_set(e1 + k*N, exp, N);

               first = 0; 
            } else
            {  /* addmul product of input poly coeffs */
               fmpz_addmul(p1 + k, poly2 + x->i, poly3 + x->j);
            }

            /* for each node in this chain */
            while ((x = x->next) != NULL)
            {
               /* addmul product of input poly coeffs */
               fmpz_addmul(p1 + k, poly2 + x->i, poly3 + x->j);

               /* take node out of heap and put into store */
               hind[x->i] |= WORD(1);
               Q[Q_len++] = x->i;
               Q[Q_len++] = x->j;
            }
         }
      }

      /* for each node temporarily stored */
      while (Q_len > 0)
      {
         /* take node from store */
         j = Q[--Q_len];
         i = Q[--Q_len];

         /* should we go right? */
         if (  (i + 1 < len2)
            && (hind[i + 1] == 2*j + 1)
            )
         {
            x = chain + i + 1;
            x->i = i + 1;
            x->j = j;
            x->next = NULL;

            hind[x->i] = 2*(x->j+1) + 0;

            if (bits <= FLINT_BITS)
                mpoly_monomial_add(exp_list[exp_next], exp2 + x->i*N, exp3 + x->j*N, N);
            else
                mpoly_monomial_add_mp(exp_list[exp_next], exp2 + x->i*N, exp3 + x->j*N, N);

            if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                      &next_loc, &heap_len, N, cmpmask))
               exp_next--;
         }

         /* should we go up? */
         if (  (j + 1 < len3)
            && ((hind[i] & 1) == 1)
            && (  (i == 0)
               || (hind[i - 1] >  2*(j + 2) + 1)
               || (hind[i - 1] == 2*(j + 2) + 1) /* gcc should fuse */
               )
            )
         {
            x = chain + i;
            x->i = i;
            x->j = j + 1;
            x->next = NULL;

            hind[x->i] = 2*(x->j+1) + 0;

            if (bits <= FLINT_BITS)
                mpoly_monomial_add(exp_list[exp_next], exp2 + x->i*N, exp3 + x->j*N, N);
            else
                mpoly_monomial_add_mp(exp_list[exp_next], exp2 + x->i*N, exp3 + x->j*N, N);

            if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                      &next_loc, &heap_len, N, cmpmask))
               exp_next--;
         }
      }

      /* set output poly coeff from temporary accumulation, if not multiprec */
      if (small)
         fmpz_set_signed_uiuiui(p1 + k, c[2], c[1], c[0]);

      if (fmpz_is_zero(p1 + k))
         k--;
   }

   k++;

   (*poly1) = p1;
   (*exp1) = e1;
   
   TMP_END;

   return k;
}

/* maxBfields gets clobbered */
void _fmpz_mpoly_mul_johnson_maxfields(
    fmpz_mpoly_t A,
    const fmpz_mpoly_t B, fmpz * maxBfields,
    const fmpz_mpoly_t C, fmpz * maxCfields,
    const fmpz_mpoly_ctx_t ctx)
{
    slong N, Alen;
    flint_bitcnt_t Abits;
    ulong * cmpmask;
    ulong * Bexp, * Cexp;
    int freeBexp, freeCexp;
    TMP_INIT;

    TMP_START;

    _fmpz_vec_add(maxBfields, maxBfields, maxCfields, ctx->minfo->nfields);

    Abits = _fmpz_vec_max_bits(maxBfields, ctx->minfo->nfields);
    Abits = FLINT_MAX(MPOLY_MIN_BITS, Abits + 1);
    Abits = FLINT_MAX(Abits, B->bits);
    Abits = FLINT_MAX(Abits, C->bits);
    Abits = mpoly_fix_bits(Abits, ctx->minfo);

    N = mpoly_words_per_exp(Abits, ctx->minfo);
    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, Abits, ctx->minfo);

    /* ensure input exponents are packed into same sized fields as output */
    freeBexp = 0;
    Bexp = B->exps;
    if (Abits > B->bits)
    {
        freeBexp = 1;
        Bexp = (ulong *) flint_malloc(N*B->length*sizeof(ulong));
        mpoly_repack_monomials(Bexp, Abits, B->exps, B->bits,
                                                        B->length, ctx->minfo);
    }

    freeCexp = 0;
    Cexp = C->exps;
    if (Abits > C->bits)
    {
        freeCexp = 1;
        Cexp = (ulong *) flint_malloc(N*C->length*sizeof(ulong));
        mpoly_repack_monomials(Cexp, Abits, C->exps, C->bits,
                                                        C->length, ctx->minfo);
    }

    /* deal with aliasing and do multiplication */
    if (A == B || A == C)
    {
        fmpz_mpoly_t T;
        fmpz_mpoly_init3(T, B->length + C->length, Abits, ctx);

        /* algorithm more efficient if smaller poly first */
        if (B->length > C->length)
        {
            Alen = _fmpz_mpoly_mul_johnson(&T->coeffs, &T->exps, &T->alloc,
                                                  C->coeffs, Cexp, C->length,
                                                  B->coeffs, Bexp, B->length,
                                                         Abits, N, cmpmask);
        }
        else
        {
            Alen = _fmpz_mpoly_mul_johnson(&T->coeffs, &T->exps, &T->alloc,
                                                  B->coeffs, Bexp, B->length,
                                                  C->coeffs, Cexp, C->length,
                                                         Abits, N, cmpmask);
        }

        fmpz_mpoly_swap(T, A, ctx);
        fmpz_mpoly_clear(T, ctx);
    }
    else
    {
        fmpz_mpoly_fit_length_reset_bits(A, B->length + C->length, Abits, ctx);

        /* algorithm more efficient if smaller poly first */
        if (B->length > C->length)
        {
            Alen = _fmpz_mpoly_mul_johnson(&A->coeffs, &A->exps, &A->alloc,
                                                  C->coeffs, Cexp, C->length,
                                                  B->coeffs, Bexp, B->length,
                                                         Abits, N, cmpmask);
        }
        else
        {
            Alen = _fmpz_mpoly_mul_johnson(&A->coeffs, &A->exps, &A->alloc,
                                                  B->coeffs, Bexp, B->length,
                                                  C->coeffs, Cexp, C->length,
                                                         Abits, N, cmpmask);
        }
    }

    if (freeBexp)
        flint_free(Bexp);

    if (freeCexp)
        flint_free(Cexp);

    _fmpz_mpoly_set_length(A, Alen, ctx);

    TMP_END;
}



void fmpz_mpoly_mul_johnson(
    fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_t C,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    fmpz * maxBfields, * maxCfields;
    TMP_INIT;

    if (B->length == 0 || C->length == 0)
    {
        fmpz_mpoly_zero(A, ctx);
        return;
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

    _fmpz_mpoly_mul_johnson_maxfields(A, B, maxBfields, C, maxCfields, ctx);

    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_clear(maxBfields + i);
        fmpz_clear(maxCfields + i);
    }

    TMP_END;
}

/* TODO: decouple the exp/coeff alloc in fmpz_mpoly and move this to mpoly */
void fmpz_mpoly_mul_monomial(
    fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_t C,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, N, Blen = B->length;
    ulong ofmask;
    flint_bitcnt_t Abits;
    ulong * Aexps, * Bexps = B->exps, * Cexps = C->exps;
    fmpz Ccoeff0 = C->coeffs[0];
    int freeCcoeff0 = 0, overflowed = 0;
    TMP_INIT;

    FLINT_ASSERT(C->length == 1);

    if (A == C)
    {
        freeCcoeff0 = 1;
        fmpz_init_set(&Ccoeff0, C->coeffs + 0);
    }

    if (C->exps[0] == 0 && mpoly_monomial_is_zero(C->exps,
                                     mpoly_words_per_exp(C->bits, ctx->minfo)))
    {
        fmpz_mpoly_scalar_mul_fmpz(A, B, &Ccoeff0, ctx);
        goto cleanup_C;
    }

    TMP_START;

    Abits = FLINT_MAX(B->bits, C->bits);
    N = mpoly_words_per_exp(Abits, ctx->minfo);

    if (A == C || Abits != C->bits)
    {
        Cexps = TMP_ARRAY_ALLOC(N, ulong);
        mpoly_repack_monomials(Cexps, Abits, C->exps, C->bits, 1, ctx->minfo);
    }

    if (A == B)
    {
        /* inplace operation on A */
        fmpz_mpoly_fit_bits(A, Abits, ctx);
        Bexps = Aexps = A->exps;
    }
    else
    {
        if (Abits != B->bits)
        {
            Bexps = TMP_ARRAY_ALLOC(N*Blen, ulong);
            mpoly_repack_monomials(Bexps, Abits, B->exps, B->bits, Blen, ctx->minfo);
        }

        fmpz_mpoly_fit_length_reset_bits(A, Blen, Abits, ctx);
        Aexps = A->exps;
    }

    if (Abits > FLINT_BITS)
    {
        for (i = 0; i < Blen; i++)
            mpoly_monomial_add_mp(Aexps + N*i, Bexps + N*i, Cexps + N*0, N);

        for (i = 0; !overflowed && i < Blen; i++)
            overflowed = mpoly_monomial_overflows_mp(Aexps + N*i, N, Abits);
    }
    else
    {
        for (i = 0; i < Blen; i++)
            mpoly_monomial_add(Aexps + N*i, Bexps + N*i, Cexps + N*0, N);

        ofmask = mpoly_overflow_mask_sp(Abits);
        for (i = 0; !overflowed && i < Blen; i++)
            overflowed = mpoly_monomial_overflows(Aexps + N*i, N, ofmask);
    }

    TMP_END;

    /* slightly dirty: repack monomials can handle 1-bit overfown fields */
    if (overflowed)
    {
        ulong * newAexps;
        flint_bitcnt_t newAbits = mpoly_fix_bits(Abits + 1, ctx->minfo);
        N = mpoly_words_per_exp(newAbits, ctx->minfo);
        newAexps = FLINT_ARRAY_ALLOC(N*A->alloc, ulong);
        mpoly_repack_monomials(newAexps, newAbits, A->exps, Abits, Blen, ctx->minfo);
        flint_free(A->exps);
        A->exps = newAexps;
        A->bits = newAbits;
    }

    _fmpz_vec_scalar_mul_fmpz(A->coeffs, B->coeffs, Blen, &Ccoeff0);
    _fmpz_mpoly_set_length(A, Blen, ctx);

cleanup_C:

    if (freeCcoeff0)
        fmpz_clear(&Ccoeff0);
}
