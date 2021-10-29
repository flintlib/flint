/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

/* improve locality */
#define BLOCK 128
#define MAX_ARRAY_SIZE (WORD(300000))
#define MAX_LEX_SIZE (WORD(300))



void _nmod_mpoly_addmul_array1_ulong1(ulong * poly1, 
                           const ulong * poly2, const ulong * exp2, slong len2,
                           const ulong * poly3, const ulong * exp3, slong len3)
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

void _nmod_mpoly_addmul_array1_ulong2(ulong * poly1, 
                           const ulong * poly2, const ulong * exp2, slong len2,
                           const ulong * poly3, const ulong * exp3, slong len3)
{
   slong ii, i, jj, j;
   ulong p[2];
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
                  umul_ppmm(p[1], p[0], poly2[i], poly3[j]);
                  add_ssaaaa(c[1], c[0], c[1], c[0], p[1], p[0]);
               }
            }
         }
      }
   }
}

void _nmod_mpoly_addmul_array1_ulong3(ulong * poly1, 
                           const ulong * poly2, const ulong * exp2, slong len2,
                           const ulong * poly3, const ulong * exp3, slong len3)
{
   slong ii, i, jj, j;
   ulong p[2];
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
                  umul_ppmm(p[1], p[0], poly2[i], poly3[j]);
                  add_sssaaaaaa(c[2], c[1], c[0], c[2], c[1], c[0], 0, p[1], p[0]);
               }
            }
         }
      }
   }
}


/****************************************************
    LEX
****************************************************/


#define LEX_UNPACK_MACRO(fxn_name, coeff_decl, nonzero_test, reduce_coeff)     \
slong fxn_name(nmod_mpoly_t P, slong Plen, coeff_decl,                         \
             const ulong * mults, slong num, slong array_size, slong top,      \
                                            const nmod_mpoly_ctx_t ctx)        \
{                                                                              \
    slong off, j;                                                              \
    slong topmult = num == 0 ? 1 : mults[num - 1];                             \
    slong lastd   = topmult - 1;                                               \
    slong reset   = array_size/topmult;                                        \
    slong counter = reset;                                                     \
    ulong startexp = (top << (P->bits*num)) + (lastd << (P->bits*(num-1)));    \
    ulong coeff;                                                               \
    for (off = array_size - 1; off >= 0; off--)                                \
    {                                                                          \
        if (nonzero_test)                                                      \
        {                                                                      \
            reduce_coeff                                                       \
            if (coeff != UWORD(0))                                             \
            {                                                                  \
                slong d = off;                                                 \
                ulong exp = startexp;                                          \
                for (j = 0; j + 1 < num; j++) {                                \
                    exp += (d % mults[j]) << (P->bits*j);                      \
                    d = d / mults[j];                                          \
                }                                                              \
                _nmod_mpoly_fit_length(&P->coeffs, &P->coeffs_alloc,           \
                                       &P->exps, &P->exps_alloc, 1, Plen + 1); \
                P->exps[Plen] = exp;                                           \
                P->coeffs[Plen] = coeff;                                       \
                Plen++;                                                        \
            }                                                                  \
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
*/

LEX_UNPACK_MACRO(
    nmod_mpoly_append_array_sm1_LEX, ulong * coeff_array
,
    coeff_array[off] != UWORD(0)
,
    NMOD_RED(coeff, coeff_array[off], ctx->mod);
    coeff_array[off] = 0;
)

LEX_UNPACK_MACRO(
    nmod_mpoly_append_array_sm2_LEX, ulong * coeff_array
,
    (coeff_array[2*off + 0] || coeff_array[2*off + 1]) != UWORD(0)
,
    NMOD2_RED2(coeff, coeff_array[2*off + 1], coeff_array[2*off + 0], ctx->mod);
    coeff_array[2*off + 0] = coeff_array[2*off + 1] = 0;
)

LEX_UNPACK_MACRO(
    nmod_mpoly_append_array_sm3_LEX, ulong * coeff_array
,
    (coeff_array[3*off + 0] || coeff_array[3*off + 1] || coeff_array[3*off + 2]) != UWORD(0)
,
    NMOD_RED3(coeff, coeff_array[3*off + 2], coeff_array[3*off + 1], coeff_array[3*off + 0], ctx->mod);
    coeff_array[3*off + 0] = coeff_array[3*off + 1] = coeff_array[3*off + 2] = UWORD(0);
)



void _nmod_mpoly_mul_array_chunked_LEX(
    nmod_mpoly_t P,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const ulong * mults,
    const nmod_mpoly_ctx_t ctx)
{
    slong num = ctx->minfo->nfields - 1;
    slong Pi, i, j, Plen, Pl, Al, Bl, array_size;
    slong * Amain, * Bmain;
    ulong * Apexp, * Bpexp;
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
    Amain = (slong *) TMP_ALLOC((Al + 1)*sizeof(slong));
    Bmain = (slong *) TMP_ALLOC((Bl + 1)*sizeof(slong));
    Apexp = (ulong *) flint_malloc(A->length*sizeof(ulong));
    Bpexp = (ulong *) flint_malloc(B->length*sizeof(ulong));
    mpoly_main_variable_split_LEX(Amain, Apexp, A->exps, Al, A->length, mults, num, A->bits);
    mpoly_main_variable_split_LEX(Bmain, Bpexp, B->exps, Bl, B->length, mults, num, B->bits);

    Pl = Al + Bl - 1;
    Plen = 0;

    {
        ulong * coeff_array = (ulong *) TMP_ALLOC(3*array_size*sizeof(ulong));
        for (j = 0; j < 3*array_size; j++)
            coeff_array[j] = 0;

        /* for each output chunk */
        for (Pi = 0; Pi < Pl; Pi++)
        {
            /* compute bound on coeffs of output chunk */
            slong len = 0;
            mp_limb_t t2, t1, t0, u1, u0;

            for (i = 0, j = Pi; i < Al && j >= 0; i++, j--)
            {
                if (j < Bl)
                {
                    len += FLINT_MIN(Amain[i + 1] - Amain[i],
                                     Bmain[j + 1] - Bmain[j]);
                }
            }

            umul_ppmm(t1, t0, ctx->mod.n - 1, ctx->mod.n - 1);
            umul_ppmm(t2, t1, t1, len);
            umul_ppmm(u1, u0, t0, len);
            add_sssaaaaaa(t2, t1, t0,  t2, t1, UWORD(0),  UWORD(0), u1, u0);

            if (t2 != UWORD(0))
            {
                /* need three words */
                for (i = 0, j = Pi; i < Al && j >= 0; i++, j--)
                {
                    if (j < Bl)
                    {
                        _nmod_mpoly_addmul_array1_ulong3(coeff_array, 
                                A->coeffs + Amain[i],
                                Apexp + Amain[i], Amain[i + 1] - Amain[i],
                                B->coeffs + Bmain[j],
                                Bpexp + Bmain[j], Bmain[j + 1] - Bmain[j]);
                    }
                }
                Plen = nmod_mpoly_append_array_sm3_LEX(P, Plen, coeff_array,
                                      mults, num, array_size, Pl - Pi - 1, ctx);

            } else if (t1 != UWORD(0))
            {
                /* fits into two words */
                for (i = 0, j = Pi; i < Al && j >= 0; i++, j--)
                {
                    if (j < Bl)
                    {
                        _nmod_mpoly_addmul_array1_ulong2(coeff_array, 
                                A->coeffs + Amain[i],
                                Apexp + Amain[i], Amain[i + 1] - Amain[i],
                                B->coeffs + Bmain[j],
                                Bpexp + Bmain[j], Bmain[j + 1] - Bmain[j]);

                    }
                }
                Plen = nmod_mpoly_append_array_sm2_LEX(P, Plen, coeff_array,
                                      mults, num, array_size, Pl - Pi - 1, ctx);

            } else
            {
                /* fits into one word */
                for (i = 0, j = Pi; i < Al && j >= 0; i++, j--)
                {
                    if (j < Bl)
                    {
                        _nmod_mpoly_addmul_array1_ulong1(coeff_array, 
                                A->coeffs + Amain[i],
                                Apexp + Amain[i], Amain[i + 1] - Amain[i],
                                B->coeffs + Bmain[j],
                                Bpexp + Bmain[j], Bmain[j + 1] - Bmain[j]);
                    }
                }
                Plen = nmod_mpoly_append_array_sm1_LEX(P, Plen, coeff_array,
                                      mults, num, array_size, Pl - Pi - 1, ctx);
            }
        }
    }

    _nmod_mpoly_set_length(P, Plen, ctx);

    flint_free(Apexp);
    flint_free(Bpexp);
    TMP_END;
}



int _nmod_mpoly_mul_array_LEX(
    nmod_mpoly_t A,
    const nmod_mpoly_t B,
    fmpz * maxBfields,
    const nmod_mpoly_t C,
    fmpz * maxCfields,
    const nmod_mpoly_ctx_t ctx)
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

    /* array multiplication assumes result fit into 1 word */
    if (1 != mpoly_words_per_exp(exp_bits, ctx->minfo))
    {
        success = 0;
        goto cleanup;
    }

    /* handle aliasing and do array multiplication */
    if (A == B || A == C)
    {
        nmod_mpoly_t T;
        nmod_mpoly_init3(T, B->length + C->length - 1, exp_bits, ctx);
        _nmod_mpoly_mul_array_chunked_LEX(T, C, B, mults, ctx);
        nmod_mpoly_swap(T, A, ctx);
        nmod_mpoly_clear(T, ctx);
    }
    else
    {
        nmod_mpoly_fit_length_reset_bits(A, B->length + C->length - 1, exp_bits, ctx);
        _nmod_mpoly_mul_array_chunked_LEX(A, C, B, mults, ctx);
    }
    success = 1;

cleanup:

    TMP_END;

    return success;
}




/****************************************************
    DEGLEX and DEGREVLEX
****************************************************/


#define DEGLEX_UNPACK_MACRO(fxn_name, coeff_decl, nonzero_test, reduce_coeff)  \
slong fxn_name(nmod_mpoly_t P, slong Plen, coeff_decl,                         \
                         slong top, slong nvars, slong degb,                   \
                                             const nmod_mpoly_ctx_t ctx)       \
{                                                                              \
    slong i;                                                                   \
    ulong exp, lomask = (UWORD(1) << (P->bits - 1)) - 1;                       \
    slong off, array_size;                                                     \
    slong * curexp, * degpow;                                                  \
    ulong * oneexp;                                                            \
    ulong coeff;                                                               \
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
            reduce_coeff                                                       \
            if (coeff != UWORD(0))                                             \
            {                                                                  \
                _nmod_mpoly_fit_length(&P->coeffs, &P->coeffs_alloc,           \
                                       &P->exps, &P->exps_alloc, 1, Plen + 1); \
                P->exps[Plen] = exp;                                           \
                P->coeffs[Plen] = coeff;                                       \
                Plen++;                                                        \
            }                                                                  \
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
    nmod_mpoly_append_array_sm1_DEGLEX, ulong * coeff_array
,
    coeff_array[off] != UWORD(0)
,
    NMOD_RED(coeff, coeff_array[off], ctx->mod);
    coeff_array[off] = 0;
)

DEGLEX_UNPACK_MACRO(
    nmod_mpoly_append_array_sm2_DEGLEX, ulong * coeff_array
,
    (coeff_array[2*off + 0] || coeff_array[2*off + 1]) != UWORD(0)
,
    NMOD2_RED2(coeff, coeff_array[2*off + 1], coeff_array[2*off + 0], ctx->mod);
    coeff_array[2*off + 0] = coeff_array[2*off + 1] = 0;
)

DEGLEX_UNPACK_MACRO(
    nmod_mpoly_append_array_sm3_DEGLEX, ulong * coeff_array
,
    (coeff_array[3*off + 0] || coeff_array[3*off + 1] || coeff_array[3*off + 2]) != UWORD(0)
,
    NMOD_RED3(coeff, coeff_array[3*off + 2], coeff_array[3*off + 1], coeff_array[3*off + 0], ctx->mod);
    coeff_array[3*off + 0] = coeff_array[3*off + 1] = coeff_array[3*off + 2] = 0;
)





#define DEGREVLEX_UNPACK_MACRO(fxn_name, coeff_decl, nonzero_test, reduce_coeff) \
slong fxn_name(nmod_mpoly_t P, slong Plen, coeff_decl,                         \
                           slong top, slong nvars, slong degb,                 \
                                            const nmod_mpoly_ctx_t ctx)        \
{                                                                              \
    slong i;                                                                   \
    ulong exp, mask = UWORD(1) << (P->bits - 1);                               \
    slong off, array_size;                                                     \
    slong * curexp, * degpow;                                                  \
    ulong * oneexp;                                                            \
    ulong coeff;                                                               \
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
            reduce_coeff                                                       \
            if (coeff != UWORD(0))                                             \
            {                                                                  \
                _nmod_mpoly_fit_length(&P->coeffs, &P->coeffs_alloc,           \
                                       &P->exps, &P->exps_alloc, 1, Plen + 1); \
                P->exps[Plen] = exp;                                           \
                P->coeffs[Plen] = coeff;                                       \
                Plen++;                                                        \
            }                                                                  \
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
    nmod_mpoly_append_array_sm1_DEGREVLEX, ulong * coeff_array
,
    coeff_array[off] != WORD(0)
,
    NMOD_RED(coeff, coeff_array[off], ctx->mod);
    coeff_array[off] = 0;
)

DEGREVLEX_UNPACK_MACRO(
    nmod_mpoly_append_array_sm2_DEGREVLEX, ulong * coeff_array
,
    (coeff_array[2*off + 0] || coeff_array[2*off + 1]) != WORD(0)
,
    NMOD2_RED2(coeff, coeff_array[2*off + 1], coeff_array[2*off + 0], ctx->mod);
    coeff_array[2*off + 0] = coeff_array[2*off + 1] = 0;
)

DEGREVLEX_UNPACK_MACRO(
    nmod_mpoly_append_array_sm3_DEGREVLEX, ulong * coeff_array
,
    (coeff_array[3*off + 0] || coeff_array[3*off + 1] || coeff_array[3*off + 2]) != WORD(0)
,
    NMOD_RED3(coeff, coeff_array[3*off + 2], coeff_array[3*off + 1], coeff_array[3*off + 0], ctx->mod);
    coeff_array[3*off + 0] = coeff_array[3*off + 1] = coeff_array[3*off + 2] = 0;
)


void _nmod_mpoly_mul_array_chunked_DEG(
    nmod_mpoly_t P,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    ulong degb,
    const nmod_mpoly_ctx_t ctx)
{
    slong nvars = ctx->minfo->nvars;
    slong Pi, i, j, Plen, Pl, Al, Bl, array_size;
    slong * Amain, * Bmain;
    ulong * Apexp, * Bpexp;
    slong (* upack_sm1)(nmod_mpoly_t, slong, ulong *, slong, slong, slong, const nmod_mpoly_ctx_t);
    slong (* upack_sm2)(nmod_mpoly_t, slong, ulong *, slong, slong, slong, const nmod_mpoly_ctx_t);
    slong (* upack_sm3)(nmod_mpoly_t, slong, ulong *, slong, slong, slong, const nmod_mpoly_ctx_t);
    TMP_INIT;

    TMP_START;

    /* compute lengths of poly2 and poly3 in chunks */
    Al = 1 + (slong) (A->exps[0] >> (A->bits*nvars));
    Bl = 1 + (slong) (B->exps[0] >> (B->bits*nvars));

    array_size = 1;
    for (i = 0; i < nvars-1; i++) {
        array_size *= degb;
    }

    upack_sm1  = &nmod_mpoly_append_array_sm1_DEGLEX;
    upack_sm2  = &nmod_mpoly_append_array_sm2_DEGLEX;
    upack_sm3  = &nmod_mpoly_append_array_sm3_DEGLEX;
    if (ctx->minfo->ord == ORD_DEGREVLEX) {
        upack_sm1  = &nmod_mpoly_append_array_sm1_DEGREVLEX;
        upack_sm2  = &nmod_mpoly_append_array_sm2_DEGREVLEX;
        upack_sm3  = &nmod_mpoly_append_array_sm3_DEGREVLEX;
    }

    /* compute indices and lengths of coefficients of polys in main variable */
    Amain = (slong *) TMP_ALLOC((Al + 1)*sizeof(slong));
    Bmain = (slong *) TMP_ALLOC((Bl + 1)*sizeof(slong));
    Apexp = (ulong *) flint_malloc(A->length*sizeof(ulong));
    Bpexp = (ulong *) flint_malloc(B->length*sizeof(ulong));
    mpoly_main_variable_split_DEG(Amain, Apexp, A->exps, Al, A->length,
                                                         degb, nvars, A->bits);
    mpoly_main_variable_split_DEG(Bmain, Bpexp, B->exps, Bl, B->length,
                                                         degb, nvars, B->bits);

    Pl = Al + Bl - 1;
    FLINT_ASSERT(Pl == degb);
    Plen = 0;

    {
        ulong * coeff_array = (ulong *) TMP_ALLOC(3*array_size*sizeof(ulong));
        for (j = 0; j < 3*array_size; j++)
            coeff_array[j] = 0;

        /* for each output chunk */
        for (Pi = 0; Pi < Pl; Pi++)
        {
            /* compute bound on coeffs of output chunk */
            slong len = 0;
            mp_limb_t t2, t1, t0, u1, u0;

            for (i = 0, j = Pi; i < Al && j >= 0; i++, j--)
            {
                if (j < Bl)
                {
                    len += FLINT_MIN(Amain[i + 1] - Amain[i],
                                     Bmain[j + 1] - Bmain[j]);
                }
            }

            umul_ppmm(t1, t0, ctx->mod.n - 1, ctx->mod.n - 1);
            umul_ppmm(t2, t1, t1, len);
            umul_ppmm(u1, u0, t0, len);
            add_sssaaaaaa(t2, t1, t0,  t2, t1, UWORD(0),  UWORD(0), u1, u0);

            if (t2 != UWORD(0))
            {
                /* need three words */
                for (i = 0, j = Pi; i < Al && j >= 0; i++, j--)
                {
                    if (j < Bl)
                    {
                        _nmod_mpoly_addmul_array1_ulong3(coeff_array, 
                                A->coeffs + Amain[i],
                                Apexp + Amain[i], Amain[i + 1] - Amain[i],
                                B->coeffs + Bmain[j],
                                Bpexp + Bmain[j], Bmain[j + 1] - Bmain[j]);
                    }
                }

                Plen = upack_sm3(P, Plen, coeff_array, Pl - Pi - 1, nvars, degb, ctx);

            } else if (t1 != UWORD(0))
            {
                /* fits into two words */
                for (i = 0, j = Pi; i < Al && j >= 0; i++, j--)
                {
                    if (j < Bl)
                    {
                        _nmod_mpoly_addmul_array1_ulong2(coeff_array, 
                                A->coeffs + Amain[i],
                                Apexp + Amain[i], Amain[i + 1] - Amain[i],
                                B->coeffs + Bmain[j],
                                Bpexp + Bmain[j], Bmain[j + 1] - Bmain[j]);

                    }
                }
                Plen = upack_sm2(P, Plen, coeff_array, Pl - Pi - 1, nvars, degb, ctx);

            } else if (t0 != UWORD(0))
            {
                /* fits into one word */
                for (i = 0, j = Pi; i < Al && j >= 0; i++, j--)
                {
                    if (j < Bl)
                    {
                        _nmod_mpoly_addmul_array1_ulong1(coeff_array, 
                                A->coeffs + Amain[i],
                                Apexp + Amain[i], Amain[i + 1] - Amain[i],
                                B->coeffs + Bmain[j],
                                Bpexp + Bmain[j], Bmain[j + 1] - Bmain[j]);
                    }
                }
                Plen = upack_sm1(P, Plen, coeff_array, Pl - Pi - 1, nvars, degb, ctx);
            }
        }
    }

    _nmod_mpoly_set_length(P, Plen, ctx);

    flint_free(Apexp);
    flint_free(Bpexp);
    TMP_END;
}



int _nmod_mpoly_mul_array_DEG(
    nmod_mpoly_t A,
    const nmod_mpoly_t B, fmpz * maxBfields,
    const nmod_mpoly_t C, fmpz * maxCfields,
    const nmod_mpoly_ctx_t ctx)
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


    /* the field of index n-1 is the one that wil be pulled out */
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
        nmod_mpoly_t T;
        nmod_mpoly_init3(T, B->length + C->length - 1, exp_bits, ctx);
        _nmod_mpoly_mul_array_chunked_DEG(T, C, B, deg, ctx);
        nmod_mpoly_swap(T, A, ctx);
        nmod_mpoly_clear(T, ctx);
    }
    else
    {
        nmod_mpoly_fit_length_reset_bits(A, B->length + C->length - 1, exp_bits, ctx);
        _nmod_mpoly_mul_array_chunked_DEG(A, C, B, deg, ctx);
    }
    success = 1;

cleanup:

    return success;
}



int nmod_mpoly_mul_array(nmod_mpoly_t A, const nmod_mpoly_t B,
                              const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    int success;
    fmpz * maxBfields, * maxCfields;
    TMP_INIT;

    if (B->length == 0 || C->length == 0)
    {
        nmod_mpoly_zero(A, ctx);
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
            success = _nmod_mpoly_mul_array_LEX(A, B, maxBfields,
                                                   C, maxCfields, ctx);
            break;
        }
        case ORD_DEGLEX:
        case ORD_DEGREVLEX:
        {
            success = _nmod_mpoly_mul_array_DEG(A, B, maxBfields,
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
