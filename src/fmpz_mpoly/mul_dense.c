/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "fmpz_mpoly.h"

static slong fmpz_mpolyd_length(const fmpz_mpolyd_t A)
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


static int fmpz_mpolyd_set_degbounds(fmpz_mpolyd_t A, slong * bounds)
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
static void fmpz_mpoly_convert_to_fmpz_mpolyd_degbound(fmpz_mpolyd_t poly1,
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
static void fmpz_mpoly_consume_fmpz_mpolyd_clear(fmpz_mpoly_t A, fmpz_mpolyd_t B,
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
            msb = flint_clz(topmask);
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
