/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"
#include "nmod_mpoly.h"
#include "assert.h"
#include "unistd.h"

void fmpz_mpolyd_swap(fmpz_mpolyd_t A, fmpz_mpolyd_t B)
{
    fmpz_mpolyd_struct t = *A;
    *A = *B;
    *B = t;
}

void fmpz_mpolyd_ctx_init(fmpz_mpolyd_ctx_t dctx, slong nvars)
{
    slong i;

    dctx->nvars = nvars;
    dctx->perm = (slong *) flint_malloc(nvars*sizeof(slong));
    for (i = 0; i < nvars; i++)
    {
        dctx->perm[i] = i;
    }
}


int fmpz_mpolyd_ctx_settle(fmpz_mpolyd_ctx_t dctx,
                            const fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong i, j, degb_prod;
    slong * Aexps, * Bexps, * deg_bounds;
    slong nvars = ctx->minfo->nvars;
    slong * perm = dctx->perm;
    TMP_INIT;

    success = 0;

    TMP_START;
    Aexps = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    Bexps = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    if (A->bits > FLINT_BITS || B->bits > FLINT_BITS)
        goto cleanup;
    fmpz_mpoly_degrees_si(Aexps, A, ctx);
    fmpz_mpoly_degrees_si(Bexps, B, ctx);

    deg_bounds = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    for (i = 0; i < nvars; i++)
    {
        dctx->perm[i] = i;
    }

    degb_prod = 1;
    for (i = 0; i < nvars; i++)
    {
        ulong hi;
        deg_bounds[i] = FLINT_MAX(Aexps[i] + 1, Bexps[i] + 1);
        umul_ppmm(hi, degb_prod, degb_prod, deg_bounds[i]);
        if (hi != WORD(0) || degb_prod < 0)
            goto cleanup;
    }

    success = 1;
    for (i = 1; i < nvars; i++)
    {
        for (j = i; (j > 0) && (deg_bounds[j-1] < deg_bounds[j-0]); j--)
        {
            slong t1, t2;
            t1 = deg_bounds[j-1];
            t2 = deg_bounds[j-0];
            deg_bounds[j-0] = t1;
            deg_bounds[j-1] = t2;
            t1 = perm[j-1];
            t2 = perm[j-0];
            perm[j-0] = t1;
            perm[j-1] = t2;
        }
    }

cleanup:
    TMP_END;
    return success;
}


void fmpz_mpolyd_ctx_clear(fmpz_mpolyd_ctx_t dctx)
{
    flint_free(dctx->perm);
}


void fmpz_mpolyd_init(fmpz_mpolyd_t poly, slong nvars)
{
    slong i;

    poly->nvars = nvars;
    poly->degb_alloc = nvars;
    poly->deg_bounds = (slong *) flint_malloc(poly->degb_alloc*sizeof(slong));
    for (i = 0; i < nvars; i++)
    {
        poly->deg_bounds[i] = WORD(1);
    }
    poly->coeff_alloc = WORD(16);
    poly->coeffs = (fmpz *) flint_malloc(poly->coeff_alloc*sizeof(fmpz));
    for (i = 0; i < poly->coeff_alloc; i++)
    {
        fmpz_init(poly->coeffs + i);
    }
}


void fmpz_mpolyd_fit_length(fmpz_mpolyd_t poly, slong len) {
    if (poly->coeff_alloc < len) {
        slong i;
        poly->coeffs = (fmpz *) flint_realloc(poly->coeffs, len*sizeof(fmpz));
        for (i = poly->coeff_alloc; i < len; i++)
        {
            fmpz_init(poly->coeffs + i);
        }
        poly->coeff_alloc = len;
    }
}


void fmpz_mpolyd_set_nvars(fmpz_mpolyd_t poly, slong nvars) {

    poly->nvars = nvars;
    if (poly->degb_alloc < nvars) {
        poly->deg_bounds = (slong *) flint_realloc(poly->deg_bounds, nvars*sizeof(slong));
        poly->degb_alloc = nvars;
    }
}


void fmpz_mpolyd_zero(fmpz_mpolyd_t poly)
{
    slong i;

    for (i = 0; i < poly->nvars; i++)
    {
        poly->deg_bounds[i] = WORD(1);
    }
    poly->coeffs[0] = UWORD(0);
}

void fmpz_mpolyd_set_fmpz(fmpz_mpolyd_t poly, fmpz_t num)
{
    slong i;

    for (i = 0; i < poly->nvars; i++)
    {
        poly->deg_bounds[i] = WORD(1);
    }
    fmpz_set(poly->coeffs + 0, num);
}


void fmpz_mpolyd_clear(fmpz_mpolyd_t poly)
{
    slong i;

    for (i = 0; i < poly->coeff_alloc; i++)
    {
        fmpz_clear(poly->coeffs + i);
    }

    flint_free(poly->deg_bounds);
    flint_free(poly->coeffs);
    poly->deg_bounds = NULL;
    poly->coeffs = NULL;
}


void fmpz_mpoly_convert_to_fmpz_mpolyd(
                            fmpz_mpolyd_t poly1, const fmpz_mpolyd_ctx_t dctx,
                          const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx)
{
    slong degb_prod;
    slong i, j, N;
    slong * exps;
    const slong * perm = dctx->perm;
    slong nvars = ctx->minfo->nvars;
    TMP_INIT;

    fmpz_mpolyd_set_nvars(poly1, ctx->minfo->nvars);

    FLINT_ASSERT(poly2->bits <= FLINT_BITS)

    if (poly2->length == 0)
    {
        fmpz_mpolyd_zero(poly1);
        return;
    }

    TMP_START;
    exps = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));

    fmpz_mpoly_degrees_si(exps, poly2, ctx);
    degb_prod = WORD(1);
    for (i = 0; i < nvars; i++)
    {
        poly1->deg_bounds[i] = exps[perm[i]] + 1;
        degb_prod *= poly1->deg_bounds[i];
    }

    fmpz_mpolyd_fit_length(poly1, degb_prod);
    for (i = 0; i < degb_prod; i++)
    {
        fmpz_zero(poly1->coeffs + i);
    }

    N = mpoly_words_per_exp(poly2->bits, ctx->minfo);
    for (i = 0; i < poly2->length; i++)
    {
        slong off = 0;

        mpoly_get_monomial_ui((ulong *)exps, poly2->exps + N*i, poly2->bits, ctx->minfo);
        for (j = 0; j < nvars; j++)
        {
            off = exps[perm[j]] + poly1->deg_bounds[j]*off;
        }

        fmpz_set(poly1->coeffs + off, poly2->coeffs + i);
    }

    TMP_END;
}

void fmpz_mpolyd_print(fmpz_mpolyd_t poly, const char ** vars,
                                                  const fmpz_mpolyd_ctx_t dctx)
{
    int first = 0;
    slong i, j;
    slong degb_prod;

    degb_prod = WORD(1);
    for (j = 0; j < poly->nvars; j++) {
        degb_prod *= poly->deg_bounds[j];
    }

    first = 1;
    for (i = 0; i < degb_prod; i++) {
        ulong k = i;

        if (fmpz_is_zero(poly->coeffs + i))
            continue;

        if (!first)
            printf(" + ");

        fmpz_print(poly->coeffs + i);

        for (j = poly->nvars - 1; j >= 0; j--)
        {
            ulong m = poly->deg_bounds[j];
            ulong e = k % m;
            k = k / m;
            flint_printf("*%s^%wd", vars[dctx->perm[j]], e);
        }
        FLINT_ASSERT(k == 0);
        first = 0;
    }

    if (first)
        flint_printf("0");
}


void fmpz_mpolyd_print_simple(fmpz_mpolyd_t poly)
{
    int first = 0;
    slong i, j;
    slong degb_prod;

    degb_prod = WORD(1);
    for (j = 0; j < poly->nvars; j++) {
        degb_prod *= poly->deg_bounds[j];
    }

    first = 1;
    for (i = 0; i < degb_prod; i++) {
        ulong k = i;

        if (fmpz_is_zero(poly->coeffs + i))
            continue;

        if (!first)
            printf(" + ");

        fmpz_print(poly->coeffs + i);

        for (j = poly->nvars - 1; j >= 0; j--)
        {
            ulong m = poly->deg_bounds[j];
            ulong e = k % m;
            k = k / m;
            flint_printf("*x%d^%wd", j, e);
        }
        FLINT_ASSERT(k == 0);
        first = 0;
    }

    if (first)
        flint_printf("0");
}


void fmpz_mpoly_convert_from_fmpz_mpolyd(
                                  fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx,
                           const fmpz_mpolyd_t B, const fmpz_mpolyd_ctx_t dctx)
{
    slong i, j;
    slong degb_prod;
    slong * perm = dctx->perm;
    ulong * exps;
    TMP_INIT;

    FLINT_ASSERT(ctx->minfo->nvars == B->nvars);

    degb_prod = WORD(1);
    for (j = 0; j < B->nvars; j++) {
        degb_prod *= B->deg_bounds[j];
    }

    TMP_START;
    exps = (ulong *) TMP_ALLOC(B->nvars*sizeof(ulong));

    fmpz_mpoly_zero(A, ctx);
    for (i = 0; i < degb_prod; i++) {
        ulong k = i;

        if (fmpz_is_zero(B->coeffs + i))
            continue;

        for (j = B->nvars - 1; j >= 0; j--) 
        {
            ulong m = B->deg_bounds[j];
            ulong e = k % m;
            k = k / m;
            exps[perm[j]] = e;
        }
        FLINT_ASSERT(k == 0);

        fmpz_mpoly_set_term_fmpz_ui(A, B->coeffs + i, exps, ctx);
    }

    TMP_END;
}



/*
void fmpz_CRT_ui(fmpz_t out, const fmpz_t r1, const fmpz_t m1,
    ulong r2, ulong m2, int sign)
*/

void fmpz_mpolyd_CRT_nmod(fmpz_mpolyd_t A,
                                 const fmpz_mpolyd_t B, const fmpz_t Bm,
                                 const nmod_mpolyd_t C, const nmodf_ctx_t fctx)
{
    int Bok, Cok;
    slong carry;
    slong Bind, Cind;
    slong i, j;
    slong * inds;
    slong nvars = B->nvars;
    slong degb_prod;
    ulong hi;
    fmpz_t zero;
    TMP_INIT;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(B->nvars == C->nvars);

    fmpz_mpolyd_set_nvars(A, B->nvars);

/*printf("CRT here0\n");*/

    degb_prod = 1;
    for (j = 0; j < nvars; j++)
    {
        A->deg_bounds[j] = FLINT_MAX(B->deg_bounds[j], C->deg_bounds[j]);
        umul_ppmm(hi, degb_prod, degb_prod, A->deg_bounds[j]);
        if (hi != WORD(0) || degb_prod < 0)
            assert(0);
    }
/*printf("CRT here1\n");*/

    fmpz_mpolyd_fit_length(A, degb_prod);

/*printf("CRT here2\n");*/

    fmpz_init_set_ui(zero, 0);


    TMP_START;
    inds = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    for (j = 0; j < nvars; j++)
        inds[j] = 0;
    Bok = 1;
    Cok = 1;
    Bind = 0;
    Cind = 0;
    for (i = 0; i < degb_prod; i++)
    {
/*printf("CRT i: %d\n",i);*/

               if (Bok && Cok) {
            fmpz_CRT_ui(A->coeffs + i, B->coeffs + Bind++, Bm, C->coeffs[Cind++], fctx->mod.n, 1);
        } else if (Bok && !Cok) {
            fmpz_CRT_ui(A->coeffs + i, B->coeffs + Bind++, Bm,    0             , fctx->mod.n, 1);
        } else if (!Bok && Cok) {
            fmpz_CRT_ui(A->coeffs + i, zero              , Bm, C->coeffs[Cind++], fctx->mod.n, 1);
        } else {
            fmpz_zero(A->coeffs + i);
        }

        Bok = 1;
        Cok = 1;
        carry = 1;
        for (j = nvars - 1; j >= 0; j--)
        {
            inds[j] += carry;
            if (inds[j] < A->deg_bounds[j])
            {
                carry = 0;
                Bok = Bok && (inds[j] < B->deg_bounds[j]);
                Cok = Cok && (inds[j] < C->deg_bounds[j]);
            } else
            {
                carry = 1;
                inds[j] = 0;
            }
        }
    }

    fmpz_clear(zero);

    TMP_END;
}


void fmpz_mpolyd_height(fmpz_t max, fmpz_mpolyd_t A)
{
    slong degb_prod, i, j;
    fmpz_t t;

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++)
    {
        degb_prod *= A->deg_bounds[j];
    }

    fmpz_init(t);
    fmpz_zero(max);
    for (i = 0; i < degb_prod; i++)
    {
        fmpz_abs(t, A->coeffs + i);
        if (fmpz_cmp(max, t) < 0)
            fmpz_set(max, t);
    }

    fmpz_clear(t);
}

void fmpz_mpolyd_heights(fmpz_t max, fmpz_t sum, fmpz_mpolyd_t A)
{
    slong degb_prod, i, j;
    fmpz_t t;

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++)
        degb_prod *= A->deg_bounds[j];

    fmpz_init(t);
    fmpz_zero(max);
    fmpz_zero(sum);
    for (i = 0; i < degb_prod; i++)
    {
        fmpz_abs(t, A->coeffs + i);
        fmpz_add(sum, sum, t);
        if (fmpz_cmp(max, t) < 0)
            fmpz_set(max, t);
    }

    fmpz_clear(t);
}

void fmpz_mpolyd_divexact_fmpz(fmpz_mpolyd_t A, fmpz_t c)
{
    slong degb_prod, i, j;

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++)
        degb_prod *= A->deg_bounds[j];

    for (i = 0; i < degb_prod; i++)
    {
        fmpz_divexact(A->coeffs + i, A->coeffs + i, c);
    }
}

void fmpz_mpolyd_mul_scalar(fmpz_mpolyd_t A, fmpz_t c)
{
    slong degb_prod, i, j;

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++)
        degb_prod *= A->deg_bounds[j];

    for (i = 0; i < degb_prod; i++)
    {
        fmpz_mul(A->coeffs + i, A->coeffs + i, c);
    }
}

void fmpz_mpolyd_content(fmpz_t c, const fmpz_mpolyd_t A)
{
    slong degb_prod, j;

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++)
        degb_prod *= A->deg_bounds[j];

    _fmpz_vec_content(c, A->coeffs, degb_prod);
}

void fmpz_mpolyd_to_nmod_mpolyd(nmod_mpolyd_t Ap, fmpz_mpolyd_t A, const nmodf_ctx_t fctx)
{
    slong j;
    slong degb_prod;

    nmod_mpolyd_set_nvars(Ap, A->nvars);

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++)
    {
        Ap->deg_bounds[j] = A->deg_bounds[j];
        degb_prod *= A->deg_bounds[j];
    }

    nmod_mpolyd_fit_length(Ap, degb_prod);
    _fmpz_vec_get_nmod_vec(Ap->coeffs, A->coeffs, degb_prod, fctx->mod);
}

void fmpz_mpolyd_set_nmod_mpolyd(fmpz_mpolyd_t A, nmod_mpolyd_t Ap, const nmodf_ctx_t fctx)
{
    slong j;
    slong degb_prod;

    fmpz_mpolyd_set_nvars(A, Ap->nvars);

    degb_prod = WORD(1);
    for (j = 0; j < Ap->nvars; j++)
    {
        A->deg_bounds[j] = Ap->deg_bounds[j];
        degb_prod *= Ap->deg_bounds[j];
    }

    fmpz_mpolyd_fit_length(A, degb_prod);
    _fmpz_vec_set_nmod_vec(A->coeffs, Ap->coeffs, degb_prod, fctx->mod);
}



void fmpz_mpolyd_set(fmpz_mpolyd_t A, const fmpz_mpolyd_t B)
{
    slong j;
    slong degb_prod;

    fmpz_mpolyd_set_nvars(A, B->nvars);

    degb_prod = WORD(1);
    for (j = 0; j < B->nvars; j++)
    {
        A->deg_bounds[j] = B->deg_bounds[j];
        degb_prod *= B->deg_bounds[j];
    }

    fmpz_mpolyd_fit_length(A, degb_prod);
    _fmpz_vec_set(A->coeffs, B->coeffs, degb_prod);
}


slong fmpz_mpolyd_leadmon(slong * exps, const fmpz_mpolyd_t A)
{

    slong i, j, k;
    slong degb_prod;

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++)
    {
        degb_prod *= A->deg_bounds[j];
    }

    for (i = degb_prod-1; i >= 0; i--)
    {
        if (!fmpz_is_zero(A->coeffs + i))
            break;
    }

    FLINT_ASSERT(i>=0);

    k = i;
    for (j = A->nvars - 1; j >= 0; j--) 
    {
        ulong m = A->deg_bounds[j];
        ulong e = k % m;
        k = k / m;
        exps[j] = e;
    }
    FLINT_ASSERT(k == 0);

    return i;
}

void fmpz_mpolyd_lc(fmpz_t a, const fmpz_mpolyd_t A)
{

    slong i, j, k;
    slong degb_prod;

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++)
    {
        degb_prod *= A->deg_bounds[j];
    }

    for (i = degb_prod-1; i >= 0; i--)
    {
        if (!fmpz_is_zero(A->coeffs + i))
        {
            fmpz_set(a, A->coeffs + i);
            return;
        }
    }
}


int fmpz_mpolyd_gcd_brown(fmpz_mpolyd_t G,
            fmpz_mpolyd_t Abar, fmpz_mpolyd_t Bbar,
                    fmpz_mpolyd_t A, fmpz_mpolyd_t B)
{
    int equal, success;
    mp_limb_t p, gammap;
    slong j, nvars;
    slong lm_idx;
    slong * exp, * texp;
    fmpz_t gamma, m;
    fmpz_t gnm, gns, anm, ans, bnm, bns;
    fmpz_t lA, lB, cA, cB, cG, bound, temp, pp;
    nmod_mpolyd_t Gp, Apbar, Bpbar, Ap, Bp;
    fmpz_mpolyd_t T;
    nmodf_ctx_t fctx;

    nmodf_ctx_init(fctx, 2);
    nvars = A->nvars;

/*printf("\n\n******* starting mpolyd_gcd_brown *********\n");*/

    fmpz_mpolyd_init(T, nvars);

    nmod_mpolyd_init(Gp, nvars);
    nmod_mpolyd_init(Apbar, nvars);
    nmod_mpolyd_init(Bpbar, nvars);
    nmod_mpolyd_init(Ap, nvars);
    nmod_mpolyd_init(Bp, nvars);


/*printf("A: "); fmpz_mpolyd_print_simple(A); printf("\n");
printf("B: "); fmpz_mpolyd_print_simple(B); printf("\n");*/

    fmpz_init(cA);
    fmpz_init(cB);
    fmpz_init(cG);
    fmpz_init(lA);
    fmpz_init(lB);
    fmpz_init(gamma);
    fmpz_init(gnm);
    fmpz_init(gns);
    fmpz_init(anm);
    fmpz_init(ans);
    fmpz_init(bnm);
    fmpz_init(bns);
    fmpz_init(bound);
    fmpz_init(temp);
    fmpz_init_set_si(m, 1);
    fmpz_init(pp);



    fmpz_mpolyd_content(cA, A);
/*printf("cA: "); fmpz_print(cA); printf("\n");*/

    fmpz_mpolyd_content(cB, B);
/*printf("cB: "); fmpz_print(cB); printf("\n");*/

    fmpz_gcd(cG, cA, cB);
    fmpz_mpolyd_divexact_fmpz(A, cA);
    fmpz_mpolyd_divexact_fmpz(B, cB);
/*printf("A: "); fmpz_mpolyd_print_simple(A); printf("\n");
printf("B: "); fmpz_mpolyd_print_simple(B); printf("\n");*/

    fmpz_mpolyd_lc(lA, A);
/*printf("lA: "); fmpz_print(lA); printf("\n");*/

    fmpz_mpolyd_lc(lB, B);
/*printf("lB: "); fmpz_print(lB); printf("\n");*/

    fmpz_gcd(gamma, lA, lB);
/*printf("gamma: "); fmpz_print(gamma); printf("\n");*/

    fmpz_mpolyd_height(bound, A);
/*printf("bound: "); fmpz_print(bound); printf("\n");*/

    fmpz_mpolyd_height(temp, B);
/*printf("temp: "); fmpz_print(temp); printf("\n");*/

    if (fmpz_cmp(bound, temp) < 0)
        fmpz_swap(bound, temp);
    fmpz_mul(bound, bound, gamma);
    fmpz_add(bound, bound, bound);
/*printf("bound: "); fmpz_print(bound); printf("\n");*/


    exp = (slong *) flint_malloc(nvars*sizeof(slong));
    texp = (slong *) flint_malloc(nvars*sizeof(slong));


    fmpz_mpolyd_leadmon(exp, A);
    fmpz_mpolyd_leadmon(texp, B);
    for (j = 0; j < nvars; j++)
        exp[j] = FLINT_MIN(exp[j], texp[j]);

    p = UWORD(1) << (FLINT_BITS - 1);

loop:
    p = n_nextprime(p, 1);

/*sleep(1);*/

    nmodf_ctx_clear(fctx);
    nmodf_ctx_init(fctx, p);

/*printf("m: "); fmpz_print(m); printf("\n");*/

    fmpz_set_ui(pp, p);
    if (fmpz_divisible(lA, pp) || fmpz_divisible(lB, pp))
        goto loop;

/*printf("pp: "); fmpz_print(pp); printf("\n");*/

    fmpz_mpolyd_to_nmod_mpolyd(Ap, A, fctx);
    fmpz_mpolyd_to_nmod_mpolyd(Bp, B, fctx);

/*
printf("Ap: "); nmod_mpolyd_print_simple(Ap); printf("\n");
printf("Bp: "); nmod_mpolyd_print_simple(Bp); printf("\n");
printf("m: "); fmpz_print(m); printf("\n");
*/


    success = nmod_mpolyd_gcd_brown(Gp, Apbar, Bpbar, Ap, Bp, fctx);

/*
printf("success: %d\n", success);
printf("   Gp: "); nmod_mpolyd_print_simple(Gp); printf("\n");
printf("Apbar: "); nmod_mpolyd_print_simple(Apbar); printf("\n");
printf("Bpbar: "); nmod_mpolyd_print_simple(Bpbar); printf("\n");
printf("m: "); fmpz_print(m); printf("\n");
*/

    if (!success)
        goto loop;

    lm_idx = nmod_mpolyd_leadmon(texp, Gp);

/*printf("lm_idx: %d\n", lm_idx);*/

    if (lm_idx <= 0)
    {
        FLINT_ASSERT(leadmon_idx == 0);
        fmpz_mpolyd_set_fmpz(G, cG);
        fmpz_mpolyd_swap(Abar, A);
        fmpz_divexact(temp, cA, cG);
        fmpz_mpolyd_mul_scalar(Abar, temp);
        fmpz_mpolyd_swap(Bbar, B);
        fmpz_divexact(temp, cB, cG);
        fmpz_mpolyd_mul_scalar(Bbar, temp);
        goto done;
    }

    equal = 1;
    for (j = 0; j < nvars; j++)
    {
        if (texp[j] > exp[j])
        {
            goto loop;
        } else if (texp[j] < exp[j])
        {
            equal = 0;
            break;
        }
    }


    gammap = fmpz_fdiv_ui(gamma, p);
    nmod_mpolyd_mul_scalar(Gp, gammap, fctx);
/*
printf("m: "); fmpz_print(m); printf("\n");
printf("   Gp: "); nmod_mpolyd_print_simple(Gp); printf("\n");
printf("Apbar: "); nmod_mpolyd_print_simple(Apbar); printf("\n");
printf("Bpbar: "); nmod_mpolyd_print_simple(Bpbar); printf("\n");
*/


    if (fmpz_is_one(m) || !equal)
    {
        fmpz_mpolyd_set_nmod_mpolyd(G, Gp, fctx);
        fmpz_mpolyd_set_nmod_mpolyd(Abar, Apbar, fctx);
        fmpz_mpolyd_set_nmod_mpolyd(Bbar, Bpbar, fctx);
/*
printf("   G: "); fmpz_mpolyd_print_simple(G); printf("\n");
printf("Abar: "); fmpz_mpolyd_print_simple(Abar); printf("\n");
printf("Bbar: "); fmpz_mpolyd_print_simple(Bbar); printf("\n");
*/
        fmpz_set_ui(m, p);
        for (j = 0; j < nvars; j++)
            exp[j] = texp[j];
        goto loop;
    }

    fmpz_mpolyd_CRT_nmod(T, G, m, Gp, fctx);
    fmpz_mpolyd_swap(T, G);
    fmpz_mpolyd_CRT_nmod(T, Abar, m, Apbar, fctx);
    fmpz_mpolyd_swap(T, Abar);
    fmpz_mpolyd_CRT_nmod(T, Bbar, m, Bpbar, fctx);
    fmpz_mpolyd_swap(T, Bbar);
    fmpz_mul(m, m, pp);

/*printf("after interp\n");
printf("   G: "); fmpz_mpolyd_print_simple(G); printf("\n");
printf("Abar: "); fmpz_mpolyd_print_simple(Abar); printf("\n");
printf("Bbar: "); fmpz_mpolyd_print_simple(Bbar); printf("\n");
printf("m: "); fmpz_print(m); printf("\n");*/

    if (fmpz_cmp(m, bound) <= 0)
        goto loop;

/*printf("more1\n");*/

    fmpz_mpolyd_heights(gnm, gns, G);
    fmpz_mpolyd_heights(anm, ans, Abar);
    fmpz_mpolyd_heights(bnm, bns, Bbar);
    fmpz_mul(ans, ans, gnm);
    fmpz_mul(anm, anm, gns);
    fmpz_mul(bns, bns, gnm);
    fmpz_mul(bnm, bnm, gns);

/*printf("more2\n");*/

    if (fmpz_cmp(ans, anm) > 0)
        fmpz_swap(ans, anm);
    if (fmpz_cmp(bns, bnm) > 0)
        fmpz_swap(bns, bnm);
    fmpz_add(ans, ans, ans);
    fmpz_add(bns, bns, bns);
    if (fmpz_cmp(ans, m) >= 0 || fmpz_cmp(bns, m) >= 0)
        goto loop;

/*printf("no more loop\n");*/

    fmpz_mpolyd_content(temp, G);
    fmpz_mpolyd_divexact_fmpz(G, temp);
    fmpz_mpolyd_lc(temp, G);
    fmpz_mpolyd_divexact_fmpz(Abar, temp);
    fmpz_mpolyd_divexact_fmpz(Bbar, temp);

    fmpz_mpolyd_mul_scalar(G, cG);
    fmpz_divexact(temp, cA, cG);
    fmpz_mpolyd_mul_scalar(Abar, temp);
    fmpz_divexact(temp, cB, cG);
    fmpz_mpolyd_mul_scalar(Bbar, temp);

done:

    flint_free(exp);
    flint_free(texp);

    fmpz_clear(cA);
    fmpz_clear(cB);
    fmpz_clear(cG);
    fmpz_clear(lA);
    fmpz_clear(lB);
    fmpz_clear(gamma);
    fmpz_clear(gnm);
    fmpz_clear(gns);
    fmpz_clear(anm);
    fmpz_clear(ans);
    fmpz_clear(bnm);
    fmpz_clear(bns);
    fmpz_clear(bound);
    fmpz_clear(temp);
    fmpz_clear(m);
    fmpz_clear(pp);

    fmpz_mpolyd_clear(T);

    nmod_mpolyd_clear(Gp);
    nmod_mpolyd_clear(Apbar);
    nmod_mpolyd_clear(Bpbar);
    nmod_mpolyd_clear(Ap);
    nmod_mpolyd_clear(Bp);

    nmodf_ctx_clear(fctx);

    return 1;
}


int fmpz_mpoly_gcd_brown(fmpz_mpoly_t G,
                              const fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    fmpz_mpolyd_t Ad, Bd, Gd, Abar, Bbar;
    fmpz_mpolyd_ctx_t dctx;
    slong nvars = ctx->minfo->nvars;

    success = 1;

    if (fmpz_mpoly_is_zero(A, ctx)) {
        fmpz_mpoly_set(G, B, ctx);
        goto cleanup_stage0;
    }
    if (fmpz_mpoly_is_zero(B, ctx)) {
        fmpz_mpoly_set(G, A, ctx);
        goto cleanup_stage0;
    }
/*printf("here1\n");*/

    fmpz_mpolyd_ctx_init(dctx, nvars);
    success = fmpz_mpolyd_ctx_settle(dctx, A, B, ctx);
    if (!success)
    {
        fmpz_mpoly_zero(G, ctx);
        goto cleanup_stage1;
    }
/*printf("here2\n");*/

    fmpz_mpolyd_init(Ad, nvars);
    fmpz_mpolyd_init(Bd, nvars);
    fmpz_mpolyd_init(Gd, nvars);
    fmpz_mpolyd_init(Abar, nvars);
    fmpz_mpolyd_init(Bbar, nvars);

/*printf("here3\n");*/

    fmpz_mpoly_convert_to_fmpz_mpolyd(Ad, dctx, A, ctx);
    fmpz_mpoly_convert_to_fmpz_mpolyd(Bd, dctx, B, ctx);

/*printf("here4\n");*/

    success = fmpz_mpolyd_gcd_brown(Gd, Abar, Bbar, Ad, Bd);
    if (!success)
    {
        fmpz_mpoly_zero(G, ctx);
    } else
    {
        fmpz_mpoly_convert_from_fmpz_mpolyd(G, ctx, Gd, dctx);
    }

    fmpz_mpolyd_clear(Bbar);
    fmpz_mpolyd_clear(Abar);
    fmpz_mpolyd_clear(Gd);
    fmpz_mpolyd_clear(Bd);
    fmpz_mpolyd_clear(Ad);

cleanup_stage1:

    fmpz_mpolyd_ctx_clear(dctx);

cleanup_stage0:

    if (success && (G->length > 0) && (fmpz_sgn(G->coeffs + 0) < 0))
        fmpz_mpoly_neg(G, G, ctx);

    return success;
}
