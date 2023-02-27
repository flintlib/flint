/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"
#include "fq_nmod_mpoly.h"

int fq_nmod_next(fq_nmod_t alpha, const fq_nmod_ctx_t fqctx)
{
    slong i;
    slong deg = nmod_poly_degree(fqctx->modulus);

    for (i = 0; i < deg; i++)
    {
        ulong c = nmod_poly_get_coeff_ui(alpha, i);
        c += UWORD(1);
        if (c < fqctx->mod.n)
        {
            nmod_poly_set_coeff_ui(alpha, i, c);
            return 1;
        }
        nmod_poly_set_coeff_ui(alpha, i, 0);
    }

    return 0;
}

void fq_nmod_next_not_zero(fq_nmod_t alpha, const fq_nmod_ctx_t fqctx)
{
    slong i;
    slong deg = fqctx->modulus->length - 1;

    for (i = 0; i < deg; i++) {
        ulong c = nmod_poly_get_coeff_ui(alpha, i);
        c += UWORD(1);
        if (c < fqctx->mod.n) {
            nmod_poly_set_coeff_ui(alpha, i, c);
            return;
        }
        nmod_poly_set_coeff_ui(alpha, i, UWORD(0));
    }

    /* we hit zero, so skip to 1 */
    nmod_poly_set_coeff_ui(alpha, 0, UWORD(1));
}


/*
    Assuming that "A" depends only on the main variable,
    convert it to a poly "a".
*/
static void fq_nmod_mpolyu_cvtto_poly(
    fq_nmod_poly_t a,
    fq_nmod_mpolyu_t A,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    fq_nmod_t at;

    FLINT_ASSERT(fq_nmod_mpolyu_is_canonical(A, ctx));

    fq_nmod_init(at, ctx->fqctx);

    fq_nmod_poly_zero(a, ctx->fqctx);
    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT((A->coeffs + i)->length == 1);
        FLINT_ASSERT(mpoly_monomial_is_zero((A->coeffs + i)->exps, 
                      mpoly_words_per_exp((A->coeffs + i)->bits, ctx->minfo)));

        n_fq_get_fq_nmod(at, (A->coeffs + i)->coeffs + 0, ctx->fqctx);
        fq_nmod_poly_set_coeff(a, A->exps[i], at, ctx->fqctx);
    }

    fq_nmod_clear(at, ctx->fqctx);
}

static void fq_nmod_mpolyu_cvtfrom_poly(
    fq_nmod_mpolyu_t A,
    fq_nmod_poly_t a,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong k;
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
    fq_nmod_t c;
    fq_nmod_init(c, ctx->fqctx);

    fq_nmod_mpolyu_zero(A, ctx);
    k = 0;
    for (i = fq_nmod_poly_length(a, ctx->fqctx) - 1; i >= 0; i--)
    {
        fq_nmod_poly_get_coeff(c, a, i, ctx->fqctx);
        if (!fq_nmod_is_zero(c, ctx->fqctx))
        {
            fq_nmod_mpolyu_fit_length(A, k + 1, ctx);
            A->exps[k] = i;
            fq_nmod_mpoly_fit_length_reset_bits(A->coeffs + k, 1, A->bits, ctx);
            n_fq_set_fq_nmod((A->coeffs + k)->coeffs + 0, c, ctx->fqctx);
            (A->coeffs + k)->length = 1;
            mpoly_monomial_zero((A->coeffs + k)->exps + N*0, N);
            k++;
        }
    }
    A->length = k;

    fq_nmod_clear(c, ctx->fqctx);
}



/* store in each coefficient the evaluation of the corresponding monomial */
static void fq_nmod_mpoly_evalsk(
    fq_nmod_mpoly_t A,
    fq_nmod_mpoly_t B,
    slong entries,
    slong * offs,
    ulong * masks,
    fq_nmod_struct * powers,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong i, j;
    slong N;

    FLINT_ASSERT(A->bits == B->bits);
    fq_nmod_mpoly_fit_length(A, B->length, ctx);
    N = mpoly_words_per_exp(B->bits, ctx->minfo);
    for (i = 0; i < B->length; i++)
    {
        _n_fq_one(A->coeffs + d*i, d);

        for (j = 0; j < entries; j++)
        {
            if ((B->exps + N*i)[offs[j]] & masks[j])
            {
                n_fq_mul_fq_nmod(A->coeffs + d*i, A->coeffs + d*i, powers + j, ctx->fqctx);
            }
        }

        mpoly_monomial_set(A->exps + N*i, B->exps + N*i, N);
    }
    A->length = B->length;
}

static void fq_nmod_mpolyu_evalsk(
    fq_nmod_mpolyu_t A,
    fq_nmod_mpolyu_t B,
    slong entries,
    slong * offs,
    ulong * masks,
    fq_nmod_struct * powers,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;

    fq_nmod_mpolyu_fit_length(A, B->length, ctx);
    for (i = 0; i < B->length; i++)
    {
        A->exps[i] = B->exps[i];
        fq_nmod_mpoly_evalsk(A->coeffs + i, B->coeffs + i,
                                            entries, offs, masks, powers, ctx);
    }
    A->length = B->length;
}


/* multiply the coefficients of A pointwise by those of B */
static void fq_nmod_mpolyu_mulsk(
    fq_nmod_mpolyu_t A,
    fq_nmod_mpolyu_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong i, j;

    FLINT_ASSERT(A->length == B->length);
    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT(A->exps[i] == B->exps[i]);

        FLINT_ASSERT((A->coeffs + i)->length == (B->coeffs + i)->length);
        for (j = 0; j < (A->coeffs + i)->length; j++)
        {
            n_fq_mul((A->coeffs + i)->coeffs + d*j,
                            (A->coeffs + i)->coeffs + d*j,
                            (B->coeffs + i)->coeffs + d*j, ctx->fqctx);
        }
    }
}

/*
    return 0 if the leading coeff of A vanishes
    else return 1
*/
static int fq_nmod_mpolyu_evalfromsk(
    fq_nmod_poly_t e,
    fq_nmod_mpolyu_t A,
    fq_nmod_mpolyu_t SK,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong i, j;
    int ret = 0;
    mp_limb_t * pp, * acc;
    fq_nmod_t acct;

    FLINT_ASSERT(A->length == SK->length);

    pp = FLINT_ARRAY_ALLOC(d, mp_limb_t);
    acc = FLINT_ARRAY_ALLOC(d, mp_limb_t);
    fq_nmod_init(acct, ctx->fqctx);

    fq_nmod_poly_zero(e, ctx->fqctx);
    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT((A->coeffs + i)->length == (SK->coeffs + i)->length);

        _n_fq_zero(acc, d);

        for (j = 0; j < A->coeffs[i].length; j++)
        {
            n_fq_mul(pp, A->coeffs[i].coeffs + d*j, SK->coeffs[i].coeffs + d*j, ctx->fqctx);
            n_fq_add(acc, acc, pp, ctx->fqctx);
        }

        n_fq_get_fq_nmod(acct, acc, ctx->fqctx);
        fq_nmod_poly_set_coeff(e, A->exps[i], acct, ctx->fqctx);
        ret |= (i == 0 && !fq_nmod_is_zero(acct, ctx->fqctx));
    }

    flint_free(acc);
    flint_free(pp);
    fq_nmod_clear(acct, ctx->fqctx);

    return ret;
}


void fq_nmod_poly_product_roots(fq_nmod_poly_t P, fq_nmod_struct * r,
                                            slong n, const fq_nmod_ctx_t fqctx)
{
    slong i;
    fq_nmod_poly_t B;
    fq_nmod_t a;
    fq_nmod_init(a, fqctx);
    fq_nmod_poly_init(B, fqctx);
    fq_nmod_poly_one(P, fqctx);
    fq_nmod_poly_gen(B, fqctx);
    for (i = 0; i < n; i++)
    {
        fq_nmod_neg(a, r + i, fqctx);
        fq_nmod_poly_set_coeff(B, 0, a, fqctx);
        fq_nmod_poly_mul(P, P, B, fqctx);
    }
    fq_nmod_clear(a, fqctx);
    fq_nmod_poly_clear(B, fqctx);
}

/*
    solve

    [ a[0]    a[1]    ... a[n-1]   ]   [ x[0]   ]    [ b[0]   ]
    [ a[0]^2  a[1]^2  ... a[n-1]^2 ]   [ x[1]   ]    [ b[1]   ]
    [  ...                 ...     ] . [ ...    ] =  [ ...    ]
    [ a[0]^n  a[1]^n  ... a[n-1]^n ]   [ x[n-1] ]    [ b[n-1] ]

    for x
*/
int fq_nmod_vandsolve(mp_limb_t * X, mp_limb_t * A, fq_nmod_struct * b,
                                            slong n, const fq_nmod_ctx_t fqctx)
{
    slong d = fq_nmod_ctx_degree(fqctx);
    int success = 0;
    slong i, j;
    fq_nmod_t t, s;
    fq_nmod_t Dinv;
    fq_nmod_poly_t Q, P, R, u;
    fq_nmod_struct * x, * a;

    x = FLINT_ARRAY_ALLOC(n, fq_nmod_struct);
    a = FLINT_ARRAY_ALLOC(n, fq_nmod_struct);

    fq_nmod_init(t, fqctx);
    fq_nmod_init(s, fqctx);
    fq_nmod_init(Dinv, fqctx);

    for (i = 0; i < n; i++)
    {
        fq_nmod_init(a + i, fqctx);
        fq_nmod_init(x + i, fqctx);
        n_fq_get_fq_nmod(a + i, A + d*i, fqctx);
        fq_nmod_zero(x + i, fqctx);
    }

    fq_nmod_poly_init(Q, fqctx);
    fq_nmod_poly_init(P, fqctx);
    fq_nmod_poly_init(R, fqctx);
    fq_nmod_poly_init(u, fqctx);
    fq_nmod_poly_gen(u, fqctx);
    fq_nmod_poly_product_roots(P, a, n, fqctx);
    for (i = 0; i < n; i++)
    {
        if (fq_nmod_is_zero(a + i, fqctx))
            goto cleanup;

        fq_nmod_neg(t, a + i, fqctx);
        fq_nmod_poly_set_coeff(u, 0, t, fqctx);
        fq_nmod_poly_divrem(Q, R, P, u, fqctx);
        fq_nmod_poly_evaluate_fq_nmod(t, Q, a + i, fqctx);
        fq_nmod_mul(t, a + i, t, fqctx);
        if (fq_nmod_is_zero(t, fqctx))
            goto cleanup;

        fq_nmod_inv(Dinv, t, fqctx);
        for (j = 0; j < n; j++)
        {
            fq_nmod_mul(t, b + j, Dinv, fqctx);
            fq_nmod_poly_get_coeff(s, Q, j, fqctx);
            fq_nmod_mul(t, t, s, fqctx);
            fq_nmod_add(x + i, x + i, t, fqctx);
        }
    }
    success = 1;

cleanup:
    fq_nmod_poly_clear(Q, fqctx);
    fq_nmod_poly_clear(P, fqctx);
    fq_nmod_poly_clear(R, fqctx);
    fq_nmod_poly_clear(u, fqctx);

    fq_nmod_clear(t, fqctx);
    fq_nmod_clear(s, fqctx);
    fq_nmod_clear(Dinv, fqctx);

    for (i = 0; i < n; i++)
    {
        n_fq_set_fq_nmod(X + d*i, x + i, fqctx);
        fq_nmod_clear(a + i, fqctx);
        fq_nmod_clear(x + i, fqctx);
    }

    flint_free(a);
    flint_free(x);

    return success;
}

/*
    Try to set G to the gcd of A and B given the form f of G.
    return codes as enumerated in nmod_mpoly.h:

    nmod_mpoly_gcds_success,
    nmod_mpoly_gcds_form_wrong,
    nmod_mpoly_gcds_no_solution,
    nmod_mpoly_gcds_scales_not_found,
    nmod_mpoly_gcds_eval_point_not_found,
    nmod_mpoly_gcds_eval_gcd_deg_too_high
*/
nmod_gcds_ret_t fq_nmod_mpolyu_gcds_zippel(
    fq_nmod_mpolyu_t G,
    fq_nmod_mpolyu_t A,
    fq_nmod_mpolyu_t B,
    fq_nmod_mpolyu_t f,
    slong var,
    const fq_nmod_mpoly_ctx_t ctx,
    flint_rand_t randstate,
    slong * degbound)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    int eval_points_tried;
    nmod_gcds_ret_t success;
    fq_nmod_mpolyu_t Aevalsk1, Bevalsk1, fevalsk1, Aevalski, Bevalski, fevalski;
    fq_nmod_poly_t Aeval, Beval, Geval;
    fq_nmod_struct * alpha, * b;
    fq_nmod_mat_struct * M, * ML;
    fq_nmod_mat_t MF, Msol;
    int lc_ok;
    int underdeterminedcount = 0;
    int exceededcount = 0;
    int * ML_is_initialized;
    slong i, j, k, s, S, nullity;
    slong * tlen;
    slong l;
    fq_nmod_struct * W;
    slong entries;
    slong * offs;
    ulong * masks;
    fq_nmod_struct * powers;
    fq_nmod_t ck;
    TMP_INIT;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(f->length > 0);

    FLINT_ASSERT(A->bits == B->bits);
    FLINT_ASSERT(A->bits == G->bits);
    FLINT_ASSERT(A->bits == f->bits);
    FLINT_ASSERT(var >= 0);

    FLINT_ASSERT(*degbound == f->exps[0]);

    FLINT_ASSERT(fq_nmod_mpolyu_is_canonical(A, ctx));
    FLINT_ASSERT(fq_nmod_mpolyu_is_canonical(B, ctx));

    FLINT_ASSERT(var > 0);

    if (f->length == 1)
    {
        if ((f->coeffs + 0)->length > 1)
        {
            /* impossible to find scale factors in this case */
            return nmod_gcds_scales_not_found;
        }
        else
        {
            /* otherwise set the coeff of the monomial to one */
            nmod_gcds_ret_t ret;
            FLINT_ASSERT((f->coeffs + 0)->length == 1);
            fq_nmod_mpolyu_set(G, f, ctx);
            _n_fq_one((G->coeffs + 0)->coeffs + d*0, d);
            fq_nmod_mpolyu_init(Aevalsk1, f->bits, ctx);
            ret = nmod_gcds_form_wrong;
            if (   fq_nmod_mpolyuu_divides(Aevalsk1, A, G, 1, ctx)
                && fq_nmod_mpolyuu_divides(Aevalsk1, B, G, 1, ctx))
            {
                ret = nmod_gcds_success;
            }
            fq_nmod_mpolyu_clear(Aevalsk1, ctx);
            return ret;
        }
    }

    TMP_START;

    fq_nmod_init(ck, ctx->fqctx);

    fq_nmod_mpolyu_init(Aevalsk1, f->bits, ctx);
    fq_nmod_mpolyu_init(Bevalsk1, f->bits, ctx);
    fq_nmod_mpolyu_init(fevalsk1, f->bits, ctx);
    fq_nmod_mpolyu_init(Aevalski, f->bits, ctx);
    fq_nmod_mpolyu_init(Bevalski, f->bits, ctx);
    fq_nmod_mpolyu_init(fevalski, f->bits, ctx);
    fq_nmod_poly_init(Aeval, ctx->fqctx);
    fq_nmod_poly_init(Beval, ctx->fqctx);
    fq_nmod_poly_init(Geval, ctx->fqctx);

    tlen = (slong *) TMP_ALLOC(f->length*sizeof(slong));
    for (i = 0; i < f->length; i++)
    {
        tlen[i] = i;
    }

    /*
        make tlen sort the coeffs so that
        (f->coeffs + tlen[j-1])->length <= (f->coeffs + tlen[j-0])->length
        for all j
    */
    for (i = 1; i<f->length; i++)
    {
        for (j=i; j > 0 && (f->coeffs + tlen[j-1])->length 
                         > (f->coeffs + tlen[j-0])->length; j--)
        {
            slong temp = tlen[j-1];
            tlen[j-1] = tlen[j-0];
            tlen[j-0] = temp;
        }
    }

    /* l is the number of images we will try to construct */
    l = f->length - 3;
    for (i = 0; i < f->length; i++)
    {
        l += (f->coeffs + i)->length;
    }
    l = l / (f->length - 1);
    l = FLINT_MAX(l, (f->coeffs + tlen[f->length - 1])->length);
    /* one extra test image */
    l += 1;

    alpha = (fq_nmod_struct *) TMP_ALLOC(var*sizeof(fq_nmod_struct));
    ML = (fq_nmod_mat_struct *) TMP_ALLOC(f->length*sizeof(fq_nmod_mat_struct));
    b = (fq_nmod_struct *) TMP_ALLOC((f->coeffs + tlen[f->length - 1])->length
                                                      *sizeof(fq_nmod_struct));

    for (i = 0; i < var; i++)
        fq_nmod_init(alpha + i, ctx->fqctx);
    for (i = 0; i < (f->coeffs + tlen[f->length - 1])->length; i++)
        fq_nmod_init(b + i, ctx->fqctx);

    fq_nmod_mat_init(MF, 0, l, ctx->fqctx);

    M = (fq_nmod_mat_struct *) TMP_ALLOC(f->length*sizeof(fq_nmod_mat_struct));
    ML_is_initialized = (int *) TMP_ALLOC(f->length*sizeof(int));
    for (i = 0; i < f->length; i++)
    {
        fq_nmod_mat_init(M + i, l, (f->coeffs + i)->length, ctx->fqctx);
        ML_is_initialized[i] = 0;
    }

    W = (fq_nmod_struct *) flint_malloc(l*f->length*sizeof(fq_nmod_struct));
    for (i = 0; i < l*f->length; i++)
        fq_nmod_init(W + i, ctx->fqctx);
    

    fq_nmod_mat_init(Msol, l, 1, ctx->fqctx);

    /* compute how many masks are needed */
    entries = f->bits * var;
    offs = (slong *) TMP_ALLOC(entries*sizeof(slong));
    masks = (ulong *) TMP_ALLOC(entries*sizeof(slong));
    powers = (fq_nmod_struct *) TMP_ALLOC(entries*sizeof(fq_nmod_struct));
    for (i = 0; i < entries; i++)
        fq_nmod_init(powers + i, ctx->fqctx);


    /***** evaluation loop head *******/
    eval_points_tried = 0;
pick_evaluation_point:

    if (++eval_points_tried > 10)
    {
        success = nmod_gcds_eval_point_not_found;
        goto finished;
    }

    /* avoid 0 for the evaluation points */
    for (i = 0; i < var; i++)
        fq_nmod_randtest_not_zero(alpha + i, randstate, ctx->fqctx);

    /* store bit masks for each power of two of the non-main variables */
    for (i = 0; i < var; i++)
    {
        slong shift, off;
        mpoly_gen_offset_shift_sp(&off, &shift, i, f->bits, ctx->minfo);
        for (j = 0; j < f->bits; j++)
        {
            offs[f->bits*i + j] = off;
            masks[f->bits*i + j] = UWORD(1) << (j + shift);
            if (j == 0)
                fq_nmod_set(powers + f->bits*i + j, alpha + i, ctx->fqctx);
            else
                fq_nmod_mul(powers + f->bits*i + j, powers + f->bits*i + j - 1,
                                                    powers + f->bits*i + j - 1, ctx->fqctx);
        }
    }

    fq_nmod_mpolyu_evalsk(Aevalsk1, A, entries, offs, masks, powers, ctx);
    fq_nmod_mpolyu_evalsk(Bevalsk1, B, entries, offs, masks, powers, ctx);
    fq_nmod_mpolyu_evalsk(fevalsk1, f, entries, offs, masks, powers, ctx);

    for (i = 0; i < l*f->length; i++)
    {
        fq_nmod_zero(W + i, ctx->fqctx);
    }


    for (i = 0; i < l; i++)
    {
        if (i == 0)
        {
            fq_nmod_mpolyu_set(Aevalski, Aevalsk1, ctx);
            fq_nmod_mpolyu_set(Bevalski, Bevalsk1, ctx);
            fq_nmod_mpolyu_set(fevalski, fevalsk1, ctx);
        } else
        {
            fq_nmod_mpolyu_mulsk(Aevalski, Aevalsk1, ctx);
            fq_nmod_mpolyu_mulsk(Bevalski, Bevalsk1, ctx);
            fq_nmod_mpolyu_mulsk(fevalski, fevalsk1, ctx);
        }

        for (j = 0; j < f->length; j++)
        {
            for (k = 0; k < (f->coeffs + j)->length; k++)
            {
                n_fq_get_fq_nmod((M + j)->rows[i] + k,
                               (fevalski->coeffs + j)->coeffs + d*k, ctx->fqctx);
            }
        }

        lc_ok = 1;
        lc_ok = lc_ok && fq_nmod_mpolyu_evalfromsk(Aeval, A, Aevalski, ctx);
        lc_ok = lc_ok && fq_nmod_mpolyu_evalfromsk(Beval, B, Bevalski, ctx);
        if (!lc_ok)
        {
            /* lc of A or B vanished */
            goto pick_evaluation_point;
        }

        fq_nmod_poly_gcd(Geval, Aeval, Beval, ctx->fqctx);

        if (f->exps[0] < fq_nmod_poly_degree(Geval, ctx->fqctx))
        {
            ++exceededcount;
            if (exceededcount < 2)
                goto pick_evaluation_point;

            success = nmod_gcds_eval_gcd_deg_too_high;
            goto finished;
        }

        if (f->exps[0] > fq_nmod_poly_degree(Geval, ctx->fqctx))
        {
            success = nmod_gcds_form_main_degree_too_high;
            *degbound = fq_nmod_poly_degree(Geval, ctx->fqctx);
            goto finished;
        }

        k = fq_nmod_poly_length(Geval, ctx->fqctx);
        j = WORD(0);
        while ((--k) >= 0)
        {
            fq_nmod_poly_get_coeff(ck, Geval, k, ctx->fqctx);
            if (!fq_nmod_is_zero(ck, ctx->fqctx))
            {
                while (j < f->length && f->exps[j] > k)
                {
                    j++;
                }
                if (j >= f->length || f->exps[j] != k)
                {
                    success = nmod_gcds_form_wrong;
                    goto finished;
                }
                fq_nmod_set(W + l*j + i, ck, ctx->fqctx);
            }
        }
    }

    nullity = -1;
    fq_nmod_mat_clear(MF, ctx->fqctx);
    fq_nmod_mat_init(MF, 0, l, ctx->fqctx);

    for (S = 0; S < f->length; S++)
    {
        s = tlen[S];

        if (!ML_is_initialized[s])
        {
            fq_nmod_mat_init(ML + s, l, (f->coeffs + s)->length + l, ctx->fqctx);
            ML_is_initialized[s] = 1;
            for (i = 0; i < l; i++)
            {
                for (j = 0; j < (f->coeffs + s)->length; j++)
                {
                    fq_nmod_set((ML + s)->rows[i] + j,
                                 (M + s)->rows[i] + j, ctx->fqctx);
                }
                fq_nmod_set((ML + s)->rows[i] + (f->coeffs + s)->length + i,
                                                      W + l*s + i, ctx->fqctx);
            }
        } else {
            for (i = 0; i < l; i++)
            {
                for (j = 0; j < (f->coeffs + s)->length; j++)
                {
                    fq_nmod_set((ML + s)->rows[i] + j,
                                             (M + s)->rows[i] + j, ctx->fqctx);
                }
                for (j = 0; j < l; j++) {
                    if (j == i)
                    {
                        fq_nmod_set((ML + s)->rows[i]
                                     + (f->coeffs + s)->length + j,
                                                      W + l*s + i, ctx->fqctx);
                    } else {
                        fq_nmod_zero((ML + s)->rows[i]
                                     + (f->coeffs + s)->length + j,
                                                                   ctx->fqctx);
                    }
                }
            }

        }
        fq_nmod_mat_rref(ML + s, ctx->fqctx);

        for (i = 0; i < (f->coeffs + s)->length; i++)
        {
            if (!fq_nmod_is_one((ML + s)->rows[i] + i, ctx->fqctx))
            {
                /* evaluation points produced a singular vandermonde matrix */
                goto pick_evaluation_point;
            }
        }

        {
            /* appends rows to MF */
            fq_nmod_mat_t MFtemp;
            fq_nmod_mat_t Mwindow;

            fq_nmod_mat_window_init(Mwindow, ML + s,
                    (f->coeffs + s)->length, (f->coeffs + s)->length,
                     l, (f->coeffs + s)->length + l, ctx->fqctx);
            fq_nmod_mat_init(MFtemp,
                             fq_nmod_mat_nrows(MF, ctx->fqctx)
                                     + l - (f->coeffs + s)->length,
                             l, ctx->fqctx);
            fq_nmod_mat_concat_vertical(MFtemp, MF, Mwindow, ctx->fqctx);
            fq_nmod_mat_swap(MFtemp, MF, ctx->fqctx);
            fq_nmod_mat_clear(MFtemp, ctx->fqctx);
            fq_nmod_mat_window_clear(Mwindow, ctx->fqctx);
        }

        nullity = l - fq_nmod_mat_rref(MF, ctx->fqctx);

        if (nullity == 0)
        {
            /* There is no solution for scale factors. Form f must be wrong */
            success = nmod_gcds_form_wrong;
            goto finished;
        }
        if (nullity == 1)
        {
            /*
                There is one solution for scale factors based on equations
                considered thus far. Accept this as a solution and perform
                checks of the remaining equations at the end.
            */
            break;
        }
    }

    if (nullity != 1)
    {
        ++underdeterminedcount;
        if (underdeterminedcount < 2)
            goto pick_evaluation_point;

        success = nmod_gcds_scales_not_found;
        goto finished;
    }

    nullity = fq_nmod_mat_nullspace(Msol, MF, ctx->fqctx);
    FLINT_ASSERT(nullity == 1);

    fq_nmod_mpolyu_setform(G, f, ctx);

    for (i = 0; i < f->length; i++)
    {
        for (j = 0; j < (f->coeffs + i)->length; j++)
        {
            FLINT_ASSERT((f->coeffs + i)->length <= l);
            FLINT_ASSERT(j < (f->coeffs + tlen[f->length - 1])->length);
            fq_nmod_mul(b + j, W + l*i + j, Msol->rows[j] + 0, ctx->fqctx);
        }
        success = fq_nmod_vandsolve((G->coeffs + i)->coeffs,
                                 (fevalsk1->coeffs + i)->coeffs, b, 
                                    (f->coeffs + i)->length, ctx->fqctx);
        if (!success)
        {
            /* evaluation points produced a singular vandermonde matrix */
            goto pick_evaluation_point;
        }
    }

    /* check solution */
    for (s = 0; s < f->length; s++)
    {
        fq_nmod_t acc, pp, u;
        fq_nmod_init(acc, ctx->fqctx);
        fq_nmod_init(pp, ctx->fqctx);
        fq_nmod_init(u, ctx->fqctx);

        for (i = 0; i < l; i++)
        {
            fq_nmod_zero(acc, ctx->fqctx);
            for (j = 0; j < (f->coeffs + s)->length; j++)
            {
                n_fq_get_fq_nmod(u, (G->coeffs + s)->coeffs + d*j, ctx->fqctx);
                fq_nmod_mul(pp, (M + s)->rows[i] + j, u, ctx->fqctx);
                fq_nmod_add(acc, acc, pp, ctx->fqctx);
            }

            fq_nmod_mul(u, W + l*s + i, Msol->rows[i] + 0, ctx->fqctx);
            if (!fq_nmod_equal(acc, u, ctx->fqctx))
            {
                fq_nmod_clear(acc, ctx->fqctx);
                fq_nmod_clear(pp, ctx->fqctx);
                fq_nmod_clear(u, ctx->fqctx);
                success = nmod_gcds_no_solution;
                goto finished;
            }
        }

        fq_nmod_clear(acc, ctx->fqctx);
        fq_nmod_clear(pp, ctx->fqctx);
        fq_nmod_clear(u, ctx->fqctx);
    }

    success = nmod_gcds_success;

finished:

    for (i = 0; i < entries; i++)
        fq_nmod_clear(powers + i, ctx->fqctx);

    for (i = 0; i < var; i++)
        fq_nmod_clear(alpha + i, ctx->fqctx);
    for (i = 0; i < (f->coeffs + tlen[f->length - 1])->length; i++)
        fq_nmod_clear(b + i, ctx->fqctx);
    for (i = 0; i < l*f->length; i++)
        fq_nmod_clear(W + i, ctx->fqctx);

    flint_free(W);
    fq_nmod_mat_clear(MF, ctx->fqctx);
    fq_nmod_mat_clear(Msol, ctx->fqctx);
    for (i = 0; i < f->length; i++)
    {
        fq_nmod_mat_clear(M + i, ctx->fqctx);
        if (ML_is_initialized[i])
        {
            fq_nmod_mat_clear(ML + i, ctx->fqctx);
        }
    }

    fq_nmod_clear(ck, ctx->fqctx);

    fq_nmod_mpolyu_clear(Aevalsk1, ctx);
    fq_nmod_mpolyu_clear(Bevalsk1, ctx);
    fq_nmod_mpolyu_clear(fevalsk1, ctx);
    fq_nmod_mpolyu_clear(Aevalski, ctx);
    fq_nmod_mpolyu_clear(Bevalski, ctx);
    fq_nmod_mpolyu_clear(fevalski, ctx);
    fq_nmod_poly_clear(Aeval, ctx->fqctx);
    fq_nmod_poly_clear(Beval, ctx->fqctx);
    fq_nmod_poly_clear(Geval, ctx->fqctx);

    TMP_END;
    return success;
}


/* setform copies the exponents and zeros the coefficients */
void fq_nmod_mpoly_setform_mpolyn(
    fq_nmod_mpoly_t A,
    fq_nmod_mpolyn_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong i;
    slong N;

    FLINT_ASSERT(A->bits == B->bits);

    fq_nmod_mpoly_fit_length(A, B->length, ctx);
    N = mpoly_words_per_exp(B->bits, ctx->minfo);
    for (i = 0; i < B->length; i++)
    {
        _n_fq_zero(A->coeffs + d*i, d);
        mpoly_monomial_set(A->exps + N*i, B->exps + N*i, N);
    }
    A->length = B->length;
}

void fq_nmod_mpolyu_setform_mpolyun(fq_nmod_mpolyu_t A, fq_nmod_mpolyun_t B,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;

    fq_nmod_mpolyu_fit_length(A, B->length, ctx);
    for (i = 0; i < B->length; i++)
    {
        FLINT_ASSERT((B->coeffs + i)->bits == B->bits);
        fq_nmod_mpoly_setform_mpolyn(A->coeffs + i, B->coeffs + i, ctx);
        A->exps[i] = B->exps[i];
    }
    A->length = B->length;
}


int fq_nmod_mpolyu_gcdp_zippel_univar(
    fq_nmod_mpolyu_t G,
    fq_nmod_mpolyu_t Abar,
    fq_nmod_mpolyu_t Bbar,
    fq_nmod_mpolyu_t A,
    fq_nmod_mpolyu_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    fq_nmod_poly_t a, b, g, t, r;
    FLINT_ASSERT(A->bits == B->bits);
    fq_nmod_poly_init(a, ctx->fqctx);
    fq_nmod_poly_init(b, ctx->fqctx);
    fq_nmod_poly_init(g, ctx->fqctx);
    fq_nmod_poly_init(t, ctx->fqctx);
    fq_nmod_poly_init(r, ctx->fqctx);
    fq_nmod_mpolyu_cvtto_poly(a, A, ctx);
    fq_nmod_mpolyu_cvtto_poly(b, B, ctx);
    fq_nmod_poly_gcd(g, a, b, ctx->fqctx);
    fq_nmod_mpolyu_cvtfrom_poly(G, g, ctx);
    fq_nmod_poly_divrem(t, r, a, g, ctx->fqctx);
    fq_nmod_mpolyu_cvtfrom_poly(Abar, t, ctx);
    fq_nmod_poly_divrem(t, r, b, g, ctx->fqctx);
    fq_nmod_mpolyu_cvtfrom_poly(Bbar, t, ctx);
    fq_nmod_poly_clear(a, ctx->fqctx);
    fq_nmod_poly_clear(b, ctx->fqctx);
    fq_nmod_poly_clear(g, ctx->fqctx);
    fq_nmod_poly_clear(t, ctx->fqctx);
    fq_nmod_poly_clear(r, ctx->fqctx);
    return 1;
}

int fq_nmod_mpolyu_gcdp_zippel_univar_no_cofactors(
    fq_nmod_mpolyu_t G,
    fq_nmod_mpolyu_t A,
    fq_nmod_mpolyu_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    fq_nmod_poly_t a, b, g;
    FLINT_ASSERT(A->bits == B->bits);

    FLINT_ASSERT(fq_nmod_mpolyu_is_canonical(A,ctx));
    FLINT_ASSERT(fq_nmod_mpolyu_is_canonical(B,ctx));

    fq_nmod_poly_init(a, ctx->fqctx);
    fq_nmod_poly_init(b, ctx->fqctx);
    fq_nmod_poly_init(g, ctx->fqctx);
    fq_nmod_mpolyu_cvtto_poly(a, A, ctx);
    fq_nmod_mpolyu_cvtto_poly(b, B, ctx);

    fq_nmod_poly_gcd(g, a, b, ctx->fqctx);

    fq_nmod_mpolyu_cvtfrom_poly(G, g, ctx);
    fq_nmod_poly_clear(a, ctx->fqctx);
    fq_nmod_poly_clear(b, ctx->fqctx);
    fq_nmod_poly_clear(g, ctx->fqctx);
    return 1;
}


int fq_nmod_mpolyu_gcdp_zippel_bivar(
    fq_nmod_mpolyu_t G,
    fq_nmod_mpolyu_t Abar,
    fq_nmod_mpolyu_t Bbar,
    fq_nmod_mpolyu_t A,
    fq_nmod_mpolyu_t B,
    const fq_nmod_mpoly_ctx_t ctx,
    flint_rand_t randstate)
{
    slong var = 0;
    slong Alastdeg;
    slong Blastdeg;
    ulong Ashift, Bshift, Gshift;
    slong lastdeg;
    slong bound;
    int success = 0, changed, have_enough;
    fq_nmod_poly_t a, b, c, g, modulus, tempmod, tmp1, tmp2;
    fq_nmod_mpolyu_t Aeval, Beval, Geval;
    fq_nmod_mpolyun_t An, Bn, H, Ht;
    fq_nmod_t geval, temp, alpha, alphastart;
    fmpz_t minusone;

    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(var >= -WORD(1));

    FLINT_ASSERT(G->bits == A->bits);
    FLINT_ASSERT(G->bits == B->bits);

    fmpz_init(minusone);
    fmpz_set_si(minusone, -WORD(1));

    fq_nmod_mpolyun_init(An, A->bits, ctx);
    fq_nmod_mpolyun_init(Bn, A->bits, ctx);
    fq_nmod_mpolyu_cvtto_mpolyun(An, A, var, ctx);
    fq_nmod_mpolyu_cvtto_mpolyun(Bn, B, var, ctx);

    FLINT_ASSERT(An->length > 0);
    FLINT_ASSERT(Bn->length > 0);

    Ashift = A->exps[A->length - 1];
    Bshift = B->exps[B->length - 1];
    Gshift = FLINT_MIN(Ashift, Bshift);
    fq_nmod_mpolyun_shift_right(An, Ashift);
    fq_nmod_mpolyun_shift_right(Bn, Bshift);

    fq_nmod_poly_init(a, ctx->fqctx);
    fq_nmod_poly_init(b, ctx->fqctx);
    fq_nmod_poly_init(c, ctx->fqctx);
    fq_nmod_poly_init(g, ctx->fqctx);
    fq_nmod_poly_init(tmp1, ctx->fqctx);
    fq_nmod_poly_init(tmp2, ctx->fqctx);

    /* if the gcd has content wrt last variable, we are going to fail */
    fq_nmod_mpolyun_content_poly(a, An, ctx);
    fq_nmod_mpolyun_content_poly(b, Bn, ctx);
    fq_nmod_mpolyun_divexact_poly(An, An, a, ctx);
    fq_nmod_mpolyun_divexact_poly(Bn, Bn, b, ctx);
    fq_nmod_poly_gcd(c, a, b, ctx->fqctx);

    n_fq_poly_get_fq_nmod_poly(tmp1, fq_nmod_mpolyun_leadcoeff_poly(An, ctx), ctx->fqctx);
    n_fq_poly_get_fq_nmod_poly(tmp2, fq_nmod_mpolyun_leadcoeff_poly(Bn, ctx), ctx->fqctx);
    fq_nmod_poly_gcd(g, tmp1, tmp2, ctx->fqctx);

    Alastdeg = fq_nmod_mpolyun_lastdeg(An, ctx);
    Blastdeg = fq_nmod_mpolyun_lastdeg(Bn, ctx);

    /* bound of the number of images required */
    bound = 1 + FLINT_MIN(Alastdeg, Blastdeg)
              + fq_nmod_poly_degree(g, ctx->fqctx);

    fq_nmod_poly_init(modulus, ctx->fqctx);
    fq_nmod_poly_init(tempmod, ctx->fqctx);
    fq_nmod_poly_set_coeff_fmpz(tempmod, 1, minusone, ctx->fqctx);
    fq_nmod_mpolyu_init(Aeval, A->bits, ctx);
    fq_nmod_mpolyu_init(Beval, A->bits, ctx);
    fq_nmod_mpolyu_init(Geval, A->bits, ctx);
    fq_nmod_mpolyun_init(H, A->bits, ctx);
    fq_nmod_mpolyun_init(Ht, A->bits, ctx);

    fq_nmod_init(temp, ctx->fqctx);
    fq_nmod_init(geval, ctx->fqctx);
    fq_nmod_init(alpha, ctx->fqctx);
    fq_nmod_init(alphastart, ctx->fqctx);

    /* fail if the gcd has content wrt last variable */
    if (fq_nmod_poly_degree(c, ctx->fqctx) > 0)
    {
        success = 0;
        goto finished;
    }

    fq_nmod_poly_one(modulus, ctx->fqctx);
    fq_nmod_mpolyun_zero(H, ctx);

    fq_nmod_randtest_not_zero(alphastart, randstate, ctx->fqctx);
    fq_nmod_set(alpha, alphastart, ctx->fqctx);

    while (1)
    {
        /* get new evaluation point */
        fq_nmod_next_not_zero(alpha, ctx->fqctx);

        if (fq_nmod_equal(alpha, alphastart, ctx->fqctx))
        {
            success = 0;
            goto finished;
        }

        /* make sure evaluation point does not kill both lc(A) and lc(B) */
        fq_nmod_poly_evaluate_fq_nmod(geval, g, alpha, ctx->fqctx);

        if (fq_nmod_is_zero(geval, ctx->fqctx))
            goto outer_continue;

        /* make sure evaluation point does not kill either A or B */
        fq_nmod_mpolyun_interp_reduce_sm_mpolyu(Aeval, An, alpha, ctx);
        fq_nmod_mpolyun_interp_reduce_sm_mpolyu(Beval, Bn, alpha, ctx);

        FLINT_ASSERT(fq_nmod_mpolyu_is_canonical(Aeval, ctx));
        FLINT_ASSERT(fq_nmod_mpolyu_is_canonical(Beval, ctx));

        if (Aeval->length == 0 || Beval->length == 0)
            goto outer_continue;

        fq_nmod_mpolyu_gcdp_zippel_univar_no_cofactors(Geval, Aeval, Beval, ctx);

        if (fq_nmod_mpolyu_is_one(Geval, ctx))
        {
            fq_nmod_mpolyu_one(G, ctx);
            fq_nmod_mpolyu_swap(Abar, A, ctx);
            fq_nmod_mpolyu_swap(Bbar, B, ctx);
            fq_nmod_mpolyu_shift_left(G, Gshift);
            fq_nmod_mpolyu_shift_left(Abar, Ashift - Gshift);
            fq_nmod_mpolyu_shift_left(Bbar, Bshift - Gshift);
            success = 1;
            goto finished;
        }

        FLINT_ASSERT(Geval->length > 0);

        if (fq_nmod_poly_degree(modulus, ctx->fqctx) > 0)
        {
            if (Geval->exps[0] > H->exps[0])
            {
                goto outer_continue;
            }
            else if (Geval->exps[0] < H->exps[0])
            {
                fq_nmod_poly_one(modulus, ctx->fqctx);                
            }
        }

        /* update interpolant H */
        n_fq_get_fq_nmod(temp, fq_nmod_mpolyu_leadcoeff(Geval, ctx), ctx->fqctx);
        fq_nmod_inv(temp, temp, ctx->fqctx);
        fq_nmod_mul(temp, geval, temp, ctx->fqctx);
        fq_nmod_mpolyu_scalar_mul_fq_nmod(Geval, temp, ctx);

        if (fq_nmod_poly_degree(modulus, ctx->fqctx) > 0)
        {
            fq_nmod_poly_evaluate_fq_nmod(temp, modulus, alpha, ctx->fqctx);
            fq_nmod_inv(temp, temp, ctx->fqctx);
            fq_nmod_poly_scalar_mul_fq_nmod(modulus, modulus, temp, ctx->fqctx);
            changed = fq_nmod_mpolyun_interp_crt_sm_mpolyu(&lastdeg, H, Ht, Geval,
                                                          modulus, alpha, ctx);
            fq_nmod_poly_set_coeff(tempmod, 0, alpha, ctx->fqctx);
            fq_nmod_poly_mul(modulus, modulus, tempmod, ctx->fqctx);

            have_enough = fq_nmod_poly_degree(modulus, ctx->fqctx) >= bound;

            if (changed && !have_enough)
            {
                goto outer_continue;
            }

            if (!changed || have_enough)
            {
                fq_nmod_mpolyun_content_poly(a, H, ctx);
                fq_nmod_mpolyun_divexact_poly(Ht, H, a, ctx);
                fq_nmod_mpolyun_shift_left(Ht, Gshift);
                fq_nmod_mpolyu_cvtfrom_mpolyun(G, Ht, var, ctx);

                if (    fq_nmod_mpolyuu_divides(Abar, A, G, 1, ctx)
                     && fq_nmod_mpolyuu_divides(Bbar, B, G, 1, ctx))
                {
                    success = 1;
                    goto finished;
                }
            }

            if (have_enough)
            {
                fq_nmod_poly_one(modulus, ctx->fqctx);
                goto outer_continue;
            }
        }
        else
        {
            fq_nmod_mpolyun_interp_lift_sm_mpolyu(H, Geval, ctx);
            lastdeg = WORD(0);
            fq_nmod_poly_set_coeff(tempmod, 0, alpha, ctx->fqctx);
            fq_nmod_poly_mul(modulus, modulus, tempmod, ctx->fqctx);
        }

outer_continue:;
    }

    success = 0;

finished:

    fmpz_clear(minusone);
    fq_nmod_clear(temp, ctx->fqctx);
    fq_nmod_clear(geval, ctx->fqctx);
    fq_nmod_clear(alpha, ctx->fqctx);
    fq_nmod_clear(alphastart, ctx->fqctx);
    fq_nmod_poly_clear(a, ctx->fqctx);
    fq_nmod_poly_clear(b, ctx->fqctx);
    fq_nmod_poly_clear(c, ctx->fqctx);
    fq_nmod_poly_clear(g, ctx->fqctx);
    fq_nmod_poly_clear(tmp1, ctx->fqctx);
    fq_nmod_poly_clear(tmp2, ctx->fqctx);
    fq_nmod_poly_clear(modulus, ctx->fqctx);
    fq_nmod_poly_clear(tempmod, ctx->fqctx);
    fq_nmod_mpolyu_clear(Aeval, ctx);
    fq_nmod_mpolyu_clear(Beval, ctx);
    fq_nmod_mpolyu_clear(Geval, ctx);
    fq_nmod_mpolyun_clear(An, ctx);
    fq_nmod_mpolyun_clear(Bn, ctx);
    fq_nmod_mpolyun_clear(H, ctx);
    fq_nmod_mpolyun_clear(Ht, ctx);

    return success;
}



int fq_nmod_mpolyu_gcdp_zippel(
    fq_nmod_mpolyu_t G,
    fq_nmod_mpolyu_t Abar,
    fq_nmod_mpolyu_t Bbar,
    fq_nmod_mpolyu_t A,
    fq_nmod_mpolyu_t B,
    slong var,
    const fq_nmod_mpoly_ctx_t ctx,
    flint_rand_t randstate)
{
    slong lastdeg;
    slong Alastdeg;
    slong Blastdeg;
    ulong Ashift, Bshift, Gshift;
    slong degbound;
    slong bound;
    int success = 0, changed, have_enough;
    fq_nmod_mpolyun_t An, Bn;
    fq_nmod_poly_t a, b, c, g, tmp1, tmp2;
    fq_nmod_poly_t modulus, tempmod;
    fq_nmod_mpolyu_t Aeval, Beval, Geval, Abareval, Bbareval, Gform;
    fq_nmod_mpolyun_t H, Ht;
    fq_nmod_t geval, temp;
    fq_nmod_t alpha, alphastart;
    fmpz_t minusone;

    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(var >= -WORD(1));

    FLINT_ASSERT(G->bits == A->bits);
    FLINT_ASSERT(G->bits == B->bits);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    FLINT_ASSERT(fq_nmod_mpolyu_is_canonical(A,ctx));
    FLINT_ASSERT(fq_nmod_mpolyu_is_canonical(B,ctx));

    if (var == -WORD(1))
    {
        /* no more variables left to interpolate */
        return fq_nmod_mpolyu_gcdp_zippel_univar(G, Abar, Bbar, A, B, ctx);
    }

    if (var == WORD(0))
    {
        /* bivariate is more comfortable separated */
        return fq_nmod_mpolyu_gcdp_zippel_bivar(G, Abar, Bbar, A, B, ctx, randstate);
    }

    fq_nmod_mpolyun_init(An, A->bits, ctx);
    fq_nmod_mpolyun_init(Bn, A->bits, ctx);
    fq_nmod_mpolyu_cvtto_mpolyun(An, A, var, ctx);
    fq_nmod_mpolyu_cvtto_mpolyun(Bn, B, var, ctx);

    Ashift = A->exps[A->length - 1];
    Bshift = B->exps[B->length - 1];
    Gshift = FLINT_MIN(Ashift, Bshift);
    fq_nmod_mpolyun_shift_right(An, Ashift);
    fq_nmod_mpolyun_shift_right(Bn, Bshift);

    fq_nmod_poly_init(a, ctx->fqctx);
    fq_nmod_poly_init(b, ctx->fqctx);
    fq_nmod_poly_init(c, ctx->fqctx);
    fq_nmod_poly_init(g, ctx->fqctx);
    fq_nmod_poly_init(tmp1, ctx->fqctx);
    fq_nmod_poly_init(tmp2, ctx->fqctx);

    /* if the gcd has content wrt last variable, we are going to fail */
    fq_nmod_mpolyun_content_poly(a, An, ctx);
    fq_nmod_mpolyun_content_poly(b, Bn, ctx);
    fq_nmod_mpolyun_divexact_poly(An, An, a, ctx);
    fq_nmod_mpolyun_divexact_poly(Bn, Bn, b, ctx);
    fq_nmod_poly_gcd(c, a, b, ctx->fqctx);

    n_fq_poly_get_fq_nmod_poly(tmp1, fq_nmod_mpolyun_leadcoeff_poly(An, ctx), ctx->fqctx);
    n_fq_poly_get_fq_nmod_poly(tmp2, fq_nmod_mpolyun_leadcoeff_poly(Bn, ctx), ctx->fqctx);
    fq_nmod_poly_gcd(g, tmp1, tmp2, ctx->fqctx);

    Alastdeg = fq_nmod_mpolyun_lastdeg(An, ctx);
    Blastdeg = fq_nmod_mpolyun_lastdeg(Bn, ctx);

    /* bound of the number of images required */
    bound = 1 + FLINT_MIN(Alastdeg, Blastdeg)
              + fq_nmod_poly_degree(g, ctx->fqctx);

    /* degree bound on the gcd */
    degbound = FLINT_MIN(A->exps[0], B->exps[0]);

    fq_nmod_poly_init(modulus, ctx->fqctx);
    fq_nmod_poly_init(tempmod, ctx->fqctx);
    fmpz_init(minusone);
    fmpz_set_si(minusone, WORD(-1));
    fq_nmod_poly_set_coeff_fmpz(tempmod, 1, minusone, ctx->fqctx);
    fq_nmod_mpolyu_init(Aeval, A->bits, ctx);
    fq_nmod_mpolyu_init(Beval, A->bits, ctx);
    fq_nmod_mpolyu_init(Geval, A->bits, ctx);
    fq_nmod_mpolyu_init(Abareval, A->bits, ctx);
    fq_nmod_mpolyu_init(Bbareval, A->bits, ctx);
    fq_nmod_mpolyu_init(Gform, A->bits, ctx);
    fq_nmod_mpolyun_init(H, A->bits, ctx);
    fq_nmod_mpolyun_init(Ht, A->bits, ctx);

    fq_nmod_init(geval, ctx->fqctx);
    fq_nmod_init(temp, ctx->fqctx);
    fq_nmod_init(alpha, ctx->fqctx);
    fq_nmod_init(alphastart, ctx->fqctx);

    /* fail if the gcd has content wrt last variable */
    if (fq_nmod_poly_degree(c, ctx->fqctx) > 0)
    {
        success = 0;
        goto finished;
    }

    /* we don't expect this function to work over F_p */
    if (nmod_poly_degree(ctx->fqctx->modulus) < WORD(2))
    {
        success = 0;
        goto finished;
    }

    fq_nmod_poly_one(modulus, ctx->fqctx);
    fq_nmod_mpolyun_zero(H, ctx);

    fq_nmod_randtest_not_zero(alphastart, randstate, ctx->fqctx);
    fq_nmod_set(alpha, alphastart, ctx->fqctx);

    while (1)
    {
        /* get new evaluation point */
        fq_nmod_next_not_zero(alpha, ctx->fqctx);
        if (fq_nmod_equal(alpha, alphastart, ctx->fqctx))
        {
            success = 0;
            goto finished;
        }

        /* make sure evaluation point does not kill both lc(A) and lc(B) */
        fq_nmod_poly_evaluate_fq_nmod(geval, g, alpha, ctx->fqctx);
        if (fq_nmod_is_zero(geval, ctx->fqctx))
        {
            goto outer_continue;
        }

        /* make sure evaluation point does not kill either A or B */
        fq_nmod_mpolyun_interp_reduce_sm_mpolyu(Aeval, An, alpha, ctx);
        fq_nmod_mpolyun_interp_reduce_sm_mpolyu(Beval, Bn, alpha, ctx);
        if (Aeval->length == 0 || Beval->length == 0)
        {
            goto outer_continue;
        }

        success = fq_nmod_mpolyu_gcdp_zippel(Geval, Abareval, Bbareval,
                                        Aeval, Beval, var - 1, ctx, randstate);
        if (!success || Geval->exps[0] > degbound)
        {
            success = 0;
            goto finished;
        }
        
        degbound = Geval->exps[0];

        if (fq_nmod_mpolyu_is_one(Geval, ctx))
        {
            fq_nmod_mpolyu_one(G, ctx);
            fq_nmod_mpolyu_swap(Abar, A, ctx);
            fq_nmod_mpolyu_swap(Bbar, B, ctx);
            fq_nmod_mpolyu_shift_left(G, Gshift);
            fq_nmod_mpolyu_shift_left(Abar, Ashift - Gshift);
            fq_nmod_mpolyu_shift_left(Bbar, Bshift - Gshift);
            success = 1;
            goto finished;
        }

        if (fq_nmod_poly_degree(modulus, ctx->fqctx) > 0)
        {
            if (Geval->exps[0] > H->exps[0])
            {
                goto outer_continue;
            }
            else if (Geval->exps[0] < H->exps[0])
            {
                fq_nmod_poly_one(modulus, ctx->fqctx);                
            }
        }

        /* update interpolant H */
        n_fq_get_fq_nmod(temp, fq_nmod_mpolyu_leadcoeff(Geval, ctx), ctx->fqctx);
        fq_nmod_inv(temp, temp, ctx->fqctx);
        fq_nmod_mul(temp, geval, temp, ctx->fqctx);
        fq_nmod_mpolyu_scalar_mul_fq_nmod(Geval, temp, ctx);
        if (fq_nmod_poly_degree(modulus, ctx->fqctx) > 0)
        {
            fq_nmod_poly_evaluate_fq_nmod(temp, modulus, alpha, ctx->fqctx);
            fq_nmod_inv(temp, temp, ctx->fqctx);
            fq_nmod_poly_scalar_mul_fq_nmod(modulus, modulus, temp, ctx->fqctx);
            changed = fq_nmod_mpolyun_interp_crt_sm_mpolyu(&lastdeg, H, Ht, Geval,
                                                          modulus, alpha, ctx);
            if (!changed)
            {
                fq_nmod_mpolyun_content_poly(a, H, ctx);
                fq_nmod_mpolyun_divexact_poly(Ht, H, a, ctx);
                fq_nmod_mpolyun_shift_left(Ht, Gshift);
                fq_nmod_mpolyu_cvtfrom_mpolyun(G, Ht, var, ctx);
                if (   !fq_nmod_mpolyuu_divides(Abar, A, G, 1, ctx)
                    || !fq_nmod_mpolyuu_divides(Bbar, B, G, 1, ctx))
                {
                    goto outer_continue;
                }
                success = 1;
                goto finished;
            }
        }
        else
        {
            fq_nmod_mpolyun_interp_lift_sm_mpolyu(H, Geval, ctx);
        }
        fq_nmod_poly_set_coeff(tempmod, 0, alpha, ctx->fqctx);
        fq_nmod_poly_mul(modulus, modulus, tempmod, ctx->fqctx);

        fq_nmod_mpolyu_setform_mpolyun(Gform, H, ctx);

        while (1)
        {
            /* get new evaluation point */
            fq_nmod_next_not_zero(alpha, ctx->fqctx);
            if (fq_nmod_equal(alpha, alphastart, ctx->fqctx))
            {
                success = 0;
                goto finished;
            }

            /* make sure evaluation point does not kill both lc(A) and lc(B) */
            fq_nmod_poly_evaluate_fq_nmod(geval, g, alpha, ctx->fqctx);
            if (fq_nmod_is_zero(geval, ctx->fqctx))
            {
                goto inner_continue;
            }

            /* make sure evaluation point does not kill either A or B */
            fq_nmod_mpolyun_interp_reduce_sm_mpolyu(Aeval, An, alpha, ctx);
            fq_nmod_mpolyun_interp_reduce_sm_mpolyu(Beval, Bn, alpha, ctx);
            if (Aeval->length == 0 || Beval->length == 0)
            {
                goto inner_continue;
            }

            switch (fq_nmod_mpolyu_gcds_zippel(Geval, Aeval, Beval, Gform, var,
                                                    ctx, randstate, &degbound))
            {
                default:
                    FLINT_ASSERT(0);
                case nmod_gcds_form_main_degree_too_high:
                    /* fq_nmod_mpolyu_gcds_zippel has updated degbound */
                    fq_nmod_poly_one(modulus, ctx->fqctx);
                    goto outer_continue;
                case nmod_gcds_form_wrong:
                case nmod_gcds_no_solution:
                    success = 0;
                    goto finished;
                case nmod_gcds_scales_not_found:
                case nmod_gcds_eval_point_not_found:
                case nmod_gcds_eval_gcd_deg_too_high:
                    goto inner_continue;
                case nmod_gcds_success:
                    (void)(NULL);
            }

            n_fq_get_fq_nmod(temp, fq_nmod_mpolyu_leadcoeff(Geval, ctx), ctx->fqctx);
            if (fq_nmod_is_zero(temp, ctx->fqctx))
            {
                goto inner_continue;
            }

            /* update interpolant H */
            FLINT_ASSERT(fq_nmod_poly_degree(modulus, ctx->fqctx) > 0);

            fq_nmod_inv(temp, temp, ctx->fqctx);
            fq_nmod_mul(temp, temp, geval, ctx->fqctx);
            fq_nmod_mpolyu_scalar_mul_fq_nmod(Geval, temp, ctx);

            fq_nmod_poly_evaluate_fq_nmod(temp, modulus, alpha, ctx->fqctx);
            fq_nmod_inv(temp, temp, ctx->fqctx);
            fq_nmod_poly_scalar_mul_fq_nmod(modulus, modulus, temp, ctx->fqctx);
            changed = fq_nmod_mpolyun_interp_crt_sm_mpolyu(&lastdeg, H, Ht, Geval, modulus,
                                                                   alpha, ctx);

            fq_nmod_poly_set_coeff(tempmod, 0, alpha, ctx->fqctx);
            fq_nmod_poly_mul(modulus, modulus, tempmod, ctx->fqctx);

            have_enough = fq_nmod_poly_degree(modulus, ctx->fqctx) >= bound;

            if (changed && !have_enough)
            {
                goto inner_continue;
            }

            if (!changed || have_enough)
            {
                fq_nmod_mpolyun_content_poly(a, H, ctx);
                fq_nmod_mpolyun_divexact_poly(Ht, H, a, ctx);
                fq_nmod_mpolyun_shift_left(Ht, Gshift);
                fq_nmod_mpolyu_cvtfrom_mpolyun(G, Ht, var, ctx);
                if (    fq_nmod_mpolyuu_divides(Abar, A, G, 1, ctx)
                     && fq_nmod_mpolyuu_divides(Bbar, B, G, 1, ctx))
                {
                    success = 1;
                    goto finished;
                }
            }

            if (have_enough)
            {
                fq_nmod_poly_one(modulus, ctx->fqctx);
                goto outer_continue;
            }

inner_continue:;
        }
        FLINT_ASSERT(0 && "not reachable");

outer_continue:;
    }
    FLINT_ASSERT(0 && "not reachable");


finished:

    fmpz_clear(minusone);
    fq_nmod_clear(geval, ctx->fqctx);
    fq_nmod_clear(temp, ctx->fqctx);
    fq_nmod_clear(alpha, ctx->fqctx);
    fq_nmod_clear(alphastart, ctx->fqctx);
    fq_nmod_poly_clear(a, ctx->fqctx);
    fq_nmod_poly_clear(b, ctx->fqctx);
    fq_nmod_poly_clear(c, ctx->fqctx);
    fq_nmod_poly_clear(g, ctx->fqctx);
    fq_nmod_poly_clear(tmp1, ctx->fqctx);
    fq_nmod_poly_clear(tmp2, ctx->fqctx);
    fq_nmod_poly_clear(modulus, ctx->fqctx);
    fq_nmod_poly_clear(tempmod, ctx->fqctx);
    fq_nmod_mpolyu_clear(Aeval, ctx);
    fq_nmod_mpolyu_clear(Beval, ctx);
    fq_nmod_mpolyu_clear(Geval, ctx);
    fq_nmod_mpolyu_clear(Abareval, ctx);
    fq_nmod_mpolyu_clear(Bbareval, ctx);
    fq_nmod_mpolyu_clear(Gform, ctx);
    fq_nmod_mpolyun_clear(An, ctx);
    fq_nmod_mpolyun_clear(Bn, ctx);
    fq_nmod_mpolyun_clear(H, ctx);
    fq_nmod_mpolyun_clear(Ht, ctx);

    return success;
}

