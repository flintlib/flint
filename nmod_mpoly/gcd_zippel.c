/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

/* store in each coefficient the evaluation of the corresponding monomial */
void nmod_mpoly_evalsk(nmod_mpoly_t A, nmod_mpoly_t B,
           slong entries, slong * offs, ulong * masks, mp_limb_t * powers,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    slong N;

    FLINT_ASSERT(A->bits == B->bits);
    nmod_mpoly_fit_length(A, B->length, ctx);
    N = mpoly_words_per_exp(B->bits, ctx->minfo);
    for (i = 0; i < B->length; i++)
    {
        mp_limb_t prod = UWORD(1);

        for (j = 0; j < entries; j++)
        {
            if ((B->exps + N*i)[offs[j]] & masks[j])
            {
                prod = nmod_mul(prod, powers[j], ctx->ffinfo->mod);
            }
        }

        A->coeffs[i] = prod;
        mpoly_monomial_set(A->exps + N*i, B->exps + N*i, N);
    }
    A->length = B->length;
}

void nmod_mpolyu_evalsk(nmod_mpolyu_t A, nmod_mpolyu_t B,
              slong entries, slong * offs, ulong * masks, mp_limb_t * powers,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i;

    nmod_mpolyu_fit_length(A, B->length, ctx);
    for (i = 0; i < B->length; i++)
    {
        A->exps[i] = B->exps[i];
        nmod_mpoly_evalsk(A->coeffs + i, B->coeffs + i,
                                            entries, offs, masks, powers, ctx);
    }
    A->length = B->length;
}

/* multiply the coefficients of A pointwise by those of B */
void nmod_mpolyu_mulsk(nmod_mpolyu_t A, nmod_mpolyu_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i, j;

    FLINT_ASSERT(A->length == B->length);
    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT(A->exps[i] == B->exps[i]);

        FLINT_ASSERT((A->coeffs + i)->length == (B->coeffs + i)->length);
        for (j = 0; j < (A->coeffs + i)->length; j++)
        {
            (A->coeffs + i)->coeffs[j] = nmod_mul((A->coeffs + i)->coeffs[j],
                                                  (B->coeffs + i)->coeffs[j],
                                                             ctx->ffinfo->mod);
        }
    }
}


/*
    return 0 if the leading coeff of A vanishes
    else return 1
*/
int nmod_mpolyu_evalfromsk(nmod_poly_t e, nmod_mpolyu_t A,
                                  nmod_mpolyu_t SK, const nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    int ret = 0;

    FLINT_ASSERT(A->length == SK->length);

    nmod_poly_zero(e);
    for (i = 0; i < A->length; i++)
    {
        mp_limb_t v, pp0, pp1, ac0 = 0, ac1 = 0, ac2 = 0;

        FLINT_ASSERT((A->coeffs + i)->length == (SK->coeffs + i)->length);

        for (j = 0; j < (A->coeffs + i)->length; j++)
        {
            umul_ppmm(pp1, pp0, (A->coeffs + i)->coeffs[j],
                                                  (SK->coeffs + i)->coeffs[j]);
            add_sssaaaaaa(ac2, ac1, ac0, ac2, ac1, ac0, WORD(0), pp1, pp0);
        }
        NMOD_RED3(v, ac2, ac1, ac0, ctx->ffinfo->mod);

        nmod_poly_set_coeff_ui(e, A->exps[i], v);
        ret |= (i == 0 && v != 0);
    }

    return ret;
}



/*
    solve

    [ a[0]    a[1]    ... a[n-1]   ]   [ x[0]   ]    [ b[0]   ]
    [ a[0]^2  a[1]^2  ... a[n-1]^2 ]   [ x[1]   ]    [ b[1]   ]
    [  ...                 ...     ] . [ ...    ] =  [ ...    ]
    [ a[0]^n  a[1]^n  ... a[n-1]^n ]   [ x[n-1] ]    [ b[n-1] ]

    for x
*/
int nmod_vandsolve(mp_limb_t * x, mp_limb_t * a, mp_limb_t * b,
                                                           slong n, nmod_t mod)
{
    int success = 0;
    slong i, j;
    mp_limb_t t;
    mp_limb_t Dinv;
    nmod_poly_t Q, P, R, u;

    for (i = 0; i < n; i++)
        x[i] = 0;

    nmod_poly_init(Q, mod.n);
    nmod_poly_init(P, mod.n);
    nmod_poly_init(R, mod.n);
    nmod_poly_init(u, mod.n);
    nmod_poly_set_coeff_ui(u, 1, 1);
    nmod_poly_product_roots_nmod_vec(P, a, n);
    for (i = 0; i < n; i++)
    {
        if (a[i] == UWORD(0))
            goto cleanup;

        nmod_poly_set_coeff_ui(u, 0, mod.n - a[i]);
        nmod_poly_divrem(Q, R, P, u);
        t = nmod_mul(a[i], nmod_poly_evaluate_nmod(Q, a[i]), mod);
        if (t == UWORD(0))
            goto cleanup;

        Dinv = nmod_inv(t, mod);
        for (j = 0; j < n; j++)
        {
            t = nmod_mul(b[j], Dinv, mod);
            t = nmod_mul(t, nmod_poly_get_coeff_ui(Q, j), mod);
            x[i] = nmod_add(x[i], t, mod);
        }
    }
    success = 1;

cleanup:
    nmod_poly_clear(Q);
    nmod_poly_clear(P);
    nmod_poly_clear(R);
    nmod_poly_clear(u);

    return success;
}


/*
    Try to set G to the gcd of A and B given the form f of G.
    return codes as enumerated in nmod_mpoly.h:

    nmod_gcds_success,
    nmod_gcds_form_wrong,
    nmod_gcds_no_solution,
    nmod_gcds_scales_not_found,
    nmod_gcds_eval_point_not_found,
    nmod_gcds_eval_gcd_deg_too_high
*/
nmod_gcds_ret_t nmod_mpolyu_gcds_zippel(nmod_mpolyu_t G,
               nmod_mpolyu_t A, nmod_mpolyu_t B,  nmod_mpolyu_t f, slong var,
          const nmod_mpoly_ctx_t ctx, flint_rand_t randstate, slong * degbound)
{
    int eval_points_tried;
    nmod_gcds_ret_t success;
    nmod_mpolyu_t Aevalsk1, Bevalsk1, fevalsk1, Aevalski, Bevalski, fevalski;
    nmod_poly_t Aeval, Beval, Geval;
    mp_limb_t * alpha, * b;
    nmod_mat_struct * M, * ML;
    nmod_mat_t MF, Msol;
    int lc_ok;
    int underdeterminedcount = 0;
    int exceededcount = 0;
    int * ML_is_initialized;
    slong i, j, k, s, S, nullity;
    slong * d;
    slong l;
    mp_limb_t * W;
    slong entries;
    slong * offs;
    ulong * masks;
    mp_limb_t * powers;
    slong N;
    TMP_INIT;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(f->length > 0);

    FLINT_ASSERT(A->bits == B->bits);
    FLINT_ASSERT(A->bits == G->bits);
    FLINT_ASSERT(A->bits == f->bits);
    FLINT_ASSERT(var >= 0);

    FLINT_ASSERT(*degbound == f->exps[0]);

    FLINT_ASSERT(var > 0);

    if (f->length == 1)
    {
        if ((f->coeffs + 0)->length > 1)
        {
            /* impossible to find scale factors in this case */
            return nmod_gcds_scales_not_found;
        } else {
            /* otherwise set the coeff of the monomial to one */
            nmod_mpolyu_set(G, f, ctx);
            (G->coeffs + 0)->coeffs[0] = UWORD(1);
            if (!nmod_mpolyu_divides(A, G, ctx))
                return nmod_gcds_form_wrong;
            if (!nmod_mpolyu_divides(B, G, ctx))
                return nmod_gcds_form_wrong;
            return nmod_gcds_success;
        }
    }

    TMP_START;

    nmod_mpolyu_init(Aevalsk1, f->bits, ctx);
    nmod_mpolyu_init(Bevalsk1, f->bits, ctx);
    nmod_mpolyu_init(fevalsk1, f->bits, ctx);
    nmod_mpolyu_init(Aevalski, f->bits, ctx);
    nmod_mpolyu_init(Bevalski, f->bits, ctx);
    nmod_mpolyu_init(fevalski, f->bits, ctx);
    nmod_poly_init(Aeval, ctx->ffinfo->mod.n);
    nmod_poly_init(Beval, ctx->ffinfo->mod.n);
    nmod_poly_init(Geval, ctx->ffinfo->mod.n);

    d = (slong *) TMP_ALLOC(f->length*sizeof(slong));
    for (i = 0; i < f->length; i++)
    {
        d[i] = i;
    }

    /*
        make d sort the coeffs so that
        (f->coeffs + d[j-1])->length <= (f->coeffs + d[j-0])->length
        for all j
    */
    for (i = 1; i<f->length; i++)
    {
        for (j=i; j > 0 && (f->coeffs + d[j-1])->length 
                         > (f->coeffs + d[j-0])->length; j--)
        {
            slong temp = d[j-1];
            d[j-1] = d[j-0];
            d[j-0] = temp;
        }
    }

    /* l is the number of images we will try to construct */
    l = f->length - 3;
    for (i = 0; i < f->length; i++)
    {
        l += (f->coeffs + i)->length;
    }
    l = l / (f->length - 1);
    l = FLINT_MAX(l, (f->coeffs + d[f->length - 1])->length);
    /* one extra test image */
    l += 1;

    alpha = (mp_limb_t *) TMP_ALLOC(var*sizeof(mp_limb_t));
    ML = (nmod_mat_struct *) TMP_ALLOC(f->length*sizeof(nmod_mat_struct));
    b = (mp_limb_t *) TMP_ALLOC((f->coeffs + d[f->length - 1])->length
                                                           *sizeof(mp_limb_t));

    nmod_mat_init(MF, 0, l, ctx->ffinfo->mod.n);

    M = (nmod_mat_struct *) TMP_ALLOC(f->length*sizeof(nmod_mat_struct));
    ML_is_initialized = (int *) TMP_ALLOC(f->length*sizeof(int));
    for (i = 0; i < f->length; i++)
    {
        nmod_mat_init(M + i, l, (f->coeffs + i)->length, ctx->ffinfo->mod.n);
        ML_is_initialized[i] = 0;
    }

    W = (mp_limb_t *) flint_malloc(l*f->length*sizeof(mp_limb_t));

    nmod_mat_init(Msol, l, 1, ctx->ffinfo->mod.n);

    /* compute how many masks are needed */
    entries = f->bits * var;
    offs = (slong *) TMP_ALLOC(entries*sizeof(slong));
    masks = (ulong *) TMP_ALLOC(entries*sizeof(slong));
    powers = (mp_limb_t *) TMP_ALLOC(entries*sizeof(mp_limb_t));
    N = mpoly_words_per_exp(f->bits, ctx->minfo);


    /***** evaluation loop head *******/
    eval_points_tried = 0;
pick_evaluation_point:

    if (++eval_points_tried > 10)
    {
        success = nmod_gcds_eval_point_not_found;
        goto finished;
    }

    /* avoid 0, 1 and -1 for the evaluation points */
    FLINT_ASSERT(ctx->ffinfo->mod.n > UWORD(3));
    for (i = 0; i < var; i++)
        alpha[i] = UWORD(2) + n_randint(randstate, ctx->ffinfo->mod.n - UWORD(3));

    /* store bit masks for each power of two of the non-main variables */
    for (i = 0; i < var; i++)
    {
        slong shift, off;
        mpoly_gen_offset_shift(&off, &shift, i, N, f->bits, ctx->minfo);
        for (j = 0; j < f->bits; j++)
        {
            offs[f->bits*i + j] = off;
            masks[f->bits*i + j] = UWORD(1) << (j + shift);
            if (j == 0)
                powers[f->bits*i + j] = alpha[i];
            else
                powers[f->bits*i + j] = nmod_mul(powers[f->bits*i + j-1],
                                                 powers[f->bits*i + j-1],
                                                             ctx->ffinfo->mod);
        }
    }

    nmod_mpolyu_evalsk(Aevalsk1, A, entries, offs, masks, powers, ctx);
    nmod_mpolyu_evalsk(Bevalsk1, B, entries, offs, masks, powers, ctx);
    nmod_mpolyu_evalsk(fevalsk1, f, entries, offs, masks, powers, ctx);

    for (i = 0; i < l*f->length; i++)
    {
        W[i] = 0;
    }


    for (i = 0; i < l; i++)
    {
        if (i == 0)
        {
            nmod_mpolyu_set(Aevalski, Aevalsk1, ctx);
            nmod_mpolyu_set(Bevalski, Bevalsk1, ctx);
            nmod_mpolyu_set(fevalski, fevalsk1, ctx);
        } else
        {
            nmod_mpolyu_mulsk(Aevalski, Aevalsk1, ctx);
            nmod_mpolyu_mulsk(Bevalski, Bevalsk1, ctx);
            nmod_mpolyu_mulsk(fevalski, fevalsk1, ctx);
        }

        for (j = 0; j < f->length; j++)
        {
            for (k = 0; k < (f->coeffs + j)->length; k++)
            {
                (M + j)->rows[i][k] = (fevalski->coeffs + j)->coeffs[k];
            }
        }

        lc_ok = 1;
        lc_ok = lc_ok && nmod_mpolyu_evalfromsk(Aeval, A, Aevalski, ctx);
        lc_ok = lc_ok && nmod_mpolyu_evalfromsk(Beval, B, Bevalski, ctx);
        if (!lc_ok)
        {
            /* lc of A or B vanished */
            goto pick_evaluation_point;
        }

        nmod_poly_gcd(Geval, Aeval, Beval);

        if (f->exps[0] < nmod_poly_degree(Geval))
        {
            ++exceededcount;
            if (exceededcount < 2)
                goto pick_evaluation_point;

            success = nmod_gcds_eval_gcd_deg_too_high;
            goto finished;
        }

        if (f->exps[0] > nmod_poly_degree(Geval))
        {
            success = nmod_gcds_form_main_degree_too_high;
            *degbound = nmod_poly_degree(Geval);
            goto finished;
        }

        k = nmod_poly_length(Geval);
        j = WORD(0);
        while ((--k) >= 0)
        {
            mp_limb_t ck = nmod_poly_get_coeff_ui(Geval, k);
            if (ck != UWORD(0))
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
                W[l*j + i] = ck;
            }
        }
    }

    nullity = -1;
    nmod_mat_clear(MF);
    nmod_mat_init(MF, 0, l, ctx->ffinfo->mod.n);

    for (S = 0; S < f->length; S++)
    {
        s = d[S];

        if (!ML_is_initialized[s])
        {
            nmod_mat_init(ML + s, l, (f->coeffs + s)->length + l,
                                                           ctx->ffinfo->mod.n);
            ML_is_initialized[s] = 1;
            for (i = 0; i < l; i++)
            {
                for (j = 0; j < (f->coeffs + s)->length; j++)
                {
                    (ML + s)->rows[i][j] = (M + s)->rows[i][j];
                }
                (ML + s)->rows[i][(f->coeffs + s)->length + i] = W[l*s + i];
            }
        } else {
            for (i = 0; i < l; i++)
            {
                for (j = 0; j < (f->coeffs + s)->length; j++)
                {
                    (ML + s)->rows[i][j] = (M + s)->rows[i][j];
                }
                for (j = 0; j < l; j++) {
                    (ML + s)->rows[i][(f->coeffs + s)->length + j]
                                             = (j==i ? W[l*s + i] : UWORD(0));
                }
            }

        }
        nmod_mat_rref(ML + s);

        for (i = 0; i < (f->coeffs + s)->length; i++)
        {
            if ((ML + s)->rows[i][i] != UWORD(1))
            {
                /* evaluation points produced a singular vandermonde matrix */
                goto pick_evaluation_point;
            }
        }

        {
            /* appends rows to MF */
            nmod_mat_t MFtemp;
            nmod_mat_t Mwindow;

            nmod_mat_window_init(Mwindow, ML + s,
                    (f->coeffs + s)->length, (f->coeffs + s)->length,
                     l, (f->coeffs + s)->length + l);
            nmod_mat_init(MFtemp,
                             nmod_mat_nrows(MF) + l - (f->coeffs + s)->length,
                             l, ctx->ffinfo->mod.n);
            nmod_mat_concat_vertical(MFtemp, MF, Mwindow);
            nmod_mat_swap(MFtemp, MF);
            nmod_mat_clear(MFtemp);
            nmod_mat_window_clear(Mwindow);
        }

        nullity = l - nmod_mat_rref(MF);

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

    nullity = nmod_mat_nullspace(Msol, MF);
    FLINT_ASSERT(nullity == 1);

    nmod_mpolyu_set(G, f, ctx);

    for (i = 0; i < f->length; i++)
    {
        for (j = 0; j < (f->coeffs + i)->length; j++)
        {
            FLINT_ASSERT((f->coeffs + i)->length <= l);
            b[j] = nmod_mul(W[l*i + j], nmod_mat_get_entry(Msol, j, 0),
                                                             ctx->ffinfo->mod);
        }
        success = nmod_vandsolve((G->coeffs + i)->coeffs,
                                 (fevalsk1->coeffs + i)->coeffs, b, 
                                    (f->coeffs + i)->length, ctx->ffinfo->mod);
        if (!success)
        {
            /* evaluation points produced a singular vandermonde matrix */
            goto pick_evaluation_point;
        }
    }

    /* check solution */
    for (s = 0; s < f->length; s++)
    {
        mp_limb_t pp0, pp1, ac0, ac1, ac2, u, v;

        for (i = 0; i < l; i++)
        {
            ac0 = ac1 = ac2 = 0;
            for (j = 0; j < (f->coeffs + s)->length; j++)
            {
                umul_ppmm(pp1, pp0, (M + s)->rows[i][j],
                                    (G->coeffs + s)->coeffs[j]);
                add_sssaaaaaa(ac2, ac1, ac0, ac2, ac1, ac0, WORD(0), pp1, pp0);
            }

            NMOD_RED3(v, ac2, ac1, ac0, ctx->ffinfo->mod);
            u = nmod_mul(W[l*s + i], nmod_mat_get_entry(Msol, i, 0),
                                                             ctx->ffinfo->mod);
            if (v != u)
            {
                success = nmod_gcds_no_solution;
                goto finished;
            }
        }
    }

    success = nmod_gcds_success;

finished:

    flint_free(W);
    nmod_mat_clear(MF);
    nmod_mat_clear(Msol);
    for (i = 0; i < f->length; i++)
    {
        nmod_mat_clear(M + i);
        if (ML_is_initialized[i])
        {
            nmod_mat_clear(ML + i);
        }
    }
    nmod_mpolyu_clear(Aevalsk1, ctx);
    nmod_mpolyu_clear(Bevalsk1, ctx);
    nmod_mpolyu_clear(fevalsk1, ctx);
    nmod_mpolyu_clear(Aevalski, ctx);
    nmod_mpolyu_clear(Bevalski, ctx);
    nmod_mpolyu_clear(fevalski, ctx);
    nmod_poly_clear(Aeval);
    nmod_poly_clear(Beval);
    nmod_poly_clear(Geval);

    TMP_END;
    return success;
}


/* setform copies the exponents and zeros the coefficients */
void nmod_mpoly_setform(nmod_mpoly_t A, nmod_mpoly_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i;

    nmod_mpoly_set(A, B, ctx);
    for (i = 0; i < A->length; i++)
    {
        A->coeffs[i] = UWORD(0);
    }
}

void nmod_mpolyu_setform(nmod_mpolyu_t A, nmod_mpolyu_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i;

    nmod_mpolyu_fit_length(A, B->length, ctx);
    for (i = 0; i < B->length; i++)
    {
        nmod_mpoly_setform(A->coeffs + i, B->coeffs + i, ctx);
        A->exps[i] = B->exps[i];
    }
    A->length = B->length;
}

void nmod_mpoly_setform_mpolyn(nmod_mpoly_t A, nmod_mpolyn_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong N;

    FLINT_ASSERT(A->bits == B->bits);

    nmod_mpoly_fit_length(A, B->length, ctx);
    N = mpoly_words_per_exp(B->bits, ctx->minfo);
    for (i = 0; i < B->length; i++)
    {
        A->coeffs[i] = UWORD(0);
        mpoly_monomial_set(A->exps + N*i, B->exps + N*i, N);
    }
    A->length = B->length;
}

void nmod_mpolyu_setform_mpolyun(nmod_mpolyu_t A, nmod_mpolyun_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i;

    nmod_mpolyu_fit_length(A, B->length, ctx);
    for (i = 0; i < B->length; i++)
    {
        FLINT_ASSERT((B->coeffs + i)->bits == B->bits);
        nmod_mpoly_setform_mpolyn(A->coeffs + i, B->coeffs + i, ctx);
        A->exps[i] = B->exps[i];
    }
    A->length = B->length;
}

void fq_nmod_mpoly_setform(fq_nmod_mpoly_t A, fq_nmod_mpoly_t B,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;

    fq_nmod_mpoly_set(A, B, ctx);
    for (i = 0; i < A->length; i++)
    {
        fq_nmod_zero(A->coeffs + i, ctx->fqctx);
    }
}

void fq_nmod_mpolyu_setform(fq_nmod_mpolyu_t A, fq_nmod_mpolyu_t B,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;

    fq_nmod_mpolyu_fit_length(A, B->length, ctx);
    for (i = 0; i < B->length; i++)
    {
        fq_nmod_mpoly_setform(A->coeffs + i, B->coeffs + i, ctx);
        A->exps[i] = B->exps[i];
    }
    A->length = B->length;
}


int nmod_mpolyu_gcdp_zippel_univar(nmod_mpolyu_t G,
                  nmod_mpolyu_t A, nmod_mpolyu_t B, const nmod_mpoly_ctx_t ctx)
{
    nmod_poly_t a, b, g;
    FLINT_ASSERT(A->bits == B->bits);
    nmod_poly_init(a, ctx->ffinfo->mod.n);
    nmod_poly_init(b, ctx->ffinfo->mod.n);
    nmod_poly_init(g, ctx->ffinfo->mod.n);
    nmod_mpolyu_cvtto_poly(a, A, ctx);
    nmod_mpolyu_cvtto_poly(b, B, ctx);
    nmod_poly_gcd(g, a, b);
    nmod_mpolyu_cvtfrom_poly(G, g, ctx);
    nmod_poly_clear(a);
    nmod_poly_clear(b);
    nmod_poly_clear(g);
    return 1;
}


int nmod_mpolyu_gcdp_zippel_bivar(nmod_mpolyu_t G,
                                  nmod_mpolyu_t A, nmod_mpolyu_t B,
                             const nmod_mpoly_ctx_t ctx, mpoly_zipinfo_t zinfo)
{
    nmod_poly_t a, b, c, g, modulus, tempmod;
    nmod_mpolyu_t Aeval, Beval, Geval;
    nmod_mpolyun_t An, Bn, H, Ht;
    mp_limb_t geval, temp, alpha;
    int success = 0, changed;
    slong Alastdeg;
    slong Blastdeg;
    slong lastdeg;
    slong var = 0;

    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(var >= -WORD(1));

    FLINT_ASSERT(G->bits == A->bits);
    FLINT_ASSERT(G->bits == B->bits);

    nmod_mpolyun_init(An, A->bits, ctx);
    nmod_mpolyun_init(Bn, A->bits, ctx);

    nmod_mpolyu_cvtto_mpolyun(An, A, var, ctx);
    nmod_mpolyu_cvtto_mpolyun(Bn, B, var, ctx);

    nmod_poly_init(a, ctx->ffinfo->mod.n);
    nmod_poly_init(b, ctx->ffinfo->mod.n);
    nmod_poly_init(c, ctx->ffinfo->mod.n);
    nmod_poly_init(g, ctx->ffinfo->mod.n);
    nmod_mpolyun_content_last(a, An, ctx);
    nmod_mpolyun_content_last(b, Bn, ctx);
    nmod_mpolyun_divexact_last(An, a, ctx);
    nmod_mpolyun_divexact_last(Bn, b, ctx);
    nmod_poly_gcd(c, a, b);
    nmod_poly_gcd(g, nmod_mpolyun_leadcoeff_ref(An, ctx),
                     nmod_mpolyun_leadcoeff_ref(Bn, ctx));

    Alastdeg = nmod_mpolyun_lastdeg(An, ctx);
    Blastdeg = nmod_mpolyun_lastdeg(Bn, ctx);


    nmod_poly_init(modulus, ctx->ffinfo->mod.n);
    nmod_poly_init(tempmod, ctx->ffinfo->mod.n);
    nmod_poly_set_coeff_ui(tempmod, 1, UWORD(1));
    nmod_mpolyu_init(Aeval, A->bits, ctx);
    nmod_mpolyu_init(Beval, A->bits, ctx);
    nmod_mpolyu_init(Geval, A->bits, ctx);
    nmod_mpolyun_init(H, A->bits, ctx);
    nmod_mpolyun_init(Ht, A->bits, ctx);

    nmod_poly_one(modulus);
    nmod_mpolyun_zero(H, ctx);

    alpha = ctx->ffinfo->mod.n;
    while (1)
    {
        if (alpha == UWORD(0))
            goto ret_fail;
        alpha--;

        /* make sure evaluation point does not kill both lc(A) and lc(B) */
        geval = nmod_poly_evaluate_nmod(g, alpha);
        if (geval == WORD(0))
            goto outer_continue;

        /* make sure evaluation point does not kill either A or B */
        nmod_mpolyun_eval_last(Aeval, An, alpha, ctx);
        nmod_mpolyun_eval_last(Beval, Bn, alpha, ctx);
        if (Aeval->length == 0 || Beval->length == 0)
            goto outer_continue;

        nmod_mpolyu_gcdp_zippel_univar(Geval, Aeval, Beval, ctx);


        if (nmod_mpolyu_is_one(Geval, ctx))
        {
            nmod_mpolyu_cvtfrom_poly_notmain(G, c, var, ctx);
            goto ret_success;
        }

        FLINT_ASSERT(Geval->length > 0);

        if (nmod_poly_degree(modulus) > 0)
        {
            if (Geval->exps[0] > H->exps[0])
            {
                goto outer_continue;

            } else if (Geval->exps[0] < H->exps[0])
            {
                nmod_poly_one(modulus);                
            }
        }

        temp = n_invmod(nmod_mpolyu_leadcoeff(Geval, ctx), ctx->ffinfo->mod.n);
        temp = nmod_mul(geval, temp, ctx->ffinfo->mod);
        nmod_mpolyu_scalar_mul_nmod(Geval, temp, ctx);

        /* update interpolant H */
        if (nmod_poly_degree(modulus) > 0)
        {
            nmod_poly_scalar_mul_nmod(modulus, modulus,
                     n_invmod(nmod_poly_evaluate_nmod(modulus, alpha),
                                                          ctx->ffinfo->mod.n));
            changed = nmod_mpolyun_addinterp(&lastdeg, H, Ht, Geval,
                                                          modulus, alpha, ctx);
            if (!changed)
            {
                nmod_mpolyun_content_last(a, H, ctx);
                nmod_mpolyun_mul_poly(Ht, H, c, ctx);
                nmod_mpolyun_divexact_last(Ht, a, ctx);

                nmod_mpolyu_cvtfrom_mpolyun(G, Ht, var, ctx);

                if (!nmod_mpolyu_divides(A, G, ctx))
                    goto outer_continue;
                if (!nmod_mpolyu_divides(B, G, ctx))
                    goto outer_continue;
                goto ret_success;
            }

        } else
        {
            nmod_mpolyun_set_mpolyu(H, Geval, ctx);
            lastdeg = WORD(0);
        }
        nmod_poly_set_coeff_ui(tempmod, 0, ctx->ffinfo->mod.n - alpha);
        nmod_poly_mul(modulus, modulus, tempmod);

        /* something is wrong if the interpolation degree is too high */
        if (lastdeg > Alastdeg || lastdeg > Blastdeg)
        {
            nmod_poly_one(modulus);
            goto outer_continue;
        }

outer_continue:
        NULL;
    }

ret_success:
    success = 1;

ret_fail:
    nmod_poly_clear(a);
    nmod_poly_clear(b);
    nmod_poly_clear(c);
    nmod_poly_clear(g);
    nmod_poly_clear(modulus);
    nmod_poly_clear(tempmod);
    nmod_mpolyu_clear(Aeval, ctx);
    nmod_mpolyu_clear(Beval, ctx);
    nmod_mpolyu_clear(Geval, ctx);
    nmod_mpolyun_clear(An, ctx);
    nmod_mpolyun_clear(Bn, ctx);
    nmod_mpolyun_clear(H, ctx);
    nmod_mpolyun_clear(Ht, ctx);

    return success;
}


int nmod_mpolyu_gcdp_zippel(nmod_mpolyu_t G,
                     nmod_mpolyu_t A, nmod_mpolyu_t B, slong var,
                           const nmod_mpoly_ctx_t ctx, mpoly_zipinfo_t zinfo,
                                                        flint_rand_t randstate)
{
    int divcheck_fail_count;
    slong lastdeg;
    slong Alastdeg;
    slong Blastdeg;
    ulong ABminshift;
    slong degbound;
    int success = 0, changed;
    nmod_mpolyun_t An, Bn;
    nmod_poly_t a, b, c, g;
    nmod_poly_t modulus, tempmod;
    nmod_mpolyu_t Aeval, Beval, Geval, Gform;
    nmod_mpolyun_t H, Ht;
    mp_limb_t geval, temp;
    mp_limb_t alpha, start_alpha;

    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(var >= -WORD(1));

    FLINT_ASSERT(G->bits == A->bits);
    FLINT_ASSERT(G->bits == B->bits);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    if (var == -WORD(1))
    {
        /* no more variables left to interpolate */
        return nmod_mpolyu_gcdp_zippel_univar(G, A, B, ctx);
    }

    nmod_mpolyun_init(An, A->bits, ctx);
    nmod_mpolyun_init(Bn, A->bits, ctx);

    nmod_mpolyu_cvtto_mpolyun(An, A, var, ctx);
    nmod_mpolyu_cvtto_mpolyun(Bn, B, var, ctx);
    ABminshift = FLINT_MIN(A->exps[A->length - 1], B->exps[B->length - 1]);
    nmod_mpolyun_shift_right(An, A->exps[A->length - 1]);
    nmod_mpolyun_shift_right(Bn, B->exps[B->length - 1]);
    Alastdeg = nmod_mpolyun_lastdeg(An, ctx);
    Blastdeg = nmod_mpolyun_lastdeg(Bn, ctx);

    nmod_poly_init(a, ctx->ffinfo->mod.n);
    nmod_poly_init(b, ctx->ffinfo->mod.n);
    nmod_poly_init(c, ctx->ffinfo->mod.n);
    nmod_poly_init(g, ctx->ffinfo->mod.n);

    nmod_mpolyun_content_last(a, An, ctx);
    nmod_mpolyun_content_last(b, Bn, ctx);

    /* if the gcd has content wrt last variable, we are going to fail */
    /*nmod_mpolyun_divexact_last(An, a, ctx);*/
    /*nmod_mpolyun_divexact_last(Bn, b, ctx);*/

    nmod_poly_gcd(c, a, b);
    nmod_poly_gcd(g, nmod_mpolyun_leadcoeff_ref(An, ctx),
                     nmod_mpolyun_leadcoeff_ref(Bn, ctx));

    nmod_poly_init(modulus, ctx->ffinfo->mod.n);
    nmod_poly_init(tempmod, ctx->ffinfo->mod.n);
    nmod_poly_set_coeff_ui(tempmod, 1, UWORD(1));
    nmod_mpolyu_init(Aeval, A->bits, ctx);
    nmod_mpolyu_init(Beval, A->bits, ctx);
    nmod_mpolyu_init(Geval, A->bits, ctx);
    nmod_mpolyu_init(Gform, A->bits, ctx);
    nmod_mpolyun_init(H, A->bits, ctx);
    nmod_mpolyun_init(Ht, A->bits, ctx);

    if (nmod_poly_degree(c) > 0)
    {
        success = 0;
        goto finished;
    }

    if (var == WORD(0))
    {
        /* bivariate is more comfortable separated */
        success = nmod_mpolyu_gcdp_zippel_bivar(G, A, B, ctx, zinfo);
        goto finished;
    }

    if (ctx->ffinfo->mod.n <= UWORD(3))
    {
        success = 0;
        goto finished;
    }

    degbound = FLINT_MIN(A->exps[0], B->exps[0]);

    nmod_poly_one(modulus);
    nmod_mpolyun_zero(H, ctx);

    start_alpha = UWORD(1) + n_randint(randstate, ctx->ffinfo->mod.n - UWORD(1));
    alpha = start_alpha;
    while (1)
    {
        /* get new evaluation point */
        if (!(--alpha))
            alpha += ctx->ffinfo->mod.n - UWORD(1);
        if (alpha == start_alpha)
        {
            success = 0;
            goto finished;
        }

        /* make sure evaluation point does not kill both lc(A) and lc(B) */
        geval = nmod_poly_evaluate_nmod(g, alpha);
        if (geval == WORD(0))
            goto outer_continue;

        /* make sure evaluation point does not kill either A or B */
        nmod_mpolyun_eval_last(Aeval, An, alpha, ctx);
        nmod_mpolyun_eval_last(Beval, Bn, alpha, ctx);
        if (Aeval->length == 0 || Beval->length == 0)
            goto outer_continue;

        success = nmod_mpolyu_gcdp_zippel(Geval, Aeval, Beval, var - 1,
                                                        ctx, zinfo, randstate);
        if (!success || Geval->exps[0] > degbound)
        {
            success = 0;
            goto finished;
        }
        
        degbound = Geval->exps[0];

        if (nmod_mpolyu_is_one(Geval, ctx))
        {
            nmod_mpolyu_cvtfrom_poly_notmain(G, c, var, ctx);
            nmod_mpolyu_shift_left(G, ABminshift);
            success = 1;
            goto finished;
        }

        if (nmod_poly_degree(modulus) > 0)
        {
            if (Geval->exps[0] > H->exps[0])
            {
                goto outer_continue;

            } else if (Geval->exps[0] < H->exps[0])
            {
                nmod_poly_one(modulus);                
            }
        }

        /* update interpolant H */
        temp = n_invmod(nmod_mpolyu_leadcoeff(Geval, ctx), ctx->ffinfo->mod.n);
        temp = nmod_mul(geval, temp, ctx->ffinfo->mod);
        nmod_mpolyu_scalar_mul_nmod(Geval, temp, ctx);
        if (nmod_poly_degree(modulus) > 0)
        {
            temp = nmod_poly_evaluate_nmod(modulus, alpha);
            temp = n_invmod(temp, ctx->ffinfo->mod.n);
            nmod_poly_scalar_mul_nmod(modulus, modulus, temp);
            changed = nmod_mpolyun_addinterp(&lastdeg, H, Ht, Geval,
                                                          modulus, alpha, ctx);
            if (!changed)
            {
                nmod_mpolyun_content_last(a, H, ctx);
                nmod_mpolyun_mul_poly(Ht, H, c, ctx);
                nmod_mpolyun_divexact_last(Ht, a, ctx);
                nmod_mpolyun_shift_left(Ht, ABminshift);
                nmod_mpolyu_cvtfrom_mpolyun(G, Ht, var, ctx);
                if (!nmod_mpolyu_divides(A, G, ctx))
                    goto outer_continue;
                if (!nmod_mpolyu_divides(B, G, ctx))
                    goto outer_continue;
                success = 1;
                goto finished;
            }

        } else
        {
            nmod_mpolyun_set_mpolyu(H, Geval, ctx);
        }
        nmod_poly_set_coeff_ui(tempmod, 0, ctx->ffinfo->mod.n - alpha);
        nmod_poly_mul(modulus, modulus, tempmod);

        nmod_mpolyu_setform_mpolyun(Gform, H, ctx);

        divcheck_fail_count = 0;
        while (1)
        {
            /* get new evaluation point */
            if (!(--alpha))
                alpha += ctx->ffinfo->mod.n - UWORD(1);
            if (alpha == start_alpha)
            {
                success = 0;
                goto finished;
            }

            /* make sure evaluation does not kill both lc(A) and lc(B) */
            geval = nmod_poly_evaluate_nmod(g, alpha);
            if (geval == WORD(0))
                goto inner_continue;

            /* make sure evaluation does not kill either A or B */
            nmod_mpolyun_eval_last(Aeval, An, alpha, ctx);
            nmod_mpolyun_eval_last(Beval, Bn, alpha, ctx);
            if (Aeval->length == 0 || Beval->length == 0)
                goto inner_continue;

            switch (nmod_mpolyu_gcds_zippel(Geval, Aeval, Beval, Gform, var,
                                                    ctx, randstate, &degbound))
            {
                default:
                    FLINT_ASSERT(0);
                case nmod_gcds_form_main_degree_too_high:
                    /* nmod_mpolyu_gcds_zippel has updated degbound */
                    nmod_poly_one(modulus);
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

            if (nmod_mpolyu_leadcoeff(Geval, ctx) == UWORD(0))
                goto inner_continue;

            /* update interpolant H */
            temp = n_invmod(nmod_mpolyu_leadcoeff(Geval, ctx), ctx->ffinfo->mod.n);
            nmod_mpolyu_scalar_mul_nmod(Geval, nmod_mul(geval, temp,
                                                       ctx->ffinfo->mod), ctx);
            FLINT_ASSERT(nmod_poly_degree(modulus) > 0);
            temp = nmod_poly_evaluate_nmod(modulus, alpha);
            temp = n_invmod(temp, ctx->ffinfo->mod.n);
            nmod_poly_scalar_mul_nmod(modulus, modulus, temp);

            changed = nmod_mpolyun_addinterp(&lastdeg, H, Ht, Geval,
                                                          modulus, alpha, ctx);
            nmod_poly_set_coeff_ui(tempmod, 0, ctx->ffinfo->mod.n - alpha);
            nmod_poly_mul(modulus, modulus, tempmod);

            if (!changed)
            {
                nmod_mpolyun_content_last(a, H, ctx);
                nmod_mpolyun_mul_poly(Ht, H, c, ctx);
                nmod_mpolyun_divexact_last(Ht, a, ctx);
                nmod_mpolyun_shift_left(Ht, ABminshift);
                nmod_mpolyu_cvtfrom_mpolyun(G, Ht, var, ctx);
                if (!nmod_mpolyu_divides(A, G, ctx)
                             || !nmod_mpolyu_divides(B, G, ctx))
                {
                    ++divcheck_fail_count;
                    if (divcheck_fail_count >= 2)
                    {
                        nmod_poly_one(modulus);
                        goto outer_continue;
                    }
                    goto inner_continue;
                }
                success = 1;
                goto finished;
            }

            /* something is wrong if the interpolation degree is too high */
            if (lastdeg > Alastdeg || lastdeg > Blastdeg)
            {
                nmod_poly_one(modulus);
                goto outer_continue;
            }
inner_continue:
            NULL;
        }
        FLINT_ASSERT(0 && "not reachable");

outer_continue:
        NULL;
    }
    FLINT_ASSERT(0 && "not reachable");


finished:

    nmod_poly_clear(a);
    nmod_poly_clear(b);
    nmod_poly_clear(c);
    nmod_poly_clear(g);
    nmod_poly_clear(modulus);
    nmod_poly_clear(tempmod);
    nmod_mpolyu_clear(Aeval, ctx);
    nmod_mpolyu_clear(Beval, ctx);
    nmod_mpolyu_clear(Geval, ctx);
    nmod_mpolyu_clear(Gform, ctx);
    nmod_mpolyun_clear(An, ctx);
    nmod_mpolyun_clear(Bn, ctx);
    nmod_mpolyun_clear(H, ctx);
    nmod_mpolyun_clear(Ht, ctx);

    return success;
}


int nmod_mpolyu_gcdm_zippel_bivar(
    nmod_mpolyu_t G,
    nmod_mpolyu_t A,
    nmod_mpolyu_t B,
    nmod_mpoly_ctx_t ctx,
    mpoly_zipinfo_t zinfo,
    flint_rand_t randstate)
{
    slong deg = 2;
    nmod_mpolyun_t An, Bn, H, Ht;
    fq_nmod_mpoly_ctx_t ffctx;
    fq_nmod_mpolyu_t Aeval, Beval, Geval;
    nmod_poly_t modulus, a,b,c,g;
    fq_nmod_t t, geval;
    int success = 0, changed;
    slong Alastdeg;
    slong Blastdeg;
    slong lastdeg;
    slong var = 0;

    FLINT_ASSERT(G->bits == A->bits);
    FLINT_ASSERT(B->bits == A->bits);

    nmod_mpolyun_init(An, A->bits, ctx);
    nmod_mpolyun_init(Bn, A->bits, ctx);

    nmod_mpolyu_cvtto_mpolyun(An, A, var, ctx);
    nmod_mpolyu_cvtto_mpolyun(Bn, B, var, ctx);

    nmod_poly_init(a, ctx->ffinfo->mod.n);
    nmod_poly_init(b, ctx->ffinfo->mod.n);
    nmod_poly_init(c, ctx->ffinfo->mod.n);
    nmod_poly_init(g, ctx->ffinfo->mod.n);
    nmod_mpolyun_content_last(a, An, ctx);
    nmod_mpolyun_content_last(b, Bn, ctx);
    nmod_mpolyun_divexact_last(An, a, ctx);
    nmod_mpolyun_divexact_last(Bn, b, ctx);
    nmod_poly_gcd(c, a, b);
    nmod_poly_gcd(g, nmod_mpolyun_leadcoeff_ref(An, ctx),
                     nmod_mpolyun_leadcoeff_ref(Bn, ctx));

    Alastdeg = nmod_mpolyun_lastdeg(An, ctx);
    Blastdeg = nmod_mpolyun_lastdeg(Bn, ctx);

    nmod_poly_init(modulus, ctx->ffinfo->mod.n);
    nmod_poly_one(modulus);
    nmod_mpolyun_init(H, A->bits, ctx);
    nmod_mpolyun_init(Ht, A->bits, ctx);

    fq_nmod_mpoly_ctx_init(ffctx, ctx->minfo->nvars, ctx->ffinfo->mod.n, deg);

    fq_nmod_mpolyu_init(Aeval, A->bits, ffctx);
    fq_nmod_mpolyu_init(Beval, A->bits, ffctx);
    fq_nmod_mpolyu_init(Geval, A->bits, ffctx);
    fq_nmod_init(geval, ffctx->fqctx);
    fq_nmod_init(t, ffctx->fqctx);


    FLINT_ASSERT(An->bits == B->bits);
    FLINT_ASSERT(An->bits == G->bits);
    FLINT_ASSERT(An->length > 0);
    FLINT_ASSERT(Bn->length > 0);
    FLINT_ASSERT(An->exps[A->length - 1] == 0);
    FLINT_ASSERT(Bn->exps[B->length - 1] == 0);

    while (1) {

        deg++;
        if (deg > 1000)
        {
            /* ran out of primes */
            success = 0;
            goto finished;
        }

        fq_nmod_mpolyu_clear(Aeval, ffctx);
        fq_nmod_mpolyu_clear(Beval, ffctx);
        fq_nmod_mpolyu_clear(Geval, ffctx);
        fq_nmod_clear(geval, ffctx->fqctx);
        fq_nmod_clear(t, ffctx->fqctx);

        fq_nmod_mpoly_ctx_change_modulus(ffctx, deg);

        fq_nmod_mpolyu_init(Aeval, A->bits, ffctx);
        fq_nmod_mpolyu_init(Beval, A->bits, ffctx);
        fq_nmod_mpolyu_init(Geval, A->bits, ffctx);
        fq_nmod_init(geval, ffctx->fqctx);
        fq_nmod_init(t, ffctx->fqctx);

        /* make sure reduction does not kill both lc(A) and lc(B) */
        nmod_poly_rem(geval, g, ffctx->fqctx->modulus);
        if (fq_nmod_is_zero(geval, ffctx->fqctx))
            goto outer_continue;

        /* make sure reduction does not kill either A or B */
        nmod_mpolyun_redto_fq_nmod_mpolyu(Aeval, An, ffctx, ctx);
        nmod_mpolyun_redto_fq_nmod_mpolyu(Beval, Bn, ffctx, ctx);
        if (Aeval->length == 0 || Beval->length == 0)
            goto outer_continue;

        FLINT_ASSERT(fq_nmod_mpolyu_is_canonical(Aeval, ffctx));
        FLINT_ASSERT(fq_nmod_mpolyu_is_canonical(Beval, ffctx));

        fq_nmod_mpolyu_gcdp_zippel_univar(Geval, Aeval, Beval, ffctx);

        if (fq_nmod_mpolyu_is_one(Geval, ffctx))
        {
            nmod_mpolyu_cvtfrom_poly_notmain(G, c, var, ctx);
            success = 1;
            goto finished;
        }

        FLINT_ASSERT(Geval->length > 0);

        if (nmod_poly_degree(modulus) > 0)
        {
            if (Geval->exps[0] > H->exps[0])
            {
                goto outer_continue;

            } else if (Geval->exps[0] < H->exps[0])
            {
                nmod_poly_one(modulus);                
            }
        }

        fq_nmod_inv(t, fq_nmod_mpolyu_leadcoeff_ref(Geval, ffctx), ffctx->fqctx);
        fq_nmod_mul(t, t, geval, ffctx->fqctx);
        fq_nmod_mpolyu_scalar_mul_fq_nmod(Geval, t, ffctx);

        if (nmod_poly_degree(modulus) > 0)
        {
            changed = nmod_mpolyun_addinterp_fq_nmod_mpolyu(&lastdeg, H, Ht,
                                                   modulus, ctx, Geval, ffctx);
            if (!changed)
            {
                nmod_mpolyun_content_last(a, H, ctx);
                nmod_mpolyun_mul_poly(Ht, H, c, ctx);
                nmod_mpolyun_divexact_last(Ht, a, ctx);

                nmod_mpolyu_cvtfrom_mpolyun(G, Ht, var, ctx);

                if (!nmod_mpolyu_divides(A, G, ctx))
                    goto outer_continue;
                if (!nmod_mpolyu_divides(B, G, ctx))
                    goto outer_continue;
                success = 1;
                goto finished;
            }

        } else
        {
            nmod_mpolyun_set_fq_nmod_mpolyu(H, ctx, Geval, ffctx);
            lastdeg = nmod_mpolyun_lastdeg(H, ctx);
        }
        nmod_poly_mul(modulus, modulus, ffctx->fqctx->modulus);

        if (lastdeg > Alastdeg || lastdeg > Blastdeg)
        {
            nmod_poly_one(modulus);
            goto outer_continue;
        }

outer_continue:
        NULL;
    }

finished:

    nmod_poly_clear(a);
    nmod_poly_clear(b);
    nmod_poly_clear(c);
    nmod_poly_clear(g);
    nmod_poly_clear(modulus);
    nmod_mpolyun_clear(An, ctx);
    nmod_mpolyun_clear(Bn, ctx);
    nmod_mpolyun_clear(H, ctx);
    nmod_mpolyun_clear(Ht, ctx);

    fq_nmod_clear(geval, ffctx->fqctx);
    fq_nmod_clear(t, ffctx->fqctx);
    fq_nmod_mpolyu_clear(Aeval, ffctx);
    fq_nmod_mpolyu_clear(Beval, ffctx);
    fq_nmod_mpolyu_clear(Geval, ffctx);
    fq_nmod_mpoly_ctx_clear(ffctx);

    return success;
}



int nmod_mpolyu_gcdm_zippel(
    nmod_mpolyu_t G,
    nmod_mpolyu_t A,
    nmod_mpolyu_t B,
    nmod_mpoly_ctx_t ctx,
    mpoly_zipinfo_t zinfo,
    flint_rand_t randstate)
{
    slong degbound;
    slong coeff_deg_bound;
    slong lastdeg;
    int success, changed;
    slong deg = 4;
    nmod_mpolyun_t An, Bn, Hn, Ht;
    fq_nmod_mpoly_ctx_t ffctx;
    fq_nmod_mpolyu_t Aff, Bff, Gff, Gform;
    nmod_poly_t modulus, gamma, hc;
    fq_nmod_t t, gammaff;

    FLINT_ASSERT(G->bits == A->bits);
    FLINT_ASSERT(B->bits == A->bits);


    success = nmod_mpolyu_gcdp_zippel(G, A, B, ctx->minfo->nvars - 1, ctx, zinfo, randstate);
    if (success)
        return 1;

    /* bivariate more comfortable separated */
    if (ctx->minfo->nvars == 1)
        return nmod_mpolyu_gcdm_zippel_bivar(G, A, B, ctx, zinfo, randstate);

    FLINT_ASSERT(ctx->minfo->nvars > 1);

    nmod_poly_init(hc, ctx->ffinfo->mod.n);
    nmod_poly_init(modulus, ctx->ffinfo->mod.n);

    nmod_mpolyun_init(An, A->bits, ctx);
    nmod_mpolyun_init(Bn, A->bits, ctx);
    nmod_mpolyu_cvtto_mpolyun(An, A, ctx->minfo->nvars - 1, ctx);
    nmod_mpolyu_cvtto_mpolyun(Bn, B, ctx->minfo->nvars - 1, ctx);

    coeff_deg_bound = FLINT_MIN(nmod_mpolyun_lastdeg(An, ctx),
                                nmod_mpolyun_lastdeg(Bn, ctx));

    FLINT_ASSERT(An->bits == B->bits);
    FLINT_ASSERT(An->bits == G->bits);
    FLINT_ASSERT(An->length > 0);
    FLINT_ASSERT(Bn->length > 0);
    FLINT_ASSERT(An->exps[A->length - 1] == 0);
    FLINT_ASSERT(Bn->exps[B->length - 1] == 0);

    nmod_poly_one(modulus);

    nmod_poly_init(gamma, ctx->ffinfo->mod.n);
    nmod_poly_gcd(gamma, nmod_mpolyun_leadcoeff_ref(An, ctx),
                         nmod_mpolyun_leadcoeff_ref(Bn, ctx));

    degbound = FLINT_MIN(A->exps[0], B->exps[0]);

    nmod_mpolyun_init(Hn, A->bits, ctx);
    nmod_mpolyun_init(Ht, A->bits, ctx);


    fq_nmod_mpoly_ctx_init(ffctx, ctx->minfo->nvars, ctx->ffinfo->mod.n, deg);

    fq_nmod_mpolyu_init(Aff, A->bits, ffctx);
    fq_nmod_mpolyu_init(Bff, A->bits, ffctx);
    fq_nmod_mpolyu_init(Gff, A->bits, ffctx);
    fq_nmod_mpolyu_init(Gform, A->bits, ffctx);
    fq_nmod_init(gammaff, ffctx->fqctx);
    fq_nmod_init(t, ffctx->fqctx);


choose_prime_outer:

    deg++;
    if (deg > 1000)
    {
        /* ran out of primes */
        success = 0;
        goto finished;
    }

    fq_nmod_mpolyu_clear(Aff, ffctx);
    fq_nmod_mpolyu_clear(Bff, ffctx);
    fq_nmod_mpolyu_clear(Gff, ffctx);
    fq_nmod_mpolyu_clear(Gform, ffctx);
    fq_nmod_clear(gammaff, ffctx->fqctx);
    fq_nmod_clear(t, ffctx->fqctx);

    fq_nmod_mpoly_ctx_change_modulus(ffctx, deg);

    fq_nmod_mpolyu_init(Aff, A->bits, ffctx);
    fq_nmod_mpolyu_init(Bff, A->bits, ffctx);
    fq_nmod_mpolyu_init(Gff, A->bits, ffctx);
    fq_nmod_mpolyu_init(Gform, A->bits, ffctx);
    fq_nmod_init(gammaff, ffctx->fqctx);
    fq_nmod_init(t, ffctx->fqctx);


    /* make sure reduction does not kill both lc(A) and lc(B) */
    nmod_poly_rem(gammaff, gamma, ffctx->fqctx->modulus);
    if (fq_nmod_is_zero(gammaff, ffctx->fqctx))
        goto choose_prime_outer;

    /* make sure reduction does not kill either A or B */
    nmod_mpolyun_redto_fq_nmod_mpolyu(Aff, An, ffctx, ctx);
    nmod_mpolyun_redto_fq_nmod_mpolyu(Bff, Bn, ffctx, ctx);
    if (Aff->length == 0 || Bff->length == 0)
        goto choose_prime_outer;

    success = fq_nmod_mpolyu_gcdp_zippel(Gff, Aff, Bff, ctx->minfo->nvars - 2,
                                                      ffctx, zinfo, randstate);
    if (!success || Gff->exps[0] > degbound)
        goto choose_prime_outer;
    degbound = Gff->exps[0];

    if (Gff->length == 1 && Gff->exps[0] == 0)
    {
        FLINT_ASSERT(fq_nmod_mpoly_is_one(Gff->coeffs + 0, ffctx));
        FLINT_ASSERT(!fq_nmod_mpoly_is_zero(Gff->coeffs + 0, ffctx));
        nmod_mpolyu_one(G, ctx);
        
        success = 1;
        goto finished;
    }

    fq_nmod_inv(t, fq_nmod_mpolyu_leadcoeff_ref(Gff, ffctx), ffctx->fqctx);
    fq_nmod_mul(t, t, gammaff, ffctx->fqctx);
    fq_nmod_mpolyu_scalar_mul_fq_nmod(Gff, t, ffctx);

    fq_nmod_mpolyu_setform(Gform, Gff, ffctx);
    nmod_mpolyun_set_fq_nmod_mpolyu(Hn, ctx, Gff, ffctx);

    nmod_poly_set(modulus, ffctx->fqctx->modulus);

choose_prime_inner:

    deg++;
    if (deg > 1000)
    {
        /* ran out of primes */
        success = 0;
        goto finished;
    }

    fq_nmod_mpolyu_clear(Aff, ffctx);
    fq_nmod_mpolyu_clear(Bff, ffctx);
    fq_nmod_mpolyu_clear(Gff, ffctx);
    fq_nmod_clear(gammaff, ffctx->fqctx);
    fq_nmod_clear(t, ffctx->fqctx);

    fq_nmod_mpoly_ctx_change_modulus(ffctx, deg);

    fq_nmod_mpolyu_init(Aff, A->bits, ffctx);
    fq_nmod_mpolyu_init(Bff, A->bits, ffctx);
    fq_nmod_mpolyu_init(Gff, A->bits, ffctx);
    fq_nmod_init(gammaff, ffctx->fqctx);
    fq_nmod_init(t, ffctx->fqctx);

    /* make sure reduction does not kill both lc(A) and lc(B) */
    nmod_poly_rem(gammaff, gamma, ffctx->fqctx->modulus);
    if (fq_nmod_is_zero(gammaff, ffctx->fqctx))
        goto choose_prime_inner;

    /* make sure reduction does not kill either A or B */
    nmod_mpolyun_redto_fq_nmod_mpolyu(Aff, An, ffctx, ctx);
    nmod_mpolyun_redto_fq_nmod_mpolyu(Bff, Bn, ffctx, ctx);
    if (Aff->length == 0 || Bff->length == 0)
        goto choose_prime_inner;

    switch (fq_nmod_mpolyu_gcds_zippel(Gff, Aff, Bff, Gform,
                           ctx->minfo->nvars - 1, ffctx, randstate, &degbound))
    {
        default:
            FLINT_ASSERT(0);
        case nmod_gcds_form_main_degree_too_high:
        case nmod_gcds_form_wrong:
        case nmod_gcds_no_solution:
            goto choose_prime_outer;
        case nmod_gcds_scales_not_found:
        case nmod_gcds_eval_point_not_found:
        case nmod_gcds_eval_gcd_deg_too_high:
            goto choose_prime_inner;
        case nmod_gcds_success:
            NULL;
    }

    if (fq_nmod_is_zero(fq_nmod_mpolyu_leadcoeff_ref(Gff, ffctx), ffctx->fqctx))
        goto choose_prime_inner;

    fq_nmod_inv(t, fq_nmod_mpolyu_leadcoeff_ref(Gff, ffctx), ffctx->fqctx);
    fq_nmod_mul(t, t, gammaff, ffctx->fqctx);
    fq_nmod_mpolyu_scalar_mul_fq_nmod(Gff, t, ffctx);

    changed = nmod_mpolyun_CRT_fq_nmod_mpolyu(&lastdeg, Hn, ctx, modulus, Gff, ffctx);
    nmod_poly_mul(modulus, modulus, ffctx->fqctx->modulus);
    if (changed)
    {
        if (lastdeg > coeff_deg_bound) {
            goto choose_prime_outer;
        }
        goto choose_prime_inner;
    }

    nmod_mpolyun_content_last(hc, Hn, ctx);
    nmod_mpolyun_set(Ht, Hn, ctx);
    nmod_mpolyun_divexact_last(Ht, hc, ctx);
    nmod_mpolyu_cvtfrom_mpolyun(G, Ht, ctx->minfo->nvars - 1, ctx);

    if (!nmod_mpolyu_divides(A, G, ctx) || !nmod_mpolyu_divides(B, G, ctx))
        goto choose_prime_inner;

    success = 1;
    goto finished;

finished:

    nmod_poly_clear(gamma);

    nmod_poly_clear(hc);
    nmod_poly_clear(modulus);

    nmod_mpolyun_clear(An, ctx);
    nmod_mpolyun_clear(Bn, ctx);
    nmod_mpolyun_clear(Hn, ctx);
    nmod_mpolyun_clear(Ht, ctx);

    fq_nmod_mpolyu_clear(Aff, ffctx);
    fq_nmod_mpolyu_clear(Bff, ffctx);
    fq_nmod_mpolyu_clear(Gff, ffctx);
    fq_nmod_mpolyu_clear(Gform, ffctx);
    fq_nmod_clear(gammaff, ffctx->fqctx);
    fq_nmod_clear(t, ffctx->fqctx);

    fq_nmod_mpoly_ctx_clear(ffctx);

    return success;
}


int nmod_mpolyu_gcd_zippel(
    nmod_mpolyu_t G,
    nmod_mpolyu_t A,
    nmod_mpolyu_t B,
    nmod_mpoly_ctx_t ctx,
    mpoly_zipinfo_t zinfo,
    flint_rand_t randstate)
{
    int success = 0;
    slong i;
    slong ABminshift;
    nmod_mpoly_t content;
    nmod_mpolyu_t Abar, Bbar, Gbar;

    FLINT_ASSERT(A->bits == B->bits);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    nmod_mpoly_init(content, ctx);
    nmod_mpolyu_init(Abar, A->bits, ctx);
    nmod_mpolyu_init(Bbar, A->bits, ctx);
    nmod_mpolyu_init(Gbar, A->bits, ctx);

    /* compute the content of GCD wrt non main variables */
    nmod_mpoly_set(content, A->coeffs + 0, ctx);
    for (i = 1; i < A->length; i++)
    {
        if (nmod_mpoly_is_one(content, ctx))
            break;
        success = _nmod_mpoly_gcd_zippel(content, content, A->coeffs + i, ctx, 1, randstate);
        if (!success)
            goto finished;
    }
    for (i = 0; i < B->length; i++)
    {
        if (nmod_mpoly_is_one(content, ctx))
            break;
        success = _nmod_mpoly_gcd_zippel(content, content, B->coeffs + i, ctx, 1, randstate);
        if (!success)
            goto finished;
    }

    nmod_mpolyu_divexact_mpoly(Abar, A, content, ctx);
    nmod_mpolyu_divexact_mpoly(Bbar, B, content, ctx);

    ABminshift = FLINT_MIN(Abar->exps[Abar->length - 1], Bbar->exps[Bbar->length - 1]);
    nmod_mpolyu_shift_right(Abar, Abar->exps[Abar->length - 1]);
    nmod_mpolyu_shift_right(Bbar, Bbar->exps[Bbar->length - 1]);

    success = nmod_mpolyu_gcdm_zippel(Gbar, Abar, Bbar, ctx, zinfo, randstate);
    if (!success)
        goto finished;

    nmod_mpolyu_shift_left(Gbar, ABminshift);
    nmod_mpolyu_mul_mpoly(G, Gbar, content, ctx);

    success = 1;

finished:

    nmod_mpolyu_clear(Abar, ctx);
    nmod_mpolyu_clear(Bbar, ctx);
    nmod_mpolyu_clear(Gbar, ctx);
    nmod_mpoly_clear(content, ctx);

    return success;
}



void nmod_mpoly_to_nmod_poly_keepbits(nmod_poly_t A, slong * Ashift,
               const nmod_mpoly_t B, slong var, const nmod_mpoly_ctx_t ctx)
{
    slong i, shift, off, N;
    slong _Ashift = 0, len = B->length;
    mp_limb_t * coeff = B->coeffs;
    ulong * exp = B->exps;
    mp_bitcnt_t bits = B->bits;

    FLINT_ASSERT(bits <= FLINT_BITS);

    N = mpoly_words_per_exp(bits, ctx->minfo);
    mpoly_gen_offset_shift(&off, &shift, var, N, bits, ctx->minfo);

    nmod_poly_zero(A);
    if (len > 0)
    {
        ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
        _Ashift = (exp[N*(len - 1)] >> shift) & mask;
        for (i = 0; i < len; i++)
        {
            ulong k = ((exp[N*i + off] >> shift) & mask) - _Ashift;
            FLINT_ASSERT(((slong)k) >= 0);
            nmod_poly_set_coeff_ui(A, k, coeff[i]);
        }
    }

    *Ashift = _Ashift;
}

void nmod_mpoly_from_nmod_poly_keepbits(nmod_mpoly_t A, const nmod_poly_t B,
                           slong Bshift, slong var, mp_bitcnt_t bits, const nmod_mpoly_ctx_t ctx)
{
    slong shift, off, N;
    slong k;
    slong Alen;
    mp_limb_t * Acoeff;
    ulong * Aexp;
    slong Aalloc;
    ulong * one;
    TMP_INIT;

    TMP_START;

    FLINT_ASSERT(!nmod_poly_is_zero(B));
    FLINT_ASSERT(Bshift >= 0);
    FLINT_ASSERT(Bshift + nmod_poly_degree(B) >= 0);
    FLINT_ASSERT(1 + FLINT_BIT_COUNT(Bshift + nmod_poly_degree(B)) <= bits);

    N = mpoly_words_per_exp(bits, ctx->minfo);
    one = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_gen_oneexp_offset_shift(one, &off, &shift, var, N, bits, ctx->minfo);

    nmod_mpoly_fit_bits(A, bits, ctx);
    A->bits = bits;

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Aalloc = A->alloc;
    Alen = 0;
    for (k = nmod_poly_degree(B); k >= 0; k--)
    {
        _nmod_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, N);
        mpoly_monomial_mul_si(Aexp + N*Alen, one, N, k + Bshift);
        Acoeff[Alen] = nmod_poly_get_coeff_ui(B, k);
        Alen += Acoeff[Alen] != UWORD(0);
    }

    A->coeffs = Acoeff;
    A->exps = Aexp;
    A->alloc = Aalloc;
    _nmod_mpoly_set_length(A, Alen, ctx);

    TMP_END;
}

int _nmod_mpoly_gcd_zippel(nmod_mpoly_t G, const nmod_mpoly_t A,
                              const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx,
                                         int keepbits, flint_rand_t randstate)
{
    slong i;
    int ret, success = 0;
    mpoly_zipinfo_t zinfo;
    nmod_mpoly_ctx_t uctx;
    nmod_mpolyu_t Au, Bu, Gu;
    mp_bitcnt_t new_bits;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT((!keepbits) || A->bits == B->bits);

    FLINT_ASSERT(!nmod_mpoly_is_zero(A, ctx));
    FLINT_ASSERT(!nmod_mpoly_is_zero(B, ctx));

    if (ctx->minfo->nvars == 1) {
        slong shiftA, shiftB;
        nmod_poly_t a, b, g;
        nmod_poly_init(a, ctx->ffinfo->mod.n);
        nmod_poly_init(b, ctx->ffinfo->mod.n);
        nmod_poly_init(g, ctx->ffinfo->mod.n);
        nmod_mpoly_to_nmod_poly_keepbits(a, &shiftA, A, 0, ctx);
        nmod_mpoly_to_nmod_poly_keepbits(b, &shiftB, B, 0, ctx);
        nmod_poly_gcd(g, a, b);
        nmod_mpoly_from_nmod_poly_keepbits(G, g, FLINT_MIN(shiftA, shiftB), 0, A->bits, ctx);
        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(g);
        return 1;
    }

    mpoly_zipinfo_init(zinfo, ctx->minfo->nvars);
    nmod_mpoly_degrees_si(zinfo->Adegs, A, ctx);
    nmod_mpoly_degrees_si(zinfo->Bdegs, B, ctx);
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        zinfo->perm[i] = i;
    }

    new_bits = FLINT_MAX(A->bits, B->bits);

    nmod_mpoly_ctx_init(uctx, ctx->minfo->nvars - 1, ORD_LEX, ctx->ffinfo->mod.n);
    nmod_mpolyu_init(Au, new_bits, uctx);
    nmod_mpolyu_init(Bu, new_bits, uctx);
    nmod_mpolyu_init(Gu, new_bits, uctx);

    nmod_mpoly_to_mpolyu_perm(Au, A, zinfo->perm, uctx, ctx);
    nmod_mpoly_to_mpolyu_perm(Bu, B, zinfo->perm, uctx, ctx);

    ret = nmod_mpolyu_gcd_zippel(Gu, Au, Bu, uctx, zinfo, randstate);
    if (ret) {
        nmod_mpoly_from_mpolyu_perm(G, Gu, keepbits, zinfo->perm, uctx, ctx);
        nmod_mpoly_make_monic(G, G, ctx);
        success = 1;
    }

    nmod_mpolyu_clear(Au, uctx);
    nmod_mpolyu_clear(Bu, uctx);
    nmod_mpolyu_clear(Gu, uctx);
    nmod_mpoly_ctx_clear(uctx);

    mpoly_zipinfo_clear(zinfo);

    return success;
}


int nmod_mpoly_gcd_zippel(nmod_mpoly_t G, const nmod_mpoly_t A,
                              const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
{
    int success;
    flint_rand_t randstate;

    if (nmod_mpoly_is_zero(A, ctx))
    {
        if (nmod_mpoly_is_zero(B, ctx))
        {
            nmod_mpoly_zero(G, ctx);
        } else
        {
            nmod_mpoly_make_monic(G, B, ctx);
        }
        return 1;
    }

    if (nmod_mpoly_is_zero(B, ctx))
    {
        nmod_mpoly_make_monic(G, A, ctx);
        return 1;
    }

    if (A->bits > FLINT_BITS || B->bits > FLINT_BITS)
        return 0;

    if (ctx->minfo->nvars == 1) {
        slong shiftA, shiftB;
        nmod_poly_t a, b, g;
        nmod_poly_init(a, ctx->ffinfo->mod.n);
        nmod_poly_init(b, ctx->ffinfo->mod.n);
        nmod_poly_init(g, ctx->ffinfo->mod.n);
        nmod_mpoly_to_nmod_poly_keepbits(a, &shiftA, A, 0, ctx);
        nmod_mpoly_to_nmod_poly_keepbits(b, &shiftB, B, 0, ctx);
        nmod_poly_gcd(g, a, b);
        nmod_mpoly_from_nmod_poly_keepbits(G, g, FLINT_MIN(shiftA, shiftB), 0, A->bits, ctx);
        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(g);
        return 1;
    }

    flint_randinit(randstate);
    success = _nmod_mpoly_gcd_zippel(G, A, B, ctx, 1, randstate);

    flint_randclear(randstate);

    return success;
}
