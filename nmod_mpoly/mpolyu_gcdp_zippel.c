/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
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
                prod = nmod_mul(prod, powers[j], ctx->mod);
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
                                                             ctx->mod);
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
        NMOD_RED3(v, ac2, ac1, ac0, ctx->mod);

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
        }
        else
        {
            /* otherwise set the coeff of the monomial to one */
            nmod_gcds_ret_t ret;
            FLINT_ASSERT((f->coeffs + 0)->length == 1);
            nmod_mpolyu_set(G, f, ctx);
            (G->coeffs + 0)->coeffs[0] = UWORD(1);
            nmod_mpolyu_init(Aevalsk1, f->bits, ctx);
            ret = nmod_gcds_form_wrong;
            if (   nmod_mpolyuu_divides(Aevalsk1, A, G, 1, ctx)
                && nmod_mpolyuu_divides(Aevalsk1, B, G, 1, ctx))
            {
                ret = nmod_gcds_success;
            }
            nmod_mpolyu_clear(Aevalsk1, ctx);
            return ret;
        }
    }

    TMP_START;

    nmod_mpolyu_init(Aevalsk1, f->bits, ctx);
    nmod_mpolyu_init(Bevalsk1, f->bits, ctx);
    nmod_mpolyu_init(fevalsk1, f->bits, ctx);
    nmod_mpolyu_init(Aevalski, f->bits, ctx);
    nmod_mpolyu_init(Bevalski, f->bits, ctx);
    nmod_mpolyu_init(fevalski, f->bits, ctx);
    nmod_poly_init(Aeval, ctx->mod.n);
    nmod_poly_init(Beval, ctx->mod.n);
    nmod_poly_init(Geval, ctx->mod.n);

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

    nmod_mat_init(MF, 0, l, ctx->mod.n);

    M = (nmod_mat_struct *) TMP_ALLOC(f->length*sizeof(nmod_mat_struct));
    ML_is_initialized = (int *) TMP_ALLOC(f->length*sizeof(int));
    for (i = 0; i < f->length; i++)
    {
        nmod_mat_init(M + i, l, (f->coeffs + i)->length, ctx->mod.n);
        ML_is_initialized[i] = 0;
    }

    W = (mp_limb_t *) flint_malloc(l*f->length*sizeof(mp_limb_t));

    nmod_mat_init(Msol, l, 1, ctx->mod.n);

    /* compute how many masks are needed */
    entries = f->bits * var;
    offs = (slong *) TMP_ALLOC(entries*sizeof(slong));
    masks = (ulong *) TMP_ALLOC(entries*sizeof(slong));
    powers = (mp_limb_t *) TMP_ALLOC(entries*sizeof(mp_limb_t));


    /***** evaluation loop head *******/
    eval_points_tried = 0;
pick_evaluation_point:

    if (++eval_points_tried > 10)
    {
        success = nmod_gcds_eval_point_not_found;
        goto finished;
    }

    /* avoid 0, 1 and -1 for the evaluation points */
    FLINT_ASSERT(ctx->mod.n > UWORD(3));
    for (i = 0; i < var; i++)
        alpha[i] = UWORD(2) + n_randint(randstate, ctx->mod.n - UWORD(3));

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
                powers[f->bits*i + j] = alpha[i];
            else
                powers[f->bits*i + j] = nmod_mul(powers[f->bits*i + j-1],
                                                 powers[f->bits*i + j-1],
                                                             ctx->mod);
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
    nmod_mat_init(MF, 0, l, ctx->mod.n);

    for (S = 0; S < f->length; S++)
    {
        s = d[S];

        if (!ML_is_initialized[s])
        {
            nmod_mat_init(ML + s, l, (f->coeffs + s)->length + l,
                                                           ctx->mod.n);
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
                             l, ctx->mod.n);
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

    nmod_mpolyu_setform(G, f, ctx);

    for (i = 0; i < f->length; i++)
    {
        for (j = 0; j < (f->coeffs + i)->length; j++)
        {
            FLINT_ASSERT((f->coeffs + i)->length <= l);
            b[j] = nmod_mul(W[l*i + j], nmod_mat_get_entry(Msol, j, 0),
                                                             ctx->mod);
        }
        success = nmod_vandsolve((G->coeffs + i)->coeffs,
                                 (fevalsk1->coeffs + i)->coeffs, b, 
                                    (f->coeffs + i)->length, ctx->mod);
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

            NMOD_RED3(v, ac2, ac1, ac0, ctx->mod);
            u = nmod_mul(W[l*s + i], nmod_mat_get_entry(Msol, i, 0),
                                                             ctx->mod);
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



static int nmod_mpolyu_gcdp_zippel_univar(
    nmod_mpolyu_t G,
    nmod_mpolyu_t Abar,
    nmod_mpolyu_t Bbar,
    nmod_mpolyu_t A,
    nmod_mpolyu_t B,
    const nmod_mpoly_ctx_t ctx)
{
    nmod_poly_t a, b, g, t;
    FLINT_ASSERT(A->bits == B->bits);
    nmod_poly_init_mod(a, ctx->mod);
    nmod_poly_init_mod(b, ctx->mod);
    nmod_poly_init_mod(g, ctx->mod);
    nmod_poly_init_mod(t, ctx->mod);
    nmod_mpolyu_cvtto_poly(a, A, ctx);
    nmod_mpolyu_cvtto_poly(b, B, ctx);
    nmod_poly_gcd(g, a, b);
    nmod_mpolyu_cvtfrom_poly(G, g, ctx);
    nmod_poly_div(t, a, g);
    nmod_mpolyu_cvtfrom_poly(Abar, t, ctx);
    nmod_poly_div(t, b, g);
    nmod_mpolyu_cvtfrom_poly(Bbar, t, ctx);
    nmod_poly_clear(a);
    nmod_poly_clear(b);
    nmod_poly_clear(g);
    nmod_poly_clear(t);
    return 1;
}

static int nmod_mpolyu_gcdp_zippel_univar_no_cofactors(
    nmod_mpolyu_t G,
    nmod_mpolyu_t A,
    nmod_mpolyu_t B,
    const nmod_mpoly_ctx_t ctx)
{
    nmod_poly_t a, b, g;
    FLINT_ASSERT(A->bits == B->bits);
    nmod_poly_init_mod(a, ctx->mod);
    nmod_poly_init_mod(b, ctx->mod);
    nmod_poly_init_mod(g, ctx->mod);
    nmod_poly_init_mod(g, ctx->mod);
    nmod_mpolyu_cvtto_poly(a, A, ctx);
    nmod_mpolyu_cvtto_poly(b, B, ctx);
    nmod_poly_gcd(g, a, b);
    nmod_mpolyu_cvtfrom_poly(G, g, ctx);
    nmod_poly_clear(a);
    nmod_poly_clear(b);
    nmod_poly_clear(g);
    return 1;
}


int nmod_mpolyu_gcdp_zippel_bivar(
    nmod_mpolyu_t G,
    nmod_mpolyu_t Abar,
    nmod_mpolyu_t Bbar,
    nmod_mpolyu_t A,
    nmod_mpolyu_t B,
    const nmod_mpoly_ctx_t ctx)
{
    slong var = 0;
    slong Alastdeg;
    slong Blastdeg;
    ulong Ashift, Bshift, Gshift;
    slong lastdeg;
    slong bound;
    int success = 0, changed, have_enough;
    n_poly_t a, b, c, g, modulus, tempmod;
    nmod_mpolyu_t Aeval, Beval, Geval;
    nmod_mpolyun_t An, Bn, H, Ht;
    mp_limb_t geval, temp, alpha;

    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(var >= -WORD(1));

    FLINT_ASSERT(G->bits == A->bits);
    FLINT_ASSERT(G->bits == B->bits);

    nmod_mpolyun_init(An, A->bits, ctx);
    nmod_mpolyun_init(Bn, A->bits, ctx);
    nmod_mpolyu_cvtto_mpolyun(An, A, var, ctx);
    nmod_mpolyu_cvtto_mpolyun(Bn, B, var, ctx);

    FLINT_ASSERT(An->length > 0);
    FLINT_ASSERT(Bn->length > 0);

    Ashift = A->exps[A->length - 1];
    Bshift = B->exps[B->length - 1];
    Gshift = FLINT_MIN(Ashift, Bshift);
    nmod_mpolyun_shift_right(An, Ashift);
    nmod_mpolyun_shift_right(Bn, Bshift);

    n_poly_init(a);
    n_poly_init(b);
    n_poly_init(c);
    n_poly_init(g);

    /* if the gcd has content wrt last variable, we are going to fail */
    nmod_mpolyun_content_last(a, An, ctx);
    nmod_mpolyun_content_last(b, Bn, ctx);
    nmod_mpolyun_divexact_last(An, a, ctx);
    nmod_mpolyun_divexact_last(Bn, b, ctx);
    n_poly_mod_gcd(c, a, b, ctx->mod);
    n_poly_mod_gcd(g, nmod_mpolyun_leadcoeff_poly(An, ctx),
                     nmod_mpolyun_leadcoeff_poly(Bn, ctx), ctx->mod);
    Alastdeg = nmod_mpolyun_lastdeg(An, ctx);
    Blastdeg = nmod_mpolyun_lastdeg(Bn, ctx);

    /* bound of the number of images required */
    bound = 1 + FLINT_MIN(Alastdeg, Blastdeg) + n_poly_degree(g);

    n_poly_init(modulus);
    n_poly_init(tempmod);
    n_poly_set_coeff(tempmod, 1, UWORD(1));
    nmod_mpolyu_init(Aeval, A->bits, ctx);
    nmod_mpolyu_init(Beval, A->bits, ctx);
    nmod_mpolyu_init(Geval, A->bits, ctx);
    nmod_mpolyun_init(H, A->bits, ctx);
    nmod_mpolyun_init(Ht, A->bits, ctx);

    /* fail if the gcd has content wrt last variable */
    if (n_poly_degree(c) > 0)
    {
        success = 0;
        goto finished;
    }

    n_poly_one(modulus);
    nmod_mpolyun_zero(H, ctx);

    alpha = ctx->mod.n;
    while (1)
    {
        if (alpha == 0)
        {
            success = 0;
            goto finished;
        }
        alpha--;

        /* make sure evaluation point does not kill both lc(A) and lc(B) */
        geval = n_poly_mod_evaluate_nmod(g, alpha, ctx->mod);
        if (geval == WORD(0))
            goto outer_continue;

        /* make sure evaluation point does not kill either A or B */
        nmod_mpolyun_interp_reduce_sm_mpolyu(Aeval, An, alpha, ctx);
        nmod_mpolyun_interp_reduce_sm_mpolyu(Beval, Bn, alpha, ctx);
        if (Aeval->length == 0 || Beval->length == 0)
            goto outer_continue;

        nmod_mpolyu_gcdp_zippel_univar_no_cofactors(Geval, Aeval, Beval, ctx);

        if (nmod_mpolyu_is_one(Geval, ctx))
        {
            nmod_mpolyu_one(G, ctx);
            nmod_mpolyu_swap(Abar, A, ctx);
            nmod_mpolyu_swap(Bbar, B, ctx);
            nmod_mpolyu_shift_left(G, Gshift);
            nmod_mpolyu_shift_left(Abar, Ashift - Gshift);
            nmod_mpolyu_shift_left(Bbar, Bshift - Gshift);
            success = 1;
            goto finished;
        }

        FLINT_ASSERT(Geval->length > 0);

        if (n_poly_degree(modulus) > 0)
        {
            if (Geval->exps[0] > H->exps[0])
            {
                goto outer_continue;
            }
            else if (Geval->exps[0] < H->exps[0])
            {
                n_poly_one(modulus);                
            }
        }

        temp = n_invmod(nmod_mpolyu_leadcoeff(Geval, ctx), ctx->mod.n);
        temp = nmod_mul(geval, temp, ctx->mod);
        nmod_mpolyu_scalar_mul_nmod(Geval, temp, ctx);

        /* update interpolant H */
        if (n_poly_degree(modulus) > 0)
        {
            mp_limb_t t = n_poly_mod_evaluate_nmod(modulus, alpha, ctx->mod);
            t = nmod_inv(t, ctx->mod);
            _n_poly_mod_scalar_mul_nmod_inplace(modulus, t, ctx->mod);

            changed = nmod_mpolyun_interp_crt_sm_mpolyu(&lastdeg, H, Ht, Geval,
                                                          modulus, alpha, ctx);
            n_poly_set_coeff(tempmod, 0, ctx->mod.n - alpha);
            n_poly_mod_mul(modulus, modulus, tempmod, ctx->mod);

            have_enough = n_poly_degree(modulus) >= bound;

            if (changed && !have_enough)
            {
                goto outer_continue;
            }

            if (!changed || have_enough)
            {
                nmod_mpolyun_content_last(a, H, ctx);
                nmod_mpolyun_mul_poly(Ht, H, c, ctx);
                nmod_mpolyun_divexact_last(Ht, a, ctx);
                nmod_mpolyun_shift_left(Ht, Gshift);
                nmod_mpolyu_cvtfrom_mpolyun(G, Ht, var, ctx);
                if (   nmod_mpolyuu_divides(Abar, A, G, 1, ctx)
                    && nmod_mpolyuu_divides(Bbar, B, G, 1, ctx))
                {
                    success = 1;
                    goto finished;
                }
            }

            if (have_enough)
            {
                n_poly_one(modulus);
                goto outer_continue;
            }
        }
        else
        {
            nmod_mpolyun_interp_lift_sm_mpolyu(H, Geval, ctx);
            n_poly_set_coeff(tempmod, 0, ctx->mod.n - alpha);
            n_poly_mod_mul(modulus, modulus, tempmod, ctx->mod);
        }

outer_continue:;
    }

finished:

    n_poly_clear(a);
    n_poly_clear(b);
    n_poly_clear(c);
    n_poly_clear(g);
    n_poly_clear(modulus);
    n_poly_clear(tempmod);
    nmod_mpolyu_clear(Aeval, ctx);
    nmod_mpolyu_clear(Beval, ctx);
    nmod_mpolyu_clear(Geval, ctx);
    nmod_mpolyun_clear(An, ctx);
    nmod_mpolyun_clear(Bn, ctx);
    nmod_mpolyun_clear(H, ctx);
    nmod_mpolyun_clear(Ht, ctx);

    return success;
}


int nmod_mpolyu_gcdp_zippel(
    nmod_mpolyu_t G,
    nmod_mpolyu_t Abar,
    nmod_mpolyu_t Bbar,
    nmod_mpolyu_t A,
    nmod_mpolyu_t B,
    slong var,
    const nmod_mpoly_ctx_t ctx,
    flint_rand_t randstate)
{
    slong lastdeg;
    slong Alastdeg;
    slong Blastdeg;
    ulong Ashift, Bshift, Gshift;
    slong degbound;
    slong bound;
    int success = 0, changed, have_enough;
    nmod_mpolyun_t An, Bn;
    n_poly_t a, b, c, g;
    n_poly_t modulus, tempmod;
    nmod_mpolyu_t Aeval, Beval, Geval, Abareval, Bbareval, Gform;
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
        return nmod_mpolyu_gcdp_zippel_univar(G, Abar, Bbar, A, B, ctx);
    }

    if (var == WORD(0))
    {
        /* bivariate is more comfortable separated */
        return nmod_mpolyu_gcdp_zippel_bivar(G, Abar, Bbar, A, B, ctx);
    }

    nmod_mpolyun_init(An, A->bits, ctx);
    nmod_mpolyun_init(Bn, A->bits, ctx);
    nmod_mpolyu_cvtto_mpolyun(An, A, var, ctx);
    nmod_mpolyu_cvtto_mpolyun(Bn, B, var, ctx);

    Ashift = A->exps[A->length - 1];
    Bshift = B->exps[B->length - 1];
    Gshift = FLINT_MIN(Ashift, Bshift);
    nmod_mpolyun_shift_right(An, Ashift);
    nmod_mpolyun_shift_right(Bn, Bshift);

    n_poly_init(a);
    n_poly_init(b);
    n_poly_init(c);
    n_poly_init(g);

    /* if the gcd has content wrt last variable, we are going to fail */
    nmod_mpolyun_content_last(a, An, ctx);
    nmod_mpolyun_content_last(b, Bn, ctx);
    nmod_mpolyun_divexact_last(An, a, ctx);
    nmod_mpolyun_divexact_last(Bn, b, ctx);
    n_poly_mod_gcd(c, a, b, ctx->mod);
    n_poly_mod_gcd(g, nmod_mpolyun_leadcoeff_poly(An, ctx),
                      nmod_mpolyun_leadcoeff_poly(Bn, ctx), ctx->mod);
    Alastdeg = nmod_mpolyun_lastdeg(An, ctx);
    Blastdeg = nmod_mpolyun_lastdeg(Bn, ctx);

    /* bound of the number of images required */
    bound = 1 + FLINT_MIN(Alastdeg, Blastdeg) + n_poly_degree(g);

    /* degree bound on the gcd */
    degbound = FLINT_MIN(A->exps[0], B->exps[0]);

    n_poly_init(modulus);
    n_poly_init(tempmod);
    n_poly_set_coeff(tempmod, 1, UWORD(1));
    nmod_mpolyu_init(Aeval, A->bits, ctx);
    nmod_mpolyu_init(Beval, A->bits, ctx);
    nmod_mpolyu_init(Geval, A->bits, ctx);
    nmod_mpolyu_init(Abareval, A->bits, ctx);
    nmod_mpolyu_init(Bbareval, A->bits, ctx);
    nmod_mpolyu_init(Gform, A->bits, ctx);
    nmod_mpolyun_init(H, A->bits, ctx);
    nmod_mpolyun_init(Ht, A->bits, ctx);

    /* fail if the gcd has content wrt last variable */
    if (n_poly_degree(c) > 0)
    {
        success = 0;
        goto finished;
    }

    if (ctx->mod.n <= UWORD(3))
    {
        success = 0;
        goto finished;
    }

    n_poly_one(modulus);
    nmod_mpolyun_zero(H, ctx);

    start_alpha = UWORD(1) + n_randint(randstate, ctx->mod.n - UWORD(1));
    alpha = start_alpha;
    while (1)
    {
        /* get new evaluation point */
        --alpha;
        if (alpha == 0)
        {
            alpha = ctx->mod.n - UWORD(1);
        }
        if (alpha == start_alpha)
        {
            success = 0;
            goto finished;
        }

        /* make sure evaluation point does not kill both lc(A) and lc(B) */
        geval = n_poly_mod_evaluate_nmod(g, alpha, ctx->mod);
        if (geval == 0)
        {
            goto outer_continue;
        }

        /* make sure evaluation point does not kill either A or B */
        nmod_mpolyun_interp_reduce_sm_mpolyu(Aeval, An, alpha, ctx);
        nmod_mpolyun_interp_reduce_sm_mpolyu(Beval, Bn, alpha, ctx);
        if (Aeval->length == 0 || Beval->length == 0)
        {
            goto outer_continue;
        }

        success = nmod_mpolyu_gcdp_zippel(Geval, Abareval, Bbareval, Aeval, Beval,
                                                      var - 1, ctx, randstate);
        if (!success || Geval->exps[0] > degbound)
        {
            success = 0;
            goto finished;
        }
        
        degbound = Geval->exps[0];

        if (nmod_mpolyu_is_one(Geval, ctx))
        {
            nmod_mpolyu_one(G, ctx);
            nmod_mpolyu_swap(Abar, A, ctx);
            nmod_mpolyu_swap(Bbar, B, ctx);
            nmod_mpolyu_shift_left(G, Gshift);
            nmod_mpolyu_shift_left(Abar, Ashift - Gshift);
            nmod_mpolyu_shift_left(Bbar, Bshift - Gshift);
            success = 1;
            goto finished;
        }

        if (n_poly_degree(modulus) > 0)
        {
            if (Geval->exps[0] > H->exps[0])
            {
                goto outer_continue;
            }
            else if (Geval->exps[0] < H->exps[0])
            {
                n_poly_one(modulus);                
            }
        }

        /* update interpolant H */
        temp = nmod_inv(nmod_mpolyu_leadcoeff(Geval, ctx), ctx->mod);
        temp = nmod_mul(geval, temp, ctx->mod);
        nmod_mpolyu_scalar_mul_nmod(Geval, temp, ctx);
        if (n_poly_degree(modulus) > 0)
        {
            temp = n_poly_mod_evaluate_nmod(modulus, alpha, ctx->mod);
            temp = n_invmod(temp, ctx->mod.n);
            _n_poly_mod_scalar_mul_nmod_inplace(modulus, temp, ctx->mod);
            changed = nmod_mpolyun_interp_crt_sm_mpolyu(&lastdeg, H, Ht, Geval,
                                                          modulus, alpha, ctx);
            if (!changed)
            {
                nmod_mpolyun_content_last(a, H, ctx);
                nmod_mpolyun_mul_poly(Ht, H, c, ctx);
                nmod_mpolyun_divexact_last(Ht, a, ctx);
                nmod_mpolyun_shift_left(Ht, Gshift);
                nmod_mpolyu_cvtfrom_mpolyun(G, Ht, var, ctx);
                if (   !nmod_mpolyuu_divides(Abar, A, G, 1, ctx)
                    || !nmod_mpolyuu_divides(Bbar, B, G, 1, ctx))
                {
                    goto outer_continue;
                }
                success = 1;
                goto finished;
            }
        }
        else
        {
            nmod_mpolyun_interp_lift_sm_mpolyu(H, Geval, ctx);
        }
        n_poly_set_coeff(tempmod, 0, ctx->mod.n - alpha);
        n_poly_mod_mul(modulus, modulus, tempmod, ctx->mod);

        nmod_mpolyu_setform_mpolyun(Gform, H, ctx);

        while (1)
        {
            /* get new evaluation point */
            --alpha;
            if (alpha == 0)
            {
                alpha = ctx->mod.n - UWORD(1);
            }
            if (alpha == start_alpha)
            {
                success = 0;
                goto finished;
            }

            /* make sure evaluation does not kill both lc(A) and lc(B) */
            geval = n_poly_mod_evaluate_nmod(g, alpha, ctx->mod);
            if (geval == WORD(0))
            {
                goto inner_continue;
            }

            /* make sure evaluation does not kill either A or B */
            nmod_mpolyun_interp_reduce_sm_mpolyu(Aeval, An, alpha, ctx);
            nmod_mpolyun_interp_reduce_sm_mpolyu(Beval, Bn, alpha, ctx);
            if (Aeval->length == 0 || Beval->length == 0)
            {
                goto inner_continue;
            }

            switch (nmod_mpolyu_gcds_zippel(Geval, Aeval, Beval, Gform, var,
                                                    ctx, randstate, &degbound))
            {
                default:
                    FLINT_ASSERT(0);
                case nmod_gcds_form_main_degree_too_high:
                    /* nmod_mpolyu_gcds_zippel has updated degbound */
                    n_poly_one(modulus);
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
            {
                goto inner_continue;
            }

            /* update interpolant H */
            temp = nmod_inv(nmod_mpolyu_leadcoeff(Geval, ctx), ctx->mod);
            nmod_mpolyu_scalar_mul_nmod(Geval, nmod_mul(geval, temp,
                                                       ctx->mod), ctx);
            FLINT_ASSERT(n_poly_degree(modulus) > 0);
            temp = n_poly_mod_evaluate_nmod(modulus, alpha, ctx->mod);
            temp = nmod_inv(temp, ctx->mod);
            _n_poly_mod_scalar_mul_nmod_inplace(modulus, temp, ctx->mod);

            changed = nmod_mpolyun_interp_crt_sm_mpolyu(&lastdeg, H, Ht, Geval,
                                                          modulus, alpha, ctx);
            n_poly_set_coeff(tempmod, 0, ctx->mod.n - alpha);
            n_poly_mod_mul(modulus, modulus, tempmod, ctx->mod);

            have_enough = n_poly_degree(modulus) >= bound;

            if (changed && !have_enough)
            {
                goto inner_continue;
            }

            if (!changed || have_enough)
            {
                nmod_mpolyun_content_last(a, H, ctx);
                nmod_mpolyun_mul_poly(Ht, H, c, ctx);
                nmod_mpolyun_divexact_last(Ht, a, ctx);
                nmod_mpolyun_shift_left(Ht, Gshift);
                nmod_mpolyu_cvtfrom_mpolyun(G, Ht, var, ctx);
                if (    nmod_mpolyuu_divides(Abar, A, G, 1, ctx)
                     && nmod_mpolyuu_divides(Bbar, B, G, 1, ctx))
                {
                    success = 1;
                    goto finished;
                }
            }

            if (have_enough)
            {
                n_poly_one(modulus);
                goto outer_continue;
            }

inner_continue:;
        }
        FLINT_ASSERT(0 && "not reachable");

outer_continue:;
    }
    FLINT_ASSERT(0 && "not reachable");

finished:

    n_poly_clear(a);
    n_poly_clear(b);
    n_poly_clear(c);
    n_poly_clear(g);
    n_poly_clear(modulus);
    n_poly_clear(tempmod);
    nmod_mpolyu_clear(Aeval, ctx);
    nmod_mpolyu_clear(Beval, ctx);
    nmod_mpolyu_clear(Geval, ctx);
    nmod_mpolyu_clear(Abareval, ctx);
    nmod_mpolyu_clear(Bbareval, ctx);
    nmod_mpolyu_clear(Gform, ctx);
    nmod_mpolyun_clear(An, ctx);
    nmod_mpolyun_clear(Bn, ctx);
    nmod_mpolyun_clear(H, ctx);
    nmod_mpolyun_clear(Ht, ctx);

    return success;
}
