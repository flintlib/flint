/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

/* todo: move out? */
static int
fmpz_mat_is_diagonal(const fmpz_mat_t A)
{
    slong r = fmpz_mat_nrows(A);
    slong c = fmpz_mat_ncols(A);
    slong j, k;

    if (r != c)
        return 0;

    for (j = 0; j < r; j++)
        for (k = 0; k < c; k++)
            if (j != k && !fmpz_is_zero(fmpz_mat_entry(A, j, k)))
                return 0;

    return 1;
}

/* todo: move out? Adapt fmpz_mat_snf algorithms to also output transform? */
/* Compute Smith normal form of A and invertible U, V s.t. S = UAV assuming A
   is square */
static void
fmpz_mat_snf_transform(fmpz_mat_t S, fmpz_mat_t U, fmpz_mat_t V, const fmpz_mat_t A)
{
    slong g = fmpz_mat_nrows(A);
    fmpz_mat_t X, M;
    fmpz_t d, u, v, p, q;
    slong j, k;

    fmpz_mat_init(X, g, g);
    fmpz_mat_init(M, g, g);
    fmpz_init(d);
    fmpz_init(u);
    fmpz_init(v);
    fmpz_init(p);
    fmpz_init(q);

    fmpz_mat_set(X, A);
    fmpz_mat_one(U);
    fmpz_mat_one(V);

    while (!fmpz_mat_is_diagonal(X))
    {
        fmpz_mat_hnf_transform(X, M, X);
        fmpz_mat_mul(U, M, U);
        fmpz_mat_transpose(X, X);
        fmpz_mat_hnf_transform(X, M, X);
        fmpz_mat_transpose(X, X);
        fmpz_mat_transpose(M, M);
        fmpz_mat_mul(V, V, M);
    }

    for (j = 0; j < g; j++)
    {
        if (fmpz_is_one(fmpz_mat_entry(X, j, j)))
        {
            continue;
        }
        for (k = j + 1; k < g; k++)
        {
            if (fmpz_is_zero(fmpz_mat_entry(X, k, k)))
            {
                continue;
            }
            fmpz_xgcd_canonical_bezout(d, u, v,
                fmpz_mat_entry(X, j, j), fmpz_mat_entry(X, k, k));
            fmpz_divexact(p, fmpz_mat_entry(X, j, j), d);
            fmpz_divexact(q, fmpz_mat_entry(X, k, k), d);

            fmpz_mat_one(M);
            fmpz_set(fmpz_mat_entry(M, j, k), v);
            fmpz_set(fmpz_mat_entry(M, k, j), q);
            fmpz_mul(fmpz_mat_entry(M, k, k), v, q);
            fmpz_add_si(fmpz_mat_entry(M, k, k), fmpz_mat_entry(M, k, k), -1);
            fmpz_mat_mul(U, M, U);
            fmpz_mat_mul(X, M, X);

            fmpz_mat_one(M);
            fmpz_set(fmpz_mat_entry(M, j, j), u);
            fmpz_one(fmpz_mat_entry(M, k, j));
            fmpz_neg(fmpz_mat_entry(M, k, k), p);
            fmpz_mul(fmpz_mat_entry(M, j, k), fmpz_mat_entry(M, k, k), u);
            fmpz_add_si(fmpz_mat_entry(M, j, k), fmpz_mat_entry(M, j, k), 1);
            fmpz_mat_mul(V, V, M);
            fmpz_mat_mul(X, X, M);
        }
    }

    fmpz_mat_set(S, X);

    fmpz_mat_clear(X);
    fmpz_mat_clear(M);
    fmpz_clear(d);
    fmpz_clear(u);
    fmpz_clear(v);
    fmpz_clear(p);
    fmpz_clear(q);
}

static fmpz_mat_struct*
sp2gz_decompose_g1(slong* nb, const fmpz_mat_t mat)
{
    fmpz_mat_struct* res;

    res = flint_malloc(1 * sizeof(fmpz_mat_struct));
    fmpz_mat_init(res, 2, 2);
    fmpz_mat_set(res, mat);
    *nb = 1;
    return res;
}

static fmpz_mat_struct*
sp2gz_decompose_nonsimplified(slong* nb, const fmpz_mat_t mat)
{
    slong g = sp2gz_dim(mat);
    fmpz_mat_t gamma, delta, last;
    fmpz_mat_t u, v, d;
    fmpz_mat_t cur, left, right, m;
    fmpz_mat_t w;
    fmpz_mat_struct* vec;
    fmpz_mat_struct* rec = NULL;
    fmpz_mat_struct* res;
    fmpz_t a;
    slong nb_rec = 0;
    slong nb_max;
    slong nb_vec = 0;
    slong r, k, j;


    if (g == 1)
    {
        return sp2gz_decompose_g1(nb, mat);
    }

    fmpz_mat_init(u, g, g);
    fmpz_mat_init(v, g, g);
    fmpz_mat_init(d, g, g);
    fmpz_mat_init(cur, 2 * g, 2 * g);
    fmpz_mat_init(left, 2 * g, 2 * g);
    fmpz_mat_init(right, 2 * g, 2 * g);
    fmpz_mat_init(m, 2 * g, 2 * g);
    fmpz_init(a);

    fmpz_mat_set(cur, mat);
    fmpz_mat_window_init(gamma, cur, g, 0, 2 * g, g);
    fmpz_mat_snf_transform(d, u, v, gamma);
    fmpz_mat_window_clear(gamma);

    r = fmpz_mat_rank(d);
    fmpz_mat_transpose(u, u);
    sp2gz_block_diag(left, u);
    sp2gz_inv(left, left);
    sp2gz_block_diag(right, v);
    fmpz_mat_mul(cur, left, cur);
    fmpz_mat_mul(cur, cur, right);

    nb_max = 3 * fmpz_bits(fmpz_mat_entry(d, 0, 0)) + 4;
    vec = flint_malloc(nb_max * sizeof(fmpz_mat_struct));
    for (k = 0; k < nb_max; k++)
    {
        fmpz_mat_init(&vec[k], 2 * g, 2 * g);
    }

    fmpz_mat_set(&vec[nb_vec], right);
    nb_vec++;

    while (r == g)
    {
        /* Set u such that delta + gamma*u is reduced, update vec */
        fmpz_mat_zero(u);
        for (j = 0; j < g; j++)
        {
            for (k = j; k < g; k++)
            {
                fmpz_smod(a, fmpz_mat_entry(cur, g + j, g + k), fmpz_mat_entry(cur, g + j, j));
                fmpz_sub(a, a, fmpz_mat_entry(cur, g + j, g + k));
                fmpz_divexact(fmpz_mat_entry(u, j, k), a, fmpz_mat_entry(cur, g + j, j));
                fmpz_set(fmpz_mat_entry(u, k, j), fmpz_mat_entry(u, j, k));
            }
        }
        sp2gz_trig(right, u);
        fmpz_mat_set(&vec[nb_vec], right);
        fmpz_mat_mul(cur, cur, right);
        nb_vec++;

        /* Swap c, d */
        sp2gz_j(&vec[nb_vec]);
        sp2gz_inv(&vec[nb_vec], &vec[nb_vec]);
        fmpz_mat_mul(cur, cur, &vec[nb_vec]);
        nb_vec++;

        /* Recompute SNF */
        fmpz_mat_window_init(gamma, cur, g, 0, 2 * g, g);
        fmpz_mat_snf_transform(d, u, v, gamma);
        fmpz_mat_window_clear(gamma);

        r = fmpz_mat_rank(d);
        fmpz_mat_transpose(u, u);
        sp2gz_block_diag(m, u);
        sp2gz_inv(m, m);
        fmpz_mat_mul(left, m, left);
        fmpz_mat_mul(cur, m, cur);

        sp2gz_block_diag(&vec[nb_vec], v);
        fmpz_mat_mul(cur, cur, &vec[nb_vec]);
        nb_vec++;
    }

    /* Now r < g: make HNF on colums for the bottom of delta and recursive call */
    fmpz_mat_init(last, g, g - r);
    for (k = 0; k < g - r; k++)
    {
        for (j = 0; j < g; j++)
        {
            fmpz_set(fmpz_mat_entry(last, j, k), fmpz_mat_entry(cur, g + r + k, g + j));
        }
    }
    fmpz_mat_hnf_transform(last, u, last);
    for (j = 0; j < g - r; j++)
    {
        fmpz_mat_swap_rows(u, NULL, g - 1 - j, g - 1 - j - r);
    }

    sp2gz_block_diag(&vec[nb_vec], u);
    sp2gz_inv(&vec[nb_vec], &vec[nb_vec]);
    fmpz_mat_mul(cur, cur, &vec[nb_vec]);
    nb_vec++;

    if (r > 0)
    {
        fmpz_mat_init(w, 2 * r, 2 * r);
        sp2gz_restrict(w, cur);
        rec = sp2gz_decompose(&nb_rec, w);

        sp2gz_embed(right, w);
        sp2gz_inv(right, right);
        fmpz_mat_mul(cur, right, cur);

        fmpz_mat_window_init(delta, cur, g, g, 2 * g, 2 * g);

        sp2gz_block_diag(&vec[nb_vec], delta);
        fmpz_mat_transpose(&vec[nb_vec], &vec[nb_vec]);
        fmpz_mat_mul(cur, cur, &vec[nb_vec]);
        nb_vec++;

        fmpz_mat_window_clear(delta);
        fmpz_mat_clear(w);
    }
    sp2gz_inv(&vec[nb_vec], cur);
    nb_vec++;

    /* Make final vector */
    *nb = 1 + nb_rec + nb_vec;
    res = flint_malloc(*nb * sizeof(fmpz_mat_struct));
    for (k = 0; k < *nb; k++)
    {
        fmpz_mat_init(&res[k], 2 * g, 2 * g);
    }

    sp2gz_inv(&res[0], left);
    for (k = 0; k < nb_rec; k++)
    {
        sp2gz_embed(&res[1 + k], &rec[k]);
    }
    for (k = 0; k < nb_vec; k++)
    {
        sp2gz_inv(&res[1 + nb_rec + k], &vec[nb_vec - 1 - k]);
    }

    fmpz_mat_clear(u);
    fmpz_mat_clear(v);
    fmpz_mat_clear(d);
    fmpz_mat_clear(cur);
    fmpz_mat_clear(left);
    fmpz_mat_clear(right);
    fmpz_mat_clear(m);
    fmpz_mat_clear(last);
    fmpz_clear(a);
    for (k = 0; k < nb_max; k++)
    {
        fmpz_mat_clear(&vec[k]);
    }
    flint_free(vec);
    if (r > 0)
    {
        for (k = 0; k < nb_rec; k++)
        {
            fmpz_mat_clear(&rec[k]);
        }
        flint_free(rec);
    }
    return res;
}

fmpz_mat_struct* sp2gz_decompose(slong* nb, const fmpz_mat_t mat)
{
    slong g = sp2gz_dim(mat);
    fmpz_mat_struct* rec;
    slong nb_rec;
    fmpz_mat_struct* res;
    fmpz_mat_t u, beta, delta;
    slong k, next_k, j;

    fmpz_mat_init(u, g, g);
    rec = sp2gz_decompose_nonsimplified(&nb_rec, mat);

    /* Move block-diagonal matrices to the left of rec as much as possible */
    k = 0;
    while (k < nb_rec)
    {
        for (j = k; j < nb_rec; j++)
            if (!sp2gz_is_block_diag(&rec[j])
                && !sp2gz_is_trig(&rec[j])
                && !sp2gz_is_j(&rec[j]))
                break;
        next_k = j + 1;

        /* Move all block-diag matrices between k and next_k to the left */
        for (j = next_k - 2; j >= k; j--)
            if (sp2gz_is_block_diag(&rec[j]))
                break;

        for (; j >= k + 1; j--)
        {
            /* Commutation of rec[j-1] and rec[j] */
            if (sp2gz_is_block_diag(&rec[j - 1]))
            {
                fmpz_mat_mul(&rec[j - 1], &rec[j - 1], &rec[j]);
                fmpz_mat_one(&rec[j]);
            }
            else if (sp2gz_is_trig(&rec[j - 1]))
            {
                fmpz_mat_window_init(beta, &rec[j - 1], 0, g, g, 2 * g);
                fmpz_mat_window_init(delta, &rec[j], g, g, 2 * g, 2 * g);
                fmpz_mat_transpose(u, delta);
                fmpz_mat_mul(u, u, beta);
                fmpz_mat_mul(u, u, delta);

                fmpz_mat_set(&rec[j - 1], &rec[j]);
                sp2gz_trig(&rec[j], u);
                fmpz_mat_window_clear(beta);
                fmpz_mat_window_clear(delta);
            }
            else if (sp2gz_is_j(&rec[j - 1]))
            {
                sp2gz_inv(&rec[j - 1], &rec[j]);
                fmpz_mat_transpose(&rec[j - 1], &rec[j - 1]);
                sp2gz_j(&rec[j]);
            }
        }
        k = next_k;
    }

    /* Move trigonal matrices to the left of rec as much as possible */
    for (k = nb_rec - 1; k >= 1; k--)
    {
        if (sp2gz_is_trig(&rec[k]) && sp2gz_is_trig(&rec[k - 1]))
        {
            fmpz_mat_mul(&rec[k - 1], &rec[k - 1], &rec[k]);
            fmpz_mat_one(&rec[k]);
        }
    }

    *nb = 0;
    for (k = 0; k < nb_rec; k++)
    {
        if (!fmpz_mat_is_one(&rec[k]))
        {
            (*nb)++;
        }
    }
    res = flint_malloc(*nb * sizeof(fmpz_mat_struct));
    for (k = 0; k < *nb; k++)
    {
        fmpz_mat_init(&res[k], 2 * g, 2 * g);
    }

    k = 0;
    for (j = 0; j < nb_rec; j++)
    {
        if (!fmpz_mat_is_one(&rec[j]))
        {
            fmpz_mat_set(&res[k], &rec[j]);
            k++;
        }
    }

    fmpz_mat_clear(u);
    for (k = 0; k < nb_rec; k++)
    {
        fmpz_mat_clear(&rec[k]);
    }
    flint_free(rec);
    return res;
}


/* static fmpz_mat_struct* */
/* sp2gz_decompose_g1_fd(slong* nb, const fmpz_mat_t mat) */
/* { */
/*     acb_mat_t tau0, tau; */
/*     arf_t tol; */
/*     fmpz_t x; */
/*     fmpz_mat_t test, m; */
/*     fmpz_mat_struct* res; */
/*     slong prec; */
/*     slong j, k; */

/*     acb_mat_init(tau, 1, 1); */
/*     acb_mat_init(tau0, 1, 1); */
/*     fmpz_mat_init(test, 2, 2); */
/*     fmpz_mat_init(m, 2, 2); */
/*     arf_init(tol); */
/*     fmpz_init(x); */

/*     flint_printf("(decompose_g1) got:\n"); */
/*     fmpz_mat_print_pretty(mat); */
/*     flint_printf("\n"); */

/*     prec = 0; */
/*     for (j = 0; j < 2; j++) */
/*     { */
/*         for (k = 0; k < 2; k++) */
/*         { */
/*             prec = FLINT_MAX(prec, fmpz_bits(fmpz_mat_entry(mat, j, k))); */
/*         } */
/*     } */
/*     prec += ACB_THETA_LOW_PREC; */

/*     acb_mat_onei(tau0); */
/*     acb_mat_scalar_mul_2exp_si(tau0, tau0, 1); */
/*     acb_siegel_transform(tau0, mat, tau0, prec); */
/*     arf_one(tol); */
/*     arf_mul_2exp_si(tol, tol, -10); */
/*     acb_mat_set(tau, tau0); */

/*     *nb = 0; */
/*     fmpz_mat_set(test, mat); */
/*     while (!acb_modular_is_in_fundamental_domain(acb_mat_entry(tau, 0, 0), tol, prec)) */
/*     { */
/*         arf_get_fmpz(x, arb_midref(acb_realref(acb_mat_entry(tau, 0, 0))), */
/*             ARF_RND_NEAR); */
/*         arb_sub_fmpz(acb_realref(acb_mat_entry(tau, 0, 0)), */
/*             acb_realref(acb_mat_entry(tau, 0, 0)), x, prec); */

/*         fmpz_mat_one(m); */
/*         fmpz_neg(fmpz_mat_entry(m, 0, 1), x); */
/*         fmpz_mat_mul(test, m, test); */
/*         (*nb)++; */

/*         if (!acb_modular_is_in_fundamental_domain(acb_mat_entry(tau, 0, 0), tol, prec)) */
/*         { */
/*             acb_mat_inv(tau, tau, prec); */
/*             acb_mat_neg(tau, tau); */

/*             sp2gz_j(m); */
/*             fmpz_mat_neg(m, m); */
/*             fmpz_mat_mul(test, m, test); */
/*             (*nb)++; */
/*         } */
/*     } */

/*     flint_printf("(decompose_g1) found nb = %wd, test:\n", *nb); */
/*     fmpz_mat_print_pretty(test); */
/*     flint_printf("\n"); */

/*     if (!fmpz_mat_is_one(test)) */
/*     { */
/*         (*nb)++; */
/*     } */

/*     res = flint_malloc(*nb * sizeof(fmpz_mat_struct)); */
/*     for (k = 0; k < *nb; k++) */
/*     { */
/*         fmpz_mat_init(&res[k], 2, 2); */
/*     } */

/*     acb_mat_set(tau, tau0); */

/*     k = 0; */
/*     if (!fmpz_mat_is_one(test)) */
/*     { */
/*         fmpz_mat_one(&res[k]); */
/*         fmpz_mat_neg(&res[k], &res[k]); */
/*         k++; */
/*     } */
/*     while (!acb_modular_is_in_fundamental_domain(acb_mat_entry(tau, 0, 0), tol, prec)) */
/*     { */
/*         arf_get_fmpz(x, arb_midref(acb_realref(acb_mat_entry(tau, 0, 0))), */
/*             ARF_RND_NEAR); */
/*         arb_sub_fmpz(acb_realref(acb_mat_entry(tau, 0, 0)), */
/*             acb_realref(acb_mat_entry(tau, 0, 0)), x, prec); */

/*         fmpz_mat_one(&res[k]); */
/*         fmpz_set(fmpz_mat_entry(&res[k], 0, 1), x); */
/*         k++; */

/*         if (!acb_modular_is_in_fundamental_domain(acb_mat_entry(tau, 0, 0), tol, prec)) */
/*         { */
/*             acb_mat_inv(tau, tau, prec); */
/*             acb_mat_neg(tau, tau); */
/*             sp2gz_j(&res[k]); */
/*             k++; */
/*         } */
/*     } */

/*     acb_mat_clear(tau); */
/*     acb_mat_clear(tau0); */
/*     arf_clear(tol); */
/*     fmpz_clear(x); */
/*     return res; */
/* } */

/* static void */
/* sp2gz_reduce_delta(fmpz_mat_t res, fmpz_mat_t w, const fmpz_mat_t mat) */
/* { */
/*     slong g = sp2gz_dim(mat); */
/*     fmpz_mat_t delta, diag, u, v; */
/*     slong k; */

/*     fmpz_mat_init(u, g, g); */
/*     fmpz_mat_init(v, g, g); */
/*     fmpz_mat_init(diag, g, g); */
/*     fmpz_mat_window_init(delta, mat, g, g, 2 * g, 2 * g); */

/*     for (k = 0; k < g; k++) */
/*     { */
/*         fmpz_one(fmpz_mat_entry(diag, k, g - 1 - k)); */
/*     } */
/*     fmpz_mat_transpose(v, delta); */
/*     fmpz_mat_mul(v, v, diag); */
/*     fmpz_mat_hnf_transform(v, u, v); */
/*     fmpz_mat_mul(u, diag, u); */

/*     sp2gz_block_diag(w, u); */
/*     sp2gz_inv(w, w); */
/*     fmpz_mat_mul(res, mat, w); */

/*     fmpz_mat_clear(u); */
/*     fmpz_mat_clear(v); */
/*     fmpz_mat_clear(diag); */
/*     fmpz_mat_window_clear(delta); */
/* } */

/* static slong */
/* sp2gz_reduce_gamma(fmpz_mat_t next, fmpz_mat_struct* vec, const fmpz_mat_t mat) */
/* { */
/*     slong g = sp2gz_dim(mat); */
/*     fmpz_mat_t cur, u; */
/*     fmpz_t d, r; */
/*     slong k; */
/*     slong res = 0; */
/*     int test = 0; */

/*     /\* Return 0 if gamma is already zero *\/ */
/*     for (k = 0; (k < g) && !test; k++) */
/*     { */
/*         if (!fmpz_is_zero(fmpz_mat_entry(mat, 2 * g - 1, k))) */
/*         { */
/*             test = 1; */
/*         } */
/*     } */
/*     if (!test) */
/*     { */
/*         fmpz_mat_set(next, mat); */
/*         return res; */
/*     } */

/*     fmpz_mat_init(cur, 2 * g, 2 * g); */
/*     fmpz_mat_init(u, g, g); */
/*     fmpz_init(d); */
/*     fmpz_init(r); */

/*     /\* Reduce last line of gamma *\/ */
/*     sp2gz_j(cur); */
/*     fmpz_mat_mul(cur, mat, cur); */
/*     sp2gz_reduce_delta(cur, &vec[res], cur); */
/*     fmpz_mat_transpose(&vec[res], &vec[res]); */
/*     sp2gz_inv(&vec[res], &vec[res]); */
/*     fmpz_mat_mul(cur, mat, &vec[res]); */
/*     if (!fmpz_mat_is_one(&vec[res])) */
/*     { */
/*         res++; */
/*     } */

/*     /\* Reduce last line of delta mod d *\/ */
/*     fmpz_set(d, fmpz_mat_entry(cur, 2 * g - 1, g - 1)); */
/*     for (k = 0; k < g; k++) */
/*     { */
/*         fmpz_smod(r, fmpz_mat_entry(cur, 2 * g - 1, g + k), d); */
/*         fmpz_sub(r, r, fmpz_mat_entry(cur, 2 * g - 1, g + k)); */
/*         fmpz_divexact(fmpz_mat_entry(u, g - 1, k), r, d); */
/*         fmpz_set(fmpz_mat_entry(u, k, g - 1), fmpz_mat_entry(u, g - 1, k)); */
/*     } */

/*     /\* Todo: faster reduction using other entries of u as well? *\/ */
/*     sp2gz_trig(&vec[res], u); */
/*     fmpz_mat_mul(cur, cur, &vec[res]); */
/*     if (!fmpz_mat_is_one(&vec[res])) */
/*     { */
/*         res++; */
/*     } */

/*     /\* Exchange c and d *\/ */
/*     sp2gz_j(&vec[res]); */
/*     sp2gz_inv(&vec[res], &vec[res]); */
/*     fmpz_mat_mul(next, cur, &vec[res]); */
/*     res++; */

/*     fmpz_mat_clear(cur); */
/*     fmpz_mat_clear(u); */
/*     fmpz_clear(d); */
/*     fmpz_clear(r); */
/*     return res; */
/* } */

/* static slong */
/* sp2gz_get_parabolic(fmpz_mat_t next, fmpz_mat_struct* vec, const fmpz_mat_t mat) */
/* { */
/*     slong g = sp2gz_dim(mat); */
/*     slong res = 0; */
/*     fmpz_mat_t u, alpha; */

/*     fmpz_mat_init(u, 2 * g, 2 * g); */

/*     sp2gz_restrict(next, mat); */
/*     sp2gz_embed(u, next); */
/*     sp2gz_inv(u, u); */
/*     fmpz_mat_mul(u, mat, u); */

/*     fmpz_mat_window_init(alpha, u, 0, 0, g, g); */
/*     if (!fmpz_mat_is_one(alpha)) */
/*     { */
/*         sp2gz_block_diag(&vec[res], alpha); */
/*         sp2gz_inv(&vec[res], &vec[res]); */
/*         fmpz_mat_mul(u, u, &vec[res]); */
/*         res++; */
/*     } */
/*     fmpz_mat_window_clear(alpha); */

/*     fmpz_mat_window_init(alpha, u, 0, g, g, 2 * g); */
/*     if (!fmpz_mat_is_zero(alpha)) */
/*     { */
/*         sp2gz_trig(&vec[res], alpha); */
/*         sp2gz_inv(&vec[res], &vec[res]); */
/*         res++; */
/*     } */
/*     fmpz_mat_window_clear(alpha); */

/*     fmpz_mat_clear(u); */
/*     return res; */
/* } */

/* fmpz_mat_struct* sp2gz_decompose(slong* nb, const fmpz_mat_t mat) */
/* { */
/*     slong g = sp2gz_dim(mat); */
/*     fmpz_mat_t cur; */
/*     fmpz_mat_struct* gamma; */
/*     fmpz_mat_struct* rec = NULL; */
/*     fmpz_mat_struct* res; */
/*     fmpz_mat_t w; */
/*     slong nb_max = 0; */
/*     slong nb_rec = 0; */
/*     slong k, nb_gamma, add; */

/*     fmpz_mat_init(cur, 2 * g, 2 * g); */

/*     /\* We need at most 5 * bits matrices to reduce gamma to zero *\/ */
/*     for (k = 0; k < g; k++) */
/*     { */
/*         nb_max = FLINT_MAX(nb_max, fmpz_bits(fmpz_mat_entry(mat, 2 * g - 1, k))); */
/*     } */
/*     nb_max = 3 * nb_max + 3; /\* for last delta reduction *\/ */

/*     gamma = flint_malloc(nb_max * sizeof(fmpz_mat_struct)); */
/*     for (k = 0; k < nb_max; k++) */
/*     { */
/*         fmpz_mat_init(&gamma[k], 2 * g, 2 * g); */
/*     } */

/*     nb_gamma = 0; */
/*     add = 1; */
/*     fmpz_mat_set(cur, mat); */
/*     while (add > 0) */
/*     { */
/*         add = sp2gz_reduce_gamma(cur, gamma + nb_gamma, cur); */
/*         nb_gamma += add; */
/*     } */

/*     /\* Reduce delta one last time and recursive call *\/ */
/*     sp2gz_reduce_delta(cur, &gamma[nb_gamma], cur); */
/*     if (!fmpz_mat_is_one(&gamma[nb_gamma])) */
/*     { */
/*         nb_gamma++; */
/*     } */
/*     if (g > 1) */
/*     { */
/*         fmpz_mat_init(w, 2 * (g - 1), 2 * (g - 1)); */
/*         add = sp2gz_get_parabolic(w, gamma + nb_gamma, cur); */
/*         rec = sp2gz_decompose(&nb_rec, w); */
/*         fmpz_mat_clear(w); */
/*     } */
/*     else */
/*     { */
/*         sp2gz_inv(&gamma[nb_gamma], cur); */
/*         add = (fmpz_mat_is_one(&gamma[nb_gamma]) ? 0 : 1); */
/*     } */

/*     /\* Make final vector *\/ */
/*     *nb = add + nb_rec + nb_gamma; */
/*     res = flint_malloc(*nb * sizeof(fmpz_mat_struct)); */
/*     for (k = 0; k < *nb; k++) */
/*     { */
/*         fmpz_mat_init(&res[k], 2 * g, 2 * g); */
/*     } */

/*     for (k = 0; k < add; k++) */
/*     { */
/*         sp2gz_inv(&res[k], &gamma[nb_gamma + add - 1 - k]); */
/*     } */
/*     for (k = 0; k < nb_rec; k++) */
/*     { */
/*         sp2gz_embed(&res[add + k], &rec[k]); */
/*     } */
/*     for (k = 0; k < nb_gamma; k++) */
/*     { */
/*         sp2gz_inv(&res[add + nb_rec + k], &gamma[nb_gamma - 1 - k]); */
/*     } */

/*     fmpz_mat_clear(cur); */
/*     for (k = 0; k < nb_max; k++) */
/*     { */
/*         fmpz_mat_clear(&gamma[k]); */
/*     } */
/*     flint_free(gamma); */
/*     if (g > 1) */
/*     { */
/*         for (k = 0; k < nb_rec; k++) */
/*         { */
/*             fmpz_mat_clear(&rec[k]); */
/*         } */
/*         flint_free(rec); */
/*     } */
/*     return res; */
/* } */
