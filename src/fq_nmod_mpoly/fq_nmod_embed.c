/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_mat.h"
#include "fmpz.h"
#include "fq_nmod.h"
#include "fq_nmod_poly.h"
#include "fq_nmod_poly_factor.h"
#include "n_poly.h"
#include "fq_nmod_mpoly.h"

/*
    F_q is the "small" field presented as F_p[theta]/f(theta)
    where q = p^m and f is a polynomial over F_p of degree m

    We would like to extend this small field to a "large" field F_q^n
    presented as F_p[phi]/g(phi) where g is a polynomial over F_p of degree m*n

    F_q^n is then isomorphic to F_p[theta][x]/(f(theta), h(theta, x)) for
    some h of degree n in x.

    we compute h and values
        x       as an element of F_p[phi]/g(phi)
        theta   as an element of F_p[phi]/g(phi)
        phi     as an element of F_p[theta][x]

Example:
    with theta printed as # as phi as $, the two embeddings for a
    degree 8 extension F_{3^16} = F_3/($^16+$^10+2) of F_9 = F_3[#]/(#^2+#+2)
    could be:

**** emb[0]:
p: 3
m: 2
n: 8
sm modulus: #^2+#+2
lg modulus: $^16+$^10+2
h: x^8+(2*#+1)*x^6+x^4+(#+1)*x^2+(2*#)
  phi $: x
theta #: 2*$^14+$^10+$^6+$^4+2*$^2+1
      x: $
**** emb[1]:
p: 3
m: 2
n: 8
sm modulus: #^2+#+2
lg modulus: $^16+$^10+2
h: x^8+(#+2)*x^6+x^4+(2*#)*x^2+(#+1)
  phi $: x
theta #: $^14+2*$^10+2*$^6+2*$^4+$^2+1
      x: $
*/


/* the polynomial h is printed with variable x */
/*
static void _embed_print(const bad_fq_nmod_embed_t emb)
{
flint_printf("p: %wu\n", emb->smctx->modulus->mod.n);
flint_printf("m: %wd\n", nmod_poly_degree(emb->smctx->modulus));
flint_printf("n: %wd\n", fq_nmod_poly_degree(emb->h, emb->smctx));

printf("sm modulus: "); nmod_poly_print_pretty(emb->smctx->modulus, emb->smctx->var); printf("\n");
printf("lg modulus: "); nmod_poly_print_pretty(emb->lgctx->modulus, emb->lgctx->var); printf("\n");
printf("h: "); fq_nmod_poly_print_pretty(emb->h, "x", emb->smctx); printf("\n");

printf("  phi %s: ",emb->lgctx->var); fq_nmod_poly_print_pretty(emb->phi_sm, "x", emb->smctx); printf("\n");
printf("theta %s: ",emb->smctx->var); fq_nmod_print_pretty(emb->theta_lg, emb->lgctx); printf("\n");
printf("      x: ");                  fq_nmod_print_pretty(emb->x_lg, emb->lgctx); printf("\n");
}
*/

void bad_fq_nmod_embed_clear(bad_fq_nmod_embed_t emb)
{
    fq_nmod_poly_clear(emb->phi_sm, emb->smctx);
    fq_nmod_poly_clear(emb->h, emb->smctx);
    n_fq_poly_clear(emb->h_as_n_fq_poly);

    fq_nmod_clear(emb->theta_lg, emb->lgctx);
    fq_nmod_clear(emb->x_lg, emb->lgctx);

    nmod_mat_clear(emb->lg_to_sm_mat);
    nmod_mat_clear(emb->sm_to_lg_mat);
}

static void bad_fq_nmod_embed_array_clear(bad_fq_nmod_embed_struct * emb, slong m)
{
    slong i;
    for (i = 0; i < m; i++)
        bad_fq_nmod_embed_clear(emb + i);
}


/* matrices that allow conversion to be done via matrix-vector products */
static void _set_matrices(bad_fq_nmod_embed_t cur)
{
    slong m = fq_nmod_ctx_degree(cur->smctx);
    slong n = fq_nmod_ctx_degree(cur->lgctx);
    slong i;
    n_fq_poly_t phi_as_n_fq_poly, phi_pow, q;
    ulong ** Mrows = cur->lg_to_sm_mat->rows;

    n_fq_poly_init(phi_as_n_fq_poly);
    n_fq_poly_init(phi_pow);
    n_fq_poly_init(q);

    n_fq_poly_set_fq_nmod_poly(phi_as_n_fq_poly, cur->phi_sm, cur->smctx);
    n_fq_poly_one(phi_pow, cur->smctx);

    for (i = 0; i < n; i++)
    {
        n_fq_poly_divrem(q, phi_pow, phi_pow, cur->h_as_n_fq_poly, cur->smctx);
        FLINT_ASSERT(phi_pow->length <= n/m);
        _nmod_vec_set(Mrows[i], phi_pow->coeffs, phi_pow->length*m);
        n_fq_poly_mul(phi_pow, phi_pow, phi_as_n_fq_poly, cur->smctx);
    }

    n_fq_poly_clear(phi_as_n_fq_poly);
    n_fq_poly_clear(phi_pow);
    n_fq_poly_clear(q);

    /* matrix for going from large to small */
    nmod_mat_transpose(cur->lg_to_sm_mat, cur->lg_to_sm_mat);

    /* matrix for going from small to large */
    if (!nmod_mat_inv(cur->sm_to_lg_mat, cur->lg_to_sm_mat))
        flint_throw(FLINT_ERROR, "bad_fq_nmod_embed_array_init: singular matrix");
}

/*
    Initialize an array of m embeddings making an extension of degree n

    This allows for the conversion between

        Fp[phi]/g(phi)  <->   Fp[theta][x]/h(x)

    with the arbitrary choice that phi is mapped to x and x is mapped to phi.
*/
void bad_fq_nmod_embed_array_init(bad_fq_nmod_embed_struct * emb,
                      const fq_nmod_ctx_t bigctx, const fq_nmod_ctx_t smallctx)
{
    slong i, j, k, l;
    bad_fq_nmod_embed_struct * cur;
    fq_nmod_poly_t poly;
    fq_nmod_t t;
    fq_nmod_poly_t poly2;
    fq_nmod_t t2, lead2;
    fq_nmod_t t3;
    fq_nmod_poly_factor_t fac2;
    nmod_mat_t M, Msol;
    fq_nmod_t biggen;
    fmpz_t P;
    ulong lc_inv;
    ulong p = smallctx->modulus->mod.n;
    slong n, m = nmod_poly_degree(smallctx->modulus);

    /* n is the degree of the extension */
    n = nmod_poly_degree(bigctx->modulus);
    FLINT_ASSERT((n%m) == 0);
    n = n/m;

    fmpz_init_set_ui(P, p);

    if (m == 1)
    {
        cur = emb + 0;
        cur->smctx = smallctx;
        cur->lgctx = bigctx;

        fq_nmod_init(cur->theta_lg, bigctx);
        FLINT_ASSERT(1 == nmod_poly_get_coeff_ui(cur->smctx->modulus, 1));
        fq_nmod_set_ui(cur->theta_lg, nmod_poly_get_coeff_ui(
                                              cur->smctx->modulus, 0), bigctx);
        fq_nmod_neg(cur->theta_lg, cur->theta_lg, bigctx);

        fq_nmod_init(cur->x_lg, bigctx);
        fq_nmod_gen(cur->x_lg, bigctx);

        fq_nmod_poly_init(cur->phi_sm, smallctx);
        fq_nmod_poly_gen(cur->phi_sm, smallctx);

        fq_nmod_poly_init(cur->h, smallctx);
        fq_nmod_init(t, smallctx);
        for (i = 0; i < nmod_poly_length(bigctx->modulus); i++)
        {
            fq_nmod_set_ui(t, nmod_poly_get_coeff_ui(bigctx->modulus, i), smallctx);
            fq_nmod_poly_set_coeff(cur->h, i, t, smallctx);
        }

        /* matrix for going from large to small - columns are powers of phi */
        n_fq_poly_init(cur->h_as_n_fq_poly);
        n_fq_poly_set_fq_nmod_poly(cur->h_as_n_fq_poly, cur->h, smallctx);

        nmod_mat_init(cur->lg_to_sm_mat, m*n, m*n, p);
        nmod_mat_init(cur->sm_to_lg_mat, m*n, m*n, p);
        _set_matrices(cur);

        fq_nmod_clear(t, smallctx);
        fmpz_clear(P);

        return;
    }

    /* poly will be bigctx->modulus as a polynomial over smallctx */
    fq_nmod_poly_init(poly, smallctx);
    fq_nmod_init(t, smallctx);
    for (i = 0; i < nmod_poly_length(bigctx->modulus); i++)
    {
        fq_nmod_set_ui(t, nmod_poly_get_coeff_ui(bigctx->modulus, i), smallctx);
        fq_nmod_poly_set_coeff(poly, i, t, smallctx);
    }

    /* poly2 will be smallctx->modulus as a polynomial over bigctx */
    fq_nmod_poly_init(poly2, bigctx);
    fq_nmod_init(t2, bigctx);
    fq_nmod_init(t3, bigctx);
    for (i = 0; i < nmod_poly_length(smallctx->modulus); i++)
    {
        fq_nmod_set_ui(t2, nmod_poly_get_coeff_ui(smallctx->modulus, i), bigctx);
        fq_nmod_poly_set_coeff(poly2, i, t2, bigctx);
    }

    /* poly2 should factor into m linear factors over bigctx*/
    fq_nmod_poly_factor_init(fac2, bigctx);
    fq_nmod_init(lead2, bigctx);
    fq_nmod_poly_factor(fac2, lead2, poly2, bigctx);
    FLINT_ASSERT(fac2->num == m);

    nmod_mat_init(M, m*n, m*n+1, p);
    nmod_mat_init(Msol, m*n+1, 1, p);

    fq_nmod_init(biggen, bigctx);
    fq_nmod_gen(biggen, bigctx);

    for (k = 0; k < m; k++)
    {
        cur = emb + k;
        cur->smctx = smallctx;
        cur->lgctx = bigctx;

        /* we will send x to phi and phi to x */
        fq_nmod_init(cur->x_lg, bigctx);
        fq_nmod_gen(cur->x_lg, bigctx);

        fq_nmod_poly_init(cur->phi_sm, smallctx);
        fq_nmod_poly_gen(cur->phi_sm, smallctx);

        /* theta is determined from a factor of poly2 */
        FLINT_ASSERT(fac2->exp[k] == 1);
        FLINT_ASSERT(fq_nmod_poly_degree(fac2->poly + k, bigctx) == 1);
        fq_nmod_init(cur->theta_lg, bigctx);
        fq_nmod_poly_get_coeff(cur->theta_lg, fac2->poly + k, 0, bigctx);
        fq_nmod_poly_get_coeff(t2, fac2->poly + k, 1, bigctx);
        fq_nmod_inv(t2, t2, bigctx);
        fq_nmod_neg(t2, t2, bigctx);
        fq_nmod_mul(cur->theta_lg, cur->theta_lg, t2, bigctx);

        /* determine h by a nullspace calculation */
        fq_nmod_one(t2, bigctx);
        for (i = 0; i < n; i++)
        {
            fq_nmod_set(t3, t2, bigctx);
            for (j = 0; j < m; j++)
            {
                FLINT_ASSERT(nmod_poly_degree(t3) < m*n);
                for (l = 0; l < m*n; l++)
                {
                    nmod_mat_entry(M, l, m*i + j) = nmod_poly_get_coeff_ui(t3, l);
                }
                fq_nmod_mul(t3, t3, cur->theta_lg, bigctx);
            }
            fq_nmod_mul(t2, t2, biggen, bigctx);
        }

        fq_nmod_pow_ui(t3, biggen, n, bigctx);
        for (l = 0; l < m*n; l++)
        {
            nmod_mat_entry(M, l, m*n) = nmod_poly_get_coeff_ui(t3, l);
        }
        i = nmod_mat_nullspace(Msol, M);
        FLINT_ASSERT(i == 1);

        /* this is the coefficient of x^n in h */
        FLINT_ASSERT(nmod_mat_entry(Msol, m*n, 0) != 0);
        lc_inv = nmod_inv(nmod_mat_entry(Msol, m*n, 0), smallctx->modulus->mod);

        /* set h now */
        fq_nmod_poly_init(cur->h, smallctx);
        for (i = 0; i < n; i++)
        {
            fq_nmod_zero(t, smallctx);
            for (j = 0; j < m; j++)
            {
                nmod_poly_set_coeff_ui(t, j,
                            nmod_mul(lc_inv, nmod_mat_entry(Msol, m*i + j, 0),
                                 smallctx->modulus->mod));
            }
            fq_nmod_poly_set_coeff(cur->h, i, t, smallctx);
        }
        fq_nmod_set_ui(t, 1, smallctx);
        fq_nmod_poly_set_coeff(cur->h, n, t, smallctx);

        /* the set of h's should be the factorization of bigctx->modulus */
        FLINT_ASSERT(fq_nmod_poly_is_irreducible(cur->h, smallctx));
        FLINT_ASSERT(fq_nmod_poly_divides(poly, poly, cur->h, smallctx));

        n_fq_poly_init(cur->h_as_n_fq_poly);
        n_fq_poly_set_fq_nmod_poly(cur->h_as_n_fq_poly, cur->h, smallctx);

        nmod_mat_init(cur->lg_to_sm_mat, m*n, m*n, p);
        nmod_mat_init(cur->sm_to_lg_mat, m*n, m*n, p);
        _set_matrices(cur);
    }
    FLINT_ASSERT(fq_nmod_poly_degree(poly, smallctx) == 0);

    nmod_mat_clear(Msol);
    nmod_mat_clear(M);

    fq_nmod_clear(biggen, bigctx);

    fq_nmod_poly_clear(poly2, bigctx);
    fq_nmod_clear(t2, bigctx);
    fq_nmod_clear(t3, bigctx);

    fq_nmod_poly_factor_clear(fac2, bigctx);
    fq_nmod_clear(lead2, bigctx);

    fq_nmod_poly_clear(poly, smallctx);
    fq_nmod_clear(t, smallctx);

    fmpz_clear(P);
}


/***************** convert Fp[phi]/g(phi) to Fp[theta][x]/h(x) ***************/

/* just matrix-vector multiplication */
void bad_n_fq_embed_lg_to_sm(
    n_fq_poly_t out,        /* poly over smctx */
    const ulong * in,   /* element of lgctx */
    const bad_fq_nmod_embed_t emb)
{
    slong smd = fq_nmod_ctx_degree(emb->smctx);
    slong lgd = fq_nmod_ctx_degree(emb->lgctx);
    slong i;
    const dot_params_t params = _nmod_vec_dot_params(lgd, emb->lgctx->mod);

    n_poly_fit_length(out, lgd);
    for (i = 0; i < lgd; i++)
        out->coeffs[i] = _nmod_vec_dot(emb->lg_to_sm_mat->rows[i], in, lgd,
                                                      emb->lgctx->mod, params);
    FLINT_ASSERT(lgd/smd == emb->h->length - 1);
    out->length = emb->h->length - 1;
    _n_fq_poly_normalise(out, smd);

#if FLINT_WANT_ASSERT
    {
        fq_nmod_t in_;
        fq_nmod_poly_t out_, out_check;

        fq_nmod_init(in_, emb->lgctx);
        fq_nmod_poly_init(out_, emb->smctx);
        fq_nmod_poly_init(out_check, emb->smctx);

        n_fq_get_fq_nmod(in_, in, emb->lgctx);
        n_fq_poly_get_fq_nmod_poly(out_, out, emb->smctx);

        bad_fq_nmod_embed_lg_to_sm(out_check, in_, emb);

        FLINT_ASSERT(fq_nmod_poly_equal(out_check, out_, emb->smctx));

        fq_nmod_clear(in_, emb->lgctx);
        fq_nmod_poly_clear(out_, emb->smctx);
        fq_nmod_poly_clear(out_check, emb->smctx);
    }
#endif
}

/* poorly-implemented old version looking at retirement */
void bad_fq_nmod_embed_lg_to_sm(
    fq_nmod_poly_t out,  /* poly over smctx */
    const fq_nmod_t in,  /* element of lgctx */
    const bad_fq_nmod_embed_t emb)
{
    slong i;
    fq_nmod_poly_t t1;
    fq_nmod_t t2;

    fq_nmod_poly_init(t1, emb->smctx);
    fq_nmod_init(t2, emb->smctx);

    fq_nmod_poly_zero(out, emb->smctx);
    for (i = 0; i < nmod_poly_length(in); i++)
    {
        fq_nmod_poly_pow(t1, emb->phi_sm, i, emb->smctx);
        fq_nmod_set_ui(t2, nmod_poly_get_coeff_ui(in, i), emb->smctx);
        fq_nmod_poly_scalar_mul_fq_nmod(t1, t1, t2, emb->smctx);
        fq_nmod_poly_add(out, out, t1, emb->smctx);
    }
    fq_nmod_poly_rem(out, out, emb->h, emb->smctx);

    fq_nmod_poly_clear(t1, emb->smctx);
    fq_nmod_clear(t2, emb->smctx);
}

void bad_fq_nmod_embed_fq_nmod_lg_to_n_fq_sm(
    n_poly_t out_,  /* poly over smctx */
    const fq_nmod_t in,  /* element of lgctx */
    const bad_fq_nmod_embed_t emb)
{
    fq_nmod_poly_t out;

    fq_nmod_poly_init(out, emb->smctx);

    bad_fq_nmod_embed_lg_to_sm(out, in, emb);
    n_fq_poly_set_fq_nmod_poly(out_, out, emb->smctx);

    fq_nmod_poly_clear(out, emb->smctx);
}

/***************** convert Fp[theta][x]/h(x) to Fp[phi]/g(phi) ***************/

/* just matrix-vector multiplication */
void bad_n_fq_embed_sm_to_lg(
    ulong * out,        /* element of lgctx */
    const n_fq_poly_t in,   /* poly over smctx */
    const bad_fq_nmod_embed_t emb)
{
    slong smd = fq_nmod_ctx_degree(emb->smctx);
    slong lgd = fq_nmod_ctx_degree(emb->lgctx);
    slong i;
    const dot_params_t params = _nmod_vec_dot_params(lgd, emb->lgctx->mod);
    n_poly_stack_t St;  /* TODO: pass the stack in */
    n_fq_poly_struct * q, * in_red;

    n_poly_stack_init(St);

    n_poly_stack_fit_request(St, 2);
    q = n_poly_stack_take_top(St);
    in_red = n_poly_stack_take_top(St);

    n_fq_poly_divrem_(q, in_red, in, emb->h_as_n_fq_poly, emb->smctx, St);

    FLINT_ASSERT(smd*in_red->length <= lgd);

    for (i = 0; i < lgd; i++)
        out[i] = _nmod_vec_dot(emb->sm_to_lg_mat->rows[i], in_red->coeffs,
                                  smd*in_red->length, emb->lgctx->mod, params);

    n_poly_stack_give_back(St, 2);

    n_poly_stack_clear(St);

#if FLINT_WANT_ASSERT
    {
        fq_nmod_t out_, out_check;
        fq_nmod_poly_t in_;

        fq_nmod_init(out_, emb->smctx);
        fq_nmod_init(out_check, emb->smctx);
        fq_nmod_poly_init(in_, emb->lgctx);

        n_fq_get_fq_nmod(out_, out, emb->lgctx);
        n_fq_poly_get_fq_nmod_poly(in_, in, emb->smctx);

        bad_fq_nmod_embed_sm_to_lg(out_check, in_, emb);

        FLINT_ASSERT(fq_nmod_equal(out_check, out_, emb->lgctx));

        fq_nmod_clear(out_, emb->smctx);
        fq_nmod_clear(out_check, emb->smctx);
        fq_nmod_poly_clear(in_, emb->lgctx);
    }
#endif
}

/* poorly-implemented old version */
void bad_fq_nmod_embed_sm_to_lg(
    fq_nmod_t out,            /* element of lgctx */
    const fq_nmod_poly_t in,  /* poly over smctx */
    const bad_fq_nmod_embed_t emb)
{
    slong i, j;
    fq_nmod_poly_t inred;
    fq_nmod_t t1, t2;

    fq_nmod_poly_init(inred, emb->smctx);
    fq_nmod_poly_rem(inred, in, emb->h, emb->smctx);

    fq_nmod_init(t1, emb->lgctx);
    fq_nmod_init(t2, emb->lgctx);

    fq_nmod_zero(out, emb->lgctx);
    for (i = 0; i < fq_nmod_poly_length(inred, emb->smctx); i++)
    {
        nmod_poly_struct * coeff = inred->coeffs + i;
        fq_nmod_zero(t1, emb->lgctx);
        for (j = 0; j < nmod_poly_length(coeff); j++)
        {
            fq_nmod_pow_ui(t2, emb->theta_lg, j, emb->lgctx);
            fq_nmod_mul_ui(t2, t2, nmod_poly_get_coeff_ui(coeff, j), emb->lgctx);
            fq_nmod_add(t1, t1, t2, emb->lgctx);
        }
        fq_nmod_pow_ui(t2, emb->x_lg, i, emb->lgctx);
        fq_nmod_mul(t1, t1, t2, emb->lgctx);
        fq_nmod_add(out, out, t1, emb->lgctx);
    }
    fq_nmod_clear(t1, emb->lgctx);
    fq_nmod_clear(t2, emb->lgctx);
    fq_nmod_poly_clear(inred, emb->smctx);
}

void bad_fq_nmod_embed_n_fq_sm_to_fq_nmod_lg(
    fq_nmod_t out,          /* element of lgctx */
    const n_poly_t in_,     /* poly over smctx */
    const bad_fq_nmod_embed_t emb)
{
    fq_nmod_poly_t in;

    fq_nmod_poly_init(in, emb->smctx);

    n_fq_poly_get_fq_nmod_poly(in, in_, emb->smctx);
    bad_fq_nmod_embed_sm_to_lg(out, in , emb);

    fq_nmod_poly_clear(in, emb->smctx);
}

/**************** convert Fp[theta]/f(theta) to Fp[phi]/g(phi) ***************/

void bad_n_fq_embed_sm_elem_to_lg(
    ulong * out,
    const ulong * in,
    const bad_fq_nmod_embed_t emb)
{
    slong smd = fq_nmod_ctx_degree(emb->smctx);
    slong lgd = fq_nmod_ctx_degree(emb->lgctx);
    slong i;
    const dot_params_t params = _nmod_vec_dot_params(smd, emb->lgctx->mod);

    for (i = 0; i < lgd; i++)
        out[i] = _nmod_vec_dot(emb->sm_to_lg_mat->rows[i], in, smd,
                                                      emb->lgctx->mod, params);
}

void bad_fq_nmod_embed_sm_elem_to_lg(
    fq_nmod_t out,
    const fq_nmod_t in,
    const bad_fq_nmod_embed_t emb)
{
    slong smd = fq_nmod_ctx_degree(emb->smctx);
    slong lgd = fq_nmod_ctx_degree(emb->lgctx);
    slong i;
    const dot_params_t params = _nmod_vec_dot_params(smd, emb->lgctx->mod);

    FLINT_ASSERT(in->length <= smd);

    nmod_poly_fit_length(out, lgd);

    for (i = 0; i < lgd; i++)
    {
        out->coeffs[i] = _nmod_vec_dot(emb->sm_to_lg_mat->rows[i],
                              in->coeffs, in->length, emb->lgctx->mod, params);
    }

    out->length = lgd;
    _nmod_poly_normalise(out);
}

/*****************************************************************************/

bad_fq_nmod_embed_struct *
bad_fq_nmod_mpoly_embed_chooser_init(bad_fq_nmod_mpoly_embed_chooser_t embc,
                  fq_nmod_mpoly_ctx_t ectx, const fq_nmod_mpoly_ctx_t ctx,
                                                        flint_rand_t randstate)
{
    nmod_poly_t ext_modulus;
    fq_nmod_ctx_t ext_fqctx;
    ulong p = ctx->fqctx->modulus->mod.n;
    slong m = nmod_poly_degree(ctx->fqctx->modulus);
    slong n;

    n = WORD(20)/(m*FLINT_BIT_COUNT(p));
    n = FLINT_MAX(n, WORD(2));

    embc->p = p;
    embc->m = m;
    embc->n = n;

    embc->embed = (bad_fq_nmod_embed_struct *) flint_malloc(m*
                                                sizeof(bad_fq_nmod_embed_struct));

    /* init ectx with modulus of degree m*n */
    nmod_poly_init2(ext_modulus, p, m*n + 1);
    nmod_poly_randtest_sparse_irreducible(ext_modulus, randstate, m*n + 1);
    fq_nmod_ctx_init_modulus(ext_fqctx, ext_modulus, "$");
    fq_nmod_mpoly_ctx_init(ectx, ctx->minfo->nvars, ORD_LEX, ext_fqctx);
    fq_nmod_ctx_clear(ext_fqctx);
    nmod_poly_clear(ext_modulus);

    bad_fq_nmod_embed_array_init(embc->embed, ectx->fqctx, ctx->fqctx);

    embc->k = 0;
    return embc->embed + embc->k;
}

void
bad_fq_nmod_mpoly_embed_chooser_clear(bad_fq_nmod_mpoly_embed_chooser_t embc,
                  fq_nmod_mpoly_ctx_t ectx, const fq_nmod_mpoly_ctx_t FLINT_UNUSED(ctx),
                                                        flint_rand_t FLINT_UNUSED(randstate))
{
    bad_fq_nmod_embed_array_clear(embc->embed, embc->m);
    fq_nmod_mpoly_ctx_clear(ectx);
    flint_free(embc->embed);
}


bad_fq_nmod_embed_struct *
bad_fq_nmod_mpoly_embed_chooser_next(bad_fq_nmod_mpoly_embed_chooser_t embc,
                  fq_nmod_mpoly_ctx_t ectx, const fq_nmod_mpoly_ctx_t ctx,
                                                        flint_rand_t randstate)
{
    nmod_poly_t ext_modulus;
    fq_nmod_ctx_t ext_fqctx;
    ulong p = embc->p;
    slong m = embc->m;
    slong n = embc->n;

    embc->k++;
    if (embc->k < m)
        return embc->embed + embc->k;

    n++;
    embc->n = n;
    if (n > 1000)
        return NULL;

    bad_fq_nmod_embed_array_clear(embc->embed, embc->m);
    fq_nmod_mpoly_ctx_clear(ectx);

    /* init ectx with modulus of degree m*n */
    nmod_poly_init2(ext_modulus, p, m*n + 1);
    nmod_poly_randtest_sparse_irreducible(ext_modulus, randstate, m*n + 1);
    fq_nmod_ctx_init_modulus(ext_fqctx, ext_modulus, "$");
    fq_nmod_mpoly_ctx_init(ectx, ctx->minfo->nvars, ORD_LEX, ext_fqctx);
    fq_nmod_ctx_clear(ext_fqctx);
    nmod_poly_clear(ext_modulus);

    bad_fq_nmod_embed_array_init(embc->embed, ectx->fqctx, ctx->fqctx);

    embc->k = 0;
    return embc->embed + embc->k;
}
