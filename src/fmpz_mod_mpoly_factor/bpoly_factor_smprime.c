/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "fmpz_poly_factor.h"
#include "fmpz_mod_vec.h"
#include "fmpz_mod_poly_factor.h"
#include "fmpz_mod_mat.h"
#include "fmpz_mod_mpoly_factor.h"

void fmpz_mod_tpoly_fit_length(
    fmpz_mod_tpoly_t A,
    slong len,
    const fmpz_mod_ctx_t ctx)
{
    slong i = A->alloc;

    if (len <= i)
        return;

    len = FLINT_MAX(len, 2*i);

    A->coeffs = FLINT_ARRAY_REALLOC(A->coeffs, len, fmpz_mod_bpoly_struct);

    for ( ; i < len; i++)
        fmpz_mod_bpoly_init(A->coeffs + i, ctx);

    A->alloc = len;
}

void fmpz_mod_tpoly_clear(fmpz_mod_tpoly_t A, const fmpz_mod_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        fmpz_mod_bpoly_clear(A->coeffs + i, ctx);
    flint_free(A->coeffs);
}



void fmpz_mod_bpoly_reverse_vars(
    fmpz_mod_bpoly_t A,
    const fmpz_mod_bpoly_t B,
    const fmpz_mod_ctx_t ctx)
{
    slong i, j;
    fmpz_mod_bpoly_zero(A, ctx);
    for (i = 0; i < B->length; i++)
    {
        fmpz_mod_poly_struct * Bi = B->coeffs + i;
        for (j = 0; j < Bi->length; j++)
        {
            if (!fmpz_is_zero(Bi->coeffs + j))
            {
                fmpz_mod_bpoly_set_coeff(A, j, i, Bi->coeffs + j, ctx);
            }
        }
    }
}

void fmpz_mod_bpoly_make_monic_series(
    fmpz_mod_bpoly_t A,
    const fmpz_mod_bpoly_t B,
    slong order,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    fmpz_mod_poly_t lcinv;

    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(fmpz_mod_bpoly_is_canonical(B, ctx));

    fmpz_mod_poly_init(lcinv, ctx);
    fmpz_mod_poly_inv_series(lcinv, B->coeffs + B->length - 1, order, ctx);

    fmpz_mod_bpoly_fit_length(A, B->length, ctx);
    for (i = 0; i < B->length; i++)
        fmpz_mod_poly_mullow(A->coeffs + i, B->coeffs + i, lcinv, order, ctx);

    A->length = B->length;
    fmpz_mod_bpoly_normalise(A, ctx);

    fmpz_mod_poly_clear(lcinv, ctx);

    FLINT_ASSERT(fmpz_mod_bpoly_is_canonical(A, ctx));
}

static void fmpz_mod_bpoly_eval(
    fmpz_mod_poly_t E,
    const fmpz_mod_bpoly_t A,
    const fmpz_t alpha,
    const fmpz_mod_ctx_t ctx)
{
    slong i;

    fmpz_mod_poly_zero(E, ctx);

    if (fmpz_is_zero(alpha))
    {
        for (i = A->length - 1; i >= 0; i--)
            if (A->coeffs[i].length > 0)
                fmpz_mod_poly_set_coeff_fmpz(E, i, A->coeffs[i].coeffs + 0, ctx);
        return;
    }

    fmpz_mod_poly_fit_length(E, A->length, ctx);

    for (i = A->length - 1; i >= 0; i--)
        fmpz_mod_poly_evaluate_fmpz(E->coeffs + i, A->coeffs + i, alpha, ctx);

    _fmpz_mod_poly_set_length(E, A->length);
    _fmpz_mod_poly_normalise(E);
}


/****************** lifting **************************************************/

static void _hensel_lift_fac(
    fmpz_mod_bpoly_t G,
    fmpz_mod_bpoly_t H,
    const fmpz_mod_bpoly_t f,
    fmpz_mod_bpoly_t g,
    fmpz_mod_bpoly_t h,
    const fmpz_mod_bpoly_t a,
    const fmpz_mod_bpoly_t b,
    slong p0,
    slong p1,
    const fmpz_mod_ctx_t ctx)
{
    slong i, j;
    fmpz_mod_bpoly_t c, t1, t2, q, r;

    fmpz_mod_bpoly_init(c, ctx);
    fmpz_mod_bpoly_init(t1, ctx);
    fmpz_mod_bpoly_init(t2, ctx);
    fmpz_mod_bpoly_init(q, ctx);
    fmpz_mod_bpoly_init(r, ctx);

    fmpz_mod_bpoly_mul(t1, g, h, ctx);
    fmpz_mod_bpoly_sub(c, f, t1, ctx);
    for (i = 0; i < c->length; i++)
    {
        for (j = 0; j < p0; j++)
            FLINT_ASSERT(j >= c->coeffs[i].length ||
                         fmpz_is_zero(c->coeffs[i].coeffs + j));

        fmpz_mod_poly_shift_right(c->coeffs + i, c->coeffs + i, p0, ctx);
        fmpz_mod_poly_truncate(c->coeffs + i, p1, ctx);
    }

    fmpz_mod_bpoly_mul_series(t1, c, b, p1, ctx);
    fmpz_mod_bpoly_divrem_series(q, r, t1, g, p1, ctx);
    for (i = 0; i < r->length; i++)
        fmpz_mod_poly_shift_left(r->coeffs + i, r->coeffs + i, p0, ctx);
    for (i = 0; i < g->length; i++)
        fmpz_mod_poly_truncate(g->coeffs + i, p0, ctx);
    fmpz_mod_bpoly_add(t1, r, g, ctx);

    fmpz_mod_bpoly_mul_series(t2, c, a, p1, ctx);
    fmpz_mod_bpoly_divrem_series(q, r, t2, h, p1, ctx);
    for (i = 0; i < r->length; i++)
        fmpz_mod_poly_shift_left(r->coeffs + i, r->coeffs + i, p0, ctx);
    for (i = 0; i < h->length; i++)
        fmpz_mod_poly_truncate(h->coeffs + i, p0, ctx);
    fmpz_mod_bpoly_add(t2, r, h, ctx);

    fmpz_mod_bpoly_swap(G, t1, ctx);
    fmpz_mod_bpoly_swap(H, t2, ctx);

#ifdef FLINT_WANT_ASSERT
    fmpz_mod_bpoly_mul(t1, G, H, ctx);
    fmpz_mod_bpoly_sub(c, f, t1, ctx);
    for (i = 0; i < c->length; i++)
        for (j = 0; j < p0 + p1; j++)
            FLINT_ASSERT(j >= c->coeffs[i].length ||
                         fmpz_is_zero(c->coeffs[i].coeffs + j));
#endif

    fmpz_mod_bpoly_clear(c, ctx);
    fmpz_mod_bpoly_clear(t1, ctx);
    fmpz_mod_bpoly_clear(t2, ctx);
    fmpz_mod_bpoly_clear(q, ctx);
    fmpz_mod_bpoly_clear(r, ctx);
}

static void _hensel_lift_inv(
    fmpz_mod_bpoly_t A,
    fmpz_mod_bpoly_t B,
    const fmpz_mod_bpoly_t G,
    const fmpz_mod_bpoly_t H,
    fmpz_mod_bpoly_t a,
    fmpz_mod_bpoly_t b,
    slong p0,
    slong p1,
    const fmpz_mod_ctx_t ctx)
{
    slong i, j;
    fmpz_mod_bpoly_t c, t1, t2, q, r;

    fmpz_mod_bpoly_init(c, ctx);
    fmpz_mod_bpoly_init(t1, ctx);
    fmpz_mod_bpoly_init(t2, ctx);
    fmpz_mod_bpoly_init(q, ctx);
    fmpz_mod_bpoly_init(r, ctx);

    for (i = 0; i < a->length; i++)
        fmpz_mod_poly_truncate(a->coeffs + i, p0, ctx);
    for (i = 0; i < b->length; i++)
        fmpz_mod_poly_truncate(b->coeffs + i, p0, ctx);

    fmpz_mod_bpoly_mul(t1, G, a, ctx);
    fmpz_mod_bpoly_mul(t2, H, b, ctx);
    fmpz_mod_bpoly_add(c, t1, t2, ctx);

    FLINT_ASSERT(c->length > 0);
    for (i = 0; i < c->length; i++)
        fmpz_mod_poly_neg(c->coeffs + i, c->coeffs + i, ctx);
    fmpz_mod_poly_add_si(c->coeffs + 0, c->coeffs + 0, 1, ctx);
    fmpz_mod_bpoly_normalise(c, ctx);

    for (i = 0; i < c->length; i++)
    {
        for (j = 0; j < p0; j++)
            FLINT_ASSERT(j >= c->coeffs[i].length ||
                         fmpz_is_zero(c->coeffs[i].coeffs + j));

        fmpz_mod_poly_shift_right(c->coeffs + i, c->coeffs + i, p0, ctx);
        fmpz_mod_poly_truncate(c->coeffs + i, p1, ctx);
    }

    fmpz_mod_bpoly_mul_series(t1, c, b, p1, ctx);
    fmpz_mod_bpoly_divrem_series(q, r, t1, G, p1, ctx);
    for (i = 0; i < r->length; i++)
        fmpz_mod_poly_shift_left(r->coeffs + i, r->coeffs + i, p0, ctx);
    fmpz_mod_bpoly_add(t1, r, b, ctx);

    fmpz_mod_bpoly_mul_series(t2, c, a, p1, ctx);
    fmpz_mod_bpoly_divrem_series(q, r, t2, H, p1, ctx);
    for (i = 0; i < r->length; i++)
        fmpz_mod_poly_shift_left(r->coeffs + i, r->coeffs + i, p0, ctx);
    fmpz_mod_bpoly_add(t2, r, a, ctx);

    fmpz_mod_bpoly_swap(t1, B, ctx);
    fmpz_mod_bpoly_swap(t2, A, ctx);

#ifdef FLINT_WANT_ASSERT
    fmpz_mod_bpoly_mul(t1, G, A, ctx);
    fmpz_mod_bpoly_mul(t2, H, B, ctx);
    fmpz_mod_bpoly_add(c, t1, t2, ctx);

    FLINT_ASSERT(c->length > 0);
    for (i = 0; i < c->length; i++)
        fmpz_mod_poly_neg(c->coeffs + i, c->coeffs + i, ctx);
    fmpz_mod_poly_add_si(c->coeffs + 0, c->coeffs + 0, 1, ctx);
    fmpz_mod_bpoly_normalise(c, ctx);

    for (i = 0; i < c->length; i++)
        for (j = 0; j < p0 + p1; j++)
            FLINT_ASSERT(j >= c->coeffs[i].length ||
                         fmpz_is_zero(c->coeffs[i].coeffs + j));
#endif

    fmpz_mod_bpoly_clear(c, ctx);
    fmpz_mod_bpoly_clear(t1, ctx);
    fmpz_mod_bpoly_clear(t2, ctx);
    fmpz_mod_bpoly_clear(q, ctx);
    fmpz_mod_bpoly_clear(r, ctx);
}

static void _hensel_lift_tree(
    int opt,
    slong * link,
    fmpz_mod_bpoly_struct * v,
    fmpz_mod_bpoly_struct * w,
    const fmpz_mod_bpoly_t f,
    slong j,
    slong p0,
    slong p1,
    const fmpz_mod_ctx_t ctx)
{
    FLINT_ASSERT(p1 <= p0);

    if (j < 0)
        return;

    if (opt >= 0)
        _hensel_lift_fac(v + j, v + j + 1,
                           f, v + j, v + j + 1, w + j, w + j + 1, p0, p1, ctx);

    if (opt <= 0)
        _hensel_lift_inv(w + j, w + j + 1,
                              v + j, v + j + 1, w + j, w + j + 1, p0, p1, ctx);

    _hensel_lift_tree(opt, link, v, w, v + j, link[j], p0, p1, ctx);
    _hensel_lift_tree(opt, link, v, w, v + j + 1, link[j + 1], p0, p1, ctx);
}

typedef struct {
    slong * link;
    fmpz_mod_bpoly_struct ** lifted_fac;
    fmpz_mod_tpoly_t tmp;
    fmpz_mod_bpoly_t bmp;
    slong r;
    slong fac_lift_order;
    slong inv_lift_order;
    int use_linear;
} fmpz_mod_bpoly_lift_struct;

typedef fmpz_mod_bpoly_lift_struct fmpz_mod_bpoly_lift_t[1];

static void fmpz_mod_bpoly_lift_init(
    fmpz_mod_bpoly_lift_t L,
    const fmpz_mod_ctx_t ctx)
{
    L->link = NULL;
    L->lifted_fac = NULL;
    fmpz_mod_tpoly_init(L->tmp, ctx);
    fmpz_mod_bpoly_init(L->bmp, ctx);
    L->r = 0;
    L->fac_lift_order = 0;
    L->inv_lift_order = 0;
}

static void fmpz_mod_bpoly_lift_clear(
    fmpz_mod_bpoly_lift_t L,
    const fmpz_mod_ctx_t ctx)
{
    flint_free(L->link);
    flint_free(L->lifted_fac);
    fmpz_mod_tpoly_clear(L->tmp, ctx);
    fmpz_mod_bpoly_clear(L->bmp, ctx);
}


static void _fmpz_mod_bpoly_lift_build_tree(
    fmpz_mod_bpoly_lift_t L,
    fmpz_mod_bpoly_struct * local_facs,
    slong r,
    const fmpz_mod_bpoly_t monicA,
    const fmpz_mod_ctx_t ctx)
{
    slong i, j;
    fmpz_mod_bpoly_struct * v, * w;
    slong e[FLINT_BITS+1];
    slong * link;
    fmpz_mod_poly_t d, g, h, a, b;

    fmpz_mod_poly_init(d, ctx);
    fmpz_mod_poly_init(g, ctx);
    fmpz_mod_poly_init(h, ctx);
    fmpz_mod_poly_init(a, ctx);
    fmpz_mod_poly_init(b, ctx);

    FLINT_ASSERT(r > 1);

    L->link = FLINT_ARRAY_REALLOC(L->link, 2*r - 2, slong);
    link = L->link;

    fmpz_mod_tpoly_clear(L->tmp, ctx);
    fmpz_mod_tpoly_init(L->tmp, ctx);
    fmpz_mod_tpoly_fit_length(L->tmp, 2*(2*r - 2), ctx);
    v = L->tmp->coeffs + 0;
    w = L->tmp->coeffs + (2*r - 2);

    for (i = 0; i < r; i++)
    {
        fmpz_mod_bpoly_swap(v + i, local_facs + i, ctx);
        link[i] = -i - 1;
    }

    for (i = r, j = 0; j < 2*r - 4; i++, j += 2)
    {
        slong s, minp, mind;

        minp = j;
        mind = fmpz_mod_bpoly_degree0(v + j, ctx);
        for (s = j + 1; s < i; s++)
        {
            if (fmpz_mod_bpoly_degree0(v + s, ctx) < mind)
            {
                minp = s;
                mind = fmpz_mod_bpoly_degree0(v + s, ctx);
            }
        }
        fmpz_mod_bpoly_swap(v + j, v + minp, ctx);
        FLINT_SWAP(slong, link[j], link[minp]);

        minp = j + 1;
        mind = fmpz_mod_bpoly_degree0(v + j + 1, ctx);
        for (s = j + 2; s < i; s++)
        {
            if (fmpz_mod_bpoly_degree0(v + s, ctx) < mind)
            {
                minp = s;
                mind = fmpz_mod_bpoly_degree0(v + s, ctx);
            }
        }
        fmpz_mod_bpoly_swap(v + j + 1, v + minp, ctx);
        FLINT_SWAP(slong, link[j + 1], link[minp]);

        fmpz_mod_bpoly_mul_series(v + i, v + j, v + j + 1, L->fac_lift_order, ctx);
        link[i] = j;
    }

    for (j = 0; j < 2*r - 2; j += 2)
    {
        fmpz zero = 0;
        fmpz_mod_bpoly_eval(g, v + j, &zero, ctx);
        fmpz_mod_bpoly_eval(h, v + j + 1, &zero, ctx);
        fmpz_mod_poly_xgcd(d, a, b, g, h, ctx);
        if (!fmpz_mod_poly_is_one(d, ctx))
            flint_throw(FLINT_IMPINV, "fmpz_mod_bpoly_lift: bad inverse");
        fmpz_mod_bpoly_set_poly_gen0(w + j, a, ctx);
        fmpz_mod_bpoly_set_poly_gen0(w + j + 1, b, ctx);
    }

    fmpz_mod_poly_clear(d, ctx);
    fmpz_mod_poly_clear(g, ctx);
    fmpz_mod_poly_clear(h, ctx);
    fmpz_mod_poly_clear(a, ctx);
    fmpz_mod_poly_clear(b, ctx);

    for (i = 0; i < 2*r - 2; i++)
        if (-L->link[i] - 1 >= 0)
            L->lifted_fac[-L->link[i] - 1] = v + i;

    for (i = 0, e[i] = L->inv_lift_order; e[i] > 1; i++)
        e[i+1] = (e[i] + 1)/2;
    e[i] = 1;

    for (i--; i >= 0; i--)
        _hensel_lift_tree(-1, L->link, v, w, monicA, 2*r-4, e[i+1], e[i]-e[i+1], ctx);
}

static void _fmpz_mod_bpoly_lift_build_steps(
    fmpz_mod_bpoly_lift_t L,
    const fmpz_mod_ctx_t ctx)
{
    slong i, j, k;
    slong r = L->r;
    slong order = L->fac_lift_order;
    fmpz_mod_bpoly_struct * A, * Bfinal, * U, * B;
    fmpz_mod_poly_struct * s, * Binv;
    fmpz_mod_poly_struct * c, * t;

    FLINT_ASSERT(L->tmp->alloc >= 4*r + 1);
    A = L->tmp->coeffs;
    Bfinal = A + 1;
    U = Bfinal + r;
    B = U + r;

    FLINT_ASSERT(L->bmp->alloc >= 2*r + 5);
    s = L->bmp->coeffs;
    Binv = s + r;
    c = Binv + r;
    t = c + 1;

    for (k = 0; k < r; k++)
    {
        /* s[k] = (prod_{i!=k} B[i].coeffs[0])^-1 (mod B[k].coeffs[0]) */
        fmpz_mod_poly_div(t, A->coeffs + 0, B[k].coeffs + 0, ctx);
        if (!fmpz_mod_poly_invmod(s + k, t, B[k].coeffs + 0, ctx))
            flint_throw(FLINT_IMPINV, "fmpz_mod_bpoly_lift: bad inverse");

        /* set up mul (mod B[k].coeffs[0]) */
        fmpz_mod_poly_reverse(t, B[k].coeffs + 0, B[k].coeffs[0].length, ctx);
        fmpz_mod_poly_inv_series(Binv + k, t, B[k].coeffs[0].length, ctx);
    }

    /* U[0], U[r-1] are not used, zero out both U and Ue */
    for (k = 1; k < r - 1; k++)
    {
        fmpz_mod_bpoly_fit_length(U + k, order, ctx);
        for (i = U[k].length; i < order; i++)
            U[k].coeffs[i].length = 0;
        U[k].length = order;
    }

    if (r > 2)
    {
        for (j = 0; j < order; j++)
        {
            k = r - 2;
            fmpz_mod_poly_zero(U[k].coeffs + j, ctx);
            for (i = 0; i <= j; i++)
            {
                if (i < B[k].length && j - i < B[k + 1].length)
                {
                    fmpz_mod_poly_mul(t, B[k].coeffs + i, B[k + 1].coeffs + j - i, ctx);
                    fmpz_mod_poly_add(U[k].coeffs + j, U[k].coeffs + j, t, ctx);
                }
            }
            for (k--; k > 0; k--)
            {
                fmpz_mod_poly_zero(U[k].coeffs + j, ctx);
                for (i = 0; i <= j; i++)
                {
                    if (i < B[k].length)
                    {
                        fmpz_mod_poly_mul(t, B[k].coeffs + i, U[k + 1].coeffs + j - i, ctx);
                        fmpz_mod_poly_add(U[k].coeffs + j, U[k].coeffs + j, t, ctx);
                    }
                }
            }
        }
    }
}

/* linear lifting has large memory requirements wrt r */
static int _use_linear_cutoff(slong r, slong degx)
{
    return r < 30 + 5*FLINT_BIT_COUNT(degx);
}

static void fmpz_mod_bpoly_lift_start(
    fmpz_mod_bpoly_lift_t L,
    fmpz_mod_poly_struct * local_facs,
    slong r,
    const fmpz_mod_bpoly_t monicA,
    const fmpz_mod_ctx_t ctx)
{
    slong i, k;
    slong degx = fmpz_mod_bpoly_degree0(monicA, ctx);
    fmpz_mod_bpoly_struct * A, * Bfinal, * U, * B;

    FLINT_ASSERT(r > 1);

    L->r = r;
    L->lifted_fac = FLINT_ARRAY_REALLOC(L->lifted_fac, r, fmpz_mod_bpoly_struct*);
    L->fac_lift_order = 1;
    L->inv_lift_order = 1;
    L->use_linear = _use_linear_cutoff(r, degx);

    if (!L->use_linear)
    {
        fmpz_mod_bpoly_struct * new_facs = FLINT_ARRAY_ALLOC(r, fmpz_mod_bpoly_struct);
        for (i = 0; i < r; i++)
        {
            fmpz_mod_bpoly_init(new_facs + i, ctx);
            fmpz_mod_bpoly_set_poly_gen0(new_facs + i, local_facs + i, ctx);
        }

        _fmpz_mod_bpoly_lift_build_tree(L, new_facs, r, monicA, ctx);

        for (i = 0; i < r; i++)
            fmpz_mod_bpoly_clear(new_facs + i, ctx);
        flint_free(new_facs);
    }
    else
    {
        fmpz_mod_tpoly_fit_length(L->tmp, 4*r + 1, ctx);
        A = L->tmp->coeffs;
        Bfinal = A + 1;
        U = Bfinal + r;
        B = U + r;

        fmpz_mod_bpoly_fit_length(L->bmp, 2*r + 5, ctx);

        fmpz_mod_bpoly_fit_length(A, 1, ctx);
        A->length = 1;
        fmpz_mod_poly_one(A->coeffs + 0, ctx);
        for (k = 0; k < r; k++)
        {
            fmpz_mod_bpoly_fit_length(B + k, 1, ctx);
            B[k].length = 1;
            fmpz_mod_poly_set(B[k].coeffs + 0, local_facs + k, ctx);
            fmpz_mod_poly_mul(A->coeffs + 0, A->coeffs + 0, B[k].coeffs + 0, ctx);

            L->lifted_fac[k] = Bfinal + k;
            fmpz_mod_bpoly_reverse_vars(L->lifted_fac[k], B + k, ctx);

            U[k].length = 0;
        }

        _fmpz_mod_bpoly_lift_build_steps(L, ctx);
    }
}

/*
    assuming N is reduced
        combine the factors in L according to the rows of N
        and then replace N by an identity matrix
*/
void fmpz_mod_bpoly_lift_combine(
    fmpz_mod_bpoly_lift_t L,
    fmpz_mod_mat_t N,
    const fmpz_mod_bpoly_t monicA,
    const fmpz_mod_ctx_t ctx)
{
    slong i, j, k, r, degx;
    fmpz_mod_bpoly_struct * A, * Bfinal, * U, * B;
    slong oldr = L->r;
    slong newr = fmpz_mod_mat_nrows(N);
    slong order = L->fac_lift_order;
    fmpz_mod_bpoly_struct * new_facs;
    fmpz_mod_bpoly_t T;

    FLINT_ASSERT(newr > 1);
    FLINT_ASSERT(newr < oldr);
    FLINT_ASSERT(oldr == fmpz_mod_mat_ncols(N));
    FLINT_ASSERT(fmpz_mod_mat_is_reduced(N));

    /* on input we should have a factorization of monicA mod y^order */
#ifdef FLINT_WANT_ASSERT
    {
        fmpz_mod_bpoly_t t1, t2;
        fmpz_mod_bpoly_init(t1, ctx);
        fmpz_mod_bpoly_init(t2, ctx);
        fmpz_mod_bpoly_set(t1, L->lifted_fac[0], ctx);
        for (k = 1; k < L->r; k++)
        {
            fmpz_mod_bpoly_mul_series(t2, t1, L->lifted_fac[k], order, ctx);
            fmpz_mod_bpoly_swap(t1, t2, ctx);
        }
        fmpz_mod_bpoly_sub(t2, monicA, t1, ctx);
        for (i = 0; i < t2->length; i++)
            for (j = 0; j < order; j++)
                FLINT_ASSERT(j >= t2->coeffs[i].length ||
                             fmpz_is_zero(t2->coeffs[i].coeffs + j));
        fmpz_mod_bpoly_clear(t1, ctx);
        fmpz_mod_bpoly_clear(t2, ctx);
    }
#endif

    fmpz_mod_bpoly_init(T, ctx);
    new_facs = FLINT_ARRAY_ALLOC(newr, fmpz_mod_bpoly_struct);
    for (i = 0; i < newr; i++)
    {
        fmpz_mod_bpoly_init(new_facs + i, ctx);
        fmpz_mod_bpoly_one(new_facs + i, ctx);
        for (j = 0; j < oldr; j++)
        {
            if (fmpz_is_zero(fmpz_mod_mat_entry(N, i, j)))
                continue;
            FLINT_ASSERT(fmpz_is_one(fmpz_mod_mat_entry(N, i, j)));
            fmpz_mod_bpoly_mul_series(T, new_facs + i, L->lifted_fac[j], order, ctx);
            fmpz_mod_bpoly_swap(new_facs + i, T, ctx);
        }
    }

    L->r = r = newr;

    degx = fmpz_mod_bpoly_degree0(monicA, ctx);

    /* do not use quadratic lifting if we were not already */
    L->use_linear = L->use_linear || _use_linear_cutoff(r, degx);

    if (!L->use_linear)
    {
        _fmpz_mod_bpoly_lift_build_tree(L, new_facs, newr, monicA, ctx);

        for (i = 0; i < newr; i++)
            fmpz_mod_bpoly_clear(new_facs + i, ctx);
        flint_free(new_facs);
        fmpz_mod_bpoly_clear(T, ctx);
    }
    else
    {
        A = L->tmp->coeffs;
        Bfinal = A + 1;

        fmpz_mod_bpoly_swap(T, A, ctx);
        fmpz_mod_tpoly_clear(L->tmp, ctx);
        fmpz_mod_tpoly_init(L->tmp, ctx);
        fmpz_mod_tpoly_fit_length(L->tmp, 4*r + 1, ctx);
        A = L->tmp->coeffs;
        fmpz_mod_bpoly_swap(A, T, ctx);
        fmpz_mod_bpoly_clear(T, ctx);
        Bfinal = A + 1;
        U = Bfinal + r;
        B = U + r;

        fmpz_mod_bpoly_clear(L->bmp, ctx);
        fmpz_mod_bpoly_init(L->bmp, ctx);
        fmpz_mod_bpoly_fit_length(L->bmp, 2*r + 5, ctx);

        for (i = 0; i < newr; i++)
        {
            L->lifted_fac[i] = Bfinal + i;
            fmpz_mod_bpoly_swap(Bfinal + i, new_facs + i, ctx);
            fmpz_mod_bpoly_clear(new_facs + i, ctx);
        }
        flint_free(new_facs);

        for (k = 0; k < r; k++)
        {
            fmpz_mod_bpoly_reverse_vars(B + k, L->lifted_fac[k], ctx);
            FLINT_ASSERT(B[k].length <= order);
            fmpz_mod_bpoly_fit_length(B + k, order, ctx);
            for (i = B[k].length; i < order; i++)
                fmpz_mod_poly_zero(B[k].coeffs + i, ctx);
        }

        _fmpz_mod_bpoly_lift_build_steps(L, ctx);
    }

    fmpz_mod_mat_clear(N);
    fmpz_mod_mat_init(N, L->r, L->r, fmpz_mod_ctx_modulus(ctx));
    for (i = 0; i < L->r;  i++)
        fmpz_one(fmpz_mod_mat_entry(N, i, i));

    /* on output we should have a factorization of monicA mod y^order */
#ifdef FLINT_WANT_ASSERT
    {
        fmpz_mod_bpoly_t t1, t2;
        fmpz_mod_bpoly_init(t1, ctx);
        fmpz_mod_bpoly_init(t2, ctx);
        fmpz_mod_bpoly_set(t1, L->lifted_fac[0], ctx);
        for (k = 1; k < L->r; k++)
        {
            fmpz_mod_bpoly_mul_series(t2, t1, L->lifted_fac[k], order, ctx);
            fmpz_mod_bpoly_swap(t1, t2, ctx);
        }
        fmpz_mod_bpoly_sub(t2, monicA, t1, ctx);
        for (i = 0; i < t2->length; i++)
            for (j = 0; j < order; j++)
                FLINT_ASSERT(j >= t2->coeffs[i].length ||
                             fmpz_is_zero(t2->coeffs[i].coeffs + j));
        fmpz_mod_bpoly_clear(t1, ctx);
        fmpz_mod_bpoly_clear(t2, ctx);
    }
#endif
}

static void fmpz_mod_bpoly_lift_continue(
    fmpz_mod_bpoly_lift_t L,
    const fmpz_mod_bpoly_t monicA,
    slong order,
    const fmpz_mod_ctx_t ctx)
{
    slong i, j, k;
    slong r = L->r;
    fmpz_mod_bpoly_struct * A, * Bfinal, * U, * Ue, * B;
    fmpz_mod_poly_struct * s, * Binv;
    fmpz_mod_poly_struct * c, * t, * ce, * vk;

    if (order <= L->fac_lift_order)
        return;

    if (!L->use_linear)
    {
        fmpz_mod_bpoly_struct * v = L->tmp->coeffs + 0;
        fmpz_mod_bpoly_struct * w = L->tmp->coeffs + (2*r - 2);
        slong e[FLINT_BITS+1];

        FLINT_ASSERT(1 <= L->inv_lift_order);
        FLINT_ASSERT(L->fac_lift_order <= 2*L->inv_lift_order);

        for (i = 0, e[i] = order; e[i] > L->fac_lift_order; i++)
            e[i+1] = (e[i] + 1)/2;
        e[i]   = L->fac_lift_order;
        e[i+1] = L->inv_lift_order;

        if (e[i+1] < e[i])
            _hensel_lift_tree(-1, L->link, v, w, monicA, 2*r-4, e[i+1], e[i]-e[i+1], ctx);

        for (i--; i > 0; i--)
            _hensel_lift_tree(0, L->link, v, w, monicA, 2*r-4, e[i+1], e[i]-e[i+1], ctx);

        _hensel_lift_tree(1, L->link, v, w, monicA, 2*r-4, e[1], e[0]-e[1], ctx);

        L->fac_lift_order = e[0];
        L->inv_lift_order = e[1];
        return;
    }

    A = L->tmp->coeffs;
    Bfinal = A + 1;
    U = Ue = Bfinal + r;
    B = U + r;

    s = L->bmp->coeffs;
    Binv = s + r;
    c = Binv + r;
    t = c + 1;
    ce = t + 1;
    vk = ce + 1;

    /* tack on reversal of monicA */
    for (i = 0; i < monicA->length; i++)
    {
        fmpz_mod_poly_struct * Bi = monicA->coeffs + i;
        j = FLINT_MIN(Bi->length, order);
        for (j--; j >= L->fac_lift_order; j--)
            fmpz_mod_bpoly_set_coeff(A, j, i, Bi->coeffs + j, ctx);
    }

    /* tack on zeros to the B[k] */
    for (k = 0; k < r; k++)
    {
        fmpz_mod_bpoly_fit_length(B + k, order, ctx);

        for (i = B[k].length; i < order; i++)
        {
            fmpz_mod_poly_zero(B[k].coeffs + i, ctx);
        }

        /* U[0] is not used, zero out both U and Ue */
        if (k > 0)
        {
            fmpz_mod_bpoly_fit_length(U + k, order, ctx);
            for (i = U[k].length; i < order; i++)
                U[k].coeffs[i].length = 0;
            U[k].length = order;
        }
    }

    if (r > 2)
    {
        for (j = L->fac_lift_order; j < order; j++)
        {
            k = r - 2;
            fmpz_mod_poly_zero(U[k].coeffs + j, ctx);
            for (i = 0; i <= j; i++)
            {
                if (i < B[k].length && j - i < B[k + 1].length)
                {
                    fmpz_mod_poly_mul(t, B[k].coeffs + i, B[k + 1].coeffs + j - i, ctx);
                    fmpz_mod_poly_add(U[k].coeffs + j, U[k].coeffs + j, t, ctx);
                }
            }
            for (k--; k > 0; k--)
            {
                fmpz_mod_poly_zero(U[k].coeffs + j, ctx);
                for (i = 0; i <= j; i++)
                {
                    if (i < B[k].length)
                    {
                        fmpz_mod_poly_mul(t, B[k].coeffs + i, U[k + 1].coeffs + j - i, ctx);
                        fmpz_mod_poly_add(U[k].coeffs + j, U[k].coeffs + j, t, ctx);
                    }
                }
            }

            if (j < A->length)
                fmpz_mod_poly_set(c, A->coeffs + j, ctx);
            else
                fmpz_mod_poly_zero(c, ctx);

            for (i = 0; i <= j; i++)
            {
                if (i < B[0].length)
                {
                    fmpz_mod_poly_mul(t, B[0].coeffs + i, U[1].coeffs + j - i, ctx);
                    fmpz_mod_poly_sub(c, c, t, ctx);
                }
            }

            if (fmpz_mod_poly_is_zero(c, ctx))
                continue;

            for (k = r - 1; k >= 0; k--)
            {
                fmpz_mod_poly_rem(t, c, B[k].coeffs + 0, ctx);
                fmpz_mod_poly_mulmod_preinv(vk, s + k, t, B[k].coeffs + 0, Binv + k, ctx);
                if (!fmpz_mod_poly_is_zero(vk, ctx))
                {
                    fmpz_mod_poly_add(B[k].coeffs + j, B[k].coeffs + j, vk, ctx);
                    if (!fmpz_mod_poly_is_zero(B[k].coeffs + j, ctx))
                        B[k].length = FLINT_MAX(B[k].length, j + 1);
                }

                /* correct the U's */
                if (k > r - 2)
                {
                    fmpz_mod_poly_swap(ce, vk, ctx);
                }
                else if (k > 0)
                {
                    fmpz_mod_poly_struct * p;
                    fmpz_mod_poly_mul(t, B[k].coeffs + 0, ce, ctx);
                    p = (k == r - 2) ? B[k + 1].coeffs : U[k + 1].coeffs;
                    fmpz_mod_poly_mul(ce, p, vk, ctx);
                    fmpz_mod_poly_add(ce, ce, t, ctx);
                    fmpz_mod_poly_add(U[k].coeffs + j, U[k].coeffs + j, ce, ctx);
                }
            }
        }
    }
    else
    {
        FLINT_ASSERT(r == 2);

        for (j = L->fac_lift_order; j < order; j++)
        {
            if (j < A->length)
                fmpz_mod_poly_set(c, A->coeffs + j, ctx);
            else
                fmpz_mod_poly_zero(c, ctx);

            for (i = FLINT_MIN(j, B[0].length - 1); i >= 0; i--)
            {
                fmpz_mod_poly_mul(t, B[0].coeffs + i, B[1].coeffs + j - i, ctx);
                fmpz_mod_poly_sub(c, c, t, ctx);
            }

            if (fmpz_mod_poly_is_zero(c, ctx))
                continue;

            for (k = 0; k < r; k++)
            {
                fmpz_mod_poly_rem(t, c, B[k].coeffs + 0, ctx);
                fmpz_mod_poly_mulmod_preinv(vk, s + k, t, B[k].coeffs + 0, Binv + k, ctx);
                if (!fmpz_mod_poly_is_zero(vk, ctx))
                {
                    fmpz_mod_poly_add(B[k].coeffs + j, B[k].coeffs + j, vk, ctx);
                    if (!fmpz_mod_poly_is_zero(B[k].coeffs + j, ctx))
                        B[k].length = FLINT_MAX(B[k].length, j + 1);
                }
            }
        }
    }

    L->fac_lift_order = order;

    for (k = 0; k < r; k++)
        fmpz_mod_bpoly_reverse_vars(Bfinal + k, B + k, ctx);

#ifdef FLINT_WANT_ASSERT
    {
        fmpz_mod_bpoly_t t1, t2;
        fmpz_mod_bpoly_init(t1, ctx);
        fmpz_mod_bpoly_init(t2, ctx);
        fmpz_mod_bpoly_set(t1, Bfinal + 0, ctx);
        for (k = 1; k < r; k++)
        {
            fmpz_mod_bpoly_mul_series(t2, t1, Bfinal + k, order, ctx);
            fmpz_mod_bpoly_swap(t1, t2, ctx);
        }
        fmpz_mod_bpoly_sub(t2, monicA, t1, ctx);

        for (i = 0; i < t2->length; i++)
        {
            for (j = 0; j < FLINT_MIN(order, t2->coeffs[i].length); j++)
            {
                FLINT_ASSERT(fmpz_is_zero(t2->coeffs[i].coeffs + j));
            }
        }
        fmpz_mod_bpoly_clear(t1, ctx);
        fmpz_mod_bpoly_clear(t2, ctx);
    }
#endif
}


/************* lattice reduction ********************************************/

/*
    The rows of N, if they are 0-1, correspond to combinations of the g that
    give factors of A. Compute {A/g[i]*g[i]'}_i and use the CLD bounds to
    try to make N smaller.
*/
static void _lattice(
    fmpz_mod_mat_t N,
    fmpz_mod_bpoly_struct * const * g,
    slong r,
    slong lift_order,
    slong * CLD,
    slong * lattice_order,
    const fmpz_mod_bpoly_t A,
    const fmpz_mod_ctx_t ctx)
{
    slong i, j, k;
    fmpz_mod_bpoly_t Q, R, dg;
    fmpz_mod_bpoly_struct * ld;
    fmpz_mod_mat_t M, T1, T2;
    fmpz * trow;

    trow = _fmpz_vec_init(r);
    fmpz_mod_bpoly_init(Q, ctx);
    fmpz_mod_bpoly_init(R, ctx);
    fmpz_mod_bpoly_init(dg, ctx);
    ld = FLINT_ARRAY_ALLOC(r, fmpz_mod_bpoly_struct);
    for (i = 0; i < r; i++)
    {
        fmpz_mod_bpoly_init(ld + i, ctx);
        fmpz_mod_bpoly_divrem_series(Q, R, A, g[i], lift_order, ctx);
        FLINT_ASSERT(R->length == 0);
        fmpz_mod_bpoly_derivative_gen0(R, g[i], ctx);
        fmpz_mod_bpoly_mul_series(ld + i, Q, R, lift_order, ctx);
    }

    for (k = 0; k + 1 < A->length; k++)
    {
        slong nrows = fmpz_mod_mat_nrows(N);
        slong lower = FLINT_MAX(CLD[k], *lattice_order);

        FLINT_ASSERT(nrows > 0);

        /*
            consider powers y^j for which
            j >= lattice_order (j < lattice_order has already been added)
            j >= CLD[k]
            j < lift_order
        */

        if (lift_order <= lower)
            continue;

        fmpz_mod_mat_init(M, lift_order - lower, nrows, fmpz_mod_ctx_modulus(ctx));

        for (j = lower; j < lift_order; j++)
        {
            for (i = 0; i < r; i++)
                fmpz_mod_bpoly_get_coeff(trow + i, ld + i, k, j, ctx);

            for (i = 0; i < nrows; i++)
                _fmpz_mod_vec_dot(fmpz_mod_mat_entry(M, j - lower, i),
                                                trow, N->mat->rows[i], r, ctx);
        }

        fmpz_mod_mat_init_nullspace_tr(T1, M, ctx);

        fmpz_mod_mat_init(T2, fmpz_mod_mat_nrows(T1), fmpz_mod_mat_ncols(N), fmpz_mod_ctx_modulus(ctx));
        fmpz_mod_mat_mul(T2, T1, N);
        fmpz_mod_mat_swap(T2, N);
        fmpz_mod_mat_rref(NULL, N);

        fmpz_mod_mat_clear(M);
        fmpz_mod_mat_clear(T1);
        fmpz_mod_mat_clear(T2);
    }

    _fmpz_vec_clear(trow, r);
    fmpz_mod_bpoly_clear(Q, ctx);
    fmpz_mod_bpoly_clear(R, ctx);
    fmpz_mod_bpoly_clear(dg, ctx);
    for (i = 0; i < r; i++)
        fmpz_mod_bpoly_clear(ld + i, ctx);
    flint_free(ld);

    *lattice_order = lift_order;
}


/**************** recombination **********************************************/

/*
    First multiply the g[i] into the gprod[i] according the the rows of N.
    Then, run zassenhaus on the gprod[i] as factors of A.
*/
static int _zassenhaus(
    const zassenhaus_prune_t zas,
    slong limit,
    fmpz_mod_tpoly_t F,
    const fmpz_t malpha,
    const fmpz_mod_mat_t N,
    fmpz_mod_bpoly_struct * const * g,
    slong r,
    slong order,
    const fmpz_mod_bpoly_t A,
    const fmpz_mod_ctx_t ctx)
{
    slong total_deg;
    int success;
    fmpz_mod_bpoly_t Q, R, t1, t2;
    fmpz_mod_poly_t cont;
    slong i, j, k, len, nrows = fmpz_mod_mat_nrows(N);
    slong * subset;
    fmpz_mod_bpoly_struct * gprod;
    fmpz_mod_bpoly_struct * f;
    fmpz_mod_bpoly_t A_copy;
    int is_simple_check = (limit == 1 && r == fmpz_mod_mat_nrows(N));

    FLINT_ASSERT(fmpz_mod_mat_is_reduced(N));
    FLINT_ASSERT(fmpz_mod_mat_ncols(N) == r);

    fmpz_mod_poly_init(cont, ctx);
    fmpz_mod_bpoly_init(Q, ctx);
    fmpz_mod_bpoly_init(R, ctx);
    fmpz_mod_bpoly_init(t1, ctx);
    fmpz_mod_bpoly_init(t2, ctx);
    fmpz_mod_bpoly_init(A_copy, ctx);
    gprod = FLINT_ARRAY_ALLOC(nrows, fmpz_mod_bpoly_struct);
    subset = FLINT_ARRAY_ALLOC(nrows, slong);
    for (i = 0; i < nrows; i++)
    {
        subset[i] = i;
        fmpz_mod_bpoly_init(gprod + i, ctx);
        fmpz_mod_bpoly_one(gprod + i, ctx);
        for (j = 0; j < r; j++)
        {
            if (fmpz_is_zero(fmpz_mod_mat_entry(N, i, j)))
                continue;
            FLINT_ASSERT(fmpz_is_one(fmpz_mod_mat_entry(N, i, j)));
            fmpz_mod_bpoly_mul_series(t1, gprod + i, g[j], order, ctx);
            fmpz_mod_bpoly_swap(gprod + i, t1, ctx);
        }
    }

    f = (fmpz_mod_bpoly_struct *) A;

    len = nrows;
    for (k = 1; k <= len/2; k++)
    {
        if (k > limit)
        {
            success = 0;
            goto cleanup;
        }

        zassenhaus_subset_first(subset, len, k);
        while (1)
        {
            total_deg = 0;
            for (i = 0; i < len; i++)
            {
                if (subset[i] >= 0)
                    total_deg += gprod[subset[i]].length - 1;
            }

            if (!zassenhaus_prune_degree_is_possible(zas, total_deg))
            {
                if (!zassenhaus_subset_next(subset, len))
                    break;
                continue;
            }

            FLINT_ASSERT(f->length > 0);
            fmpz_mod_bpoly_set_poly_gen1(t1, f->coeffs + f->length - 1, ctx);
            for (i = 0; i < len; i++)
            {
                if (subset[i] >= 0)
                {
                    fmpz_mod_bpoly_mul_series(t2, t1, gprod + subset[i], order, ctx);
                    fmpz_mod_bpoly_swap(t1, t2, ctx);
                }
            }

            fmpz_mod_bpoly_make_primitive(cont, t1, ctx);
            if (fmpz_mod_bpoly_divides(Q, f, t1, ctx))
            {
                fmpz_mod_bpoly_taylor_shift_gen1(t1, t1, malpha, ctx);
                fmpz_mod_tpoly_fit_length(F, F->length + 1, ctx);
                fmpz_mod_bpoly_swap(F->coeffs + F->length, t1, ctx);
                F->length++;
                f = A_copy;
                fmpz_mod_bpoly_swap(f, Q, ctx);

                len -= k;
                if (!zassenhaus_subset_next_disjoint(subset, len + k))
                    break;
            }
            else if (is_simple_check)
            {
                success = 0;
                goto cleanup;
            }
            else
            {
                if (!zassenhaus_subset_next(subset, len))
                    break;
            }
        }
    }

    if (f->length > 1)
    {
        fmpz_mod_tpoly_fit_length(F, F->length + 1, ctx);
        fmpz_mod_bpoly_taylor_shift_gen1(F->coeffs + F->length, f, malpha, ctx);
        F->length++;
    }
    else
    {
        FLINT_ASSERT(f->length == 1);
        FLINT_ASSERT(fmpz_mod_poly_is_one(f->coeffs + 0, ctx));
    }

    success = 1;

cleanup:

    for (i = 0; i < nrows; i++)
        fmpz_mod_bpoly_clear(gprod + i, ctx);
    flint_free(gprod);

    flint_free(subset);

    fmpz_mod_poly_clear(cont, ctx);
    fmpz_mod_bpoly_clear(Q, ctx);
    fmpz_mod_bpoly_clear(R, ctx);
    fmpz_mod_bpoly_clear(t1, ctx);
    fmpz_mod_bpoly_clear(t2, ctx);
    fmpz_mod_bpoly_clear(A_copy, ctx);

    return success;
}


/*****************************************************************************/

/*
    x = gen(0), y = gen(1).
    A is supposed to be separable wrt x.
    Put the content of A wrt x in c, and the factors in F.
    Return 1 for success, i.e. a good small prime (y - alpha) was found.
    If allow_shift is false, only (y - 0) is tried.

    TODO: copy this precision strategy to the other bpoly factorers

    The helpers fmpz_mod_bpoly_lift_{start|continue|combine} are used:
        start:      start the lift mod y^0
        continue:   lift up to mod y^n
        combine:    when the lattice work has proven several factors to be
                    grouped together, combine these and start over with fewer
                    local factors. ex: if N = [1 1 0 1 0],
                                              [0 0 1 0 1]
                                       then combine the five factors f1,...,f5
                                       into two f1*f2*f4, f3*f5
*/
int fmpz_mod_bpoly_factor_smprime(
    fmpz_mod_poly_t c,     /* poly in y */
    fmpz_mod_tpoly_t F,
    fmpz_mod_bpoly_t A,    /* clobbered */
    int allow_shift,
    const fmpz_mod_ctx_t ctx)
{
    int success;
    slong i, r;
    slong Alenx, Aleny;
    slong final_order, lift_order, lattice_order;
    slong * CLD;
    fmpz_mod_poly_t Aeval;
    fmpz_t alpha_best, alpha_tmp, malpha_best;
    fmpz_mod_poly_factor_t local_fac_best, local_fac_tmp;
    int local_fac_tries = 0;
    fmpz_mod_bpoly_t monicA;
    fmpz_mod_mat_t N;
    slong zas_limit;
    zassenhaus_prune_t zas;
    fmpz_mod_bpoly_lift_t L;

    fmpz_mod_bpoly_make_primitive(c, A, ctx);

    Alenx = A->length;

    FLINT_ASSERT(Alenx > 1);

    fmpz_init(alpha_tmp);
    fmpz_init(alpha_best);
    fmpz_init(malpha_best);

    fmpz_mod_poly_init(Aeval, ctx);
    fmpz_mod_poly_factor_init(local_fac_best, ctx);
    fmpz_mod_poly_factor_init(local_fac_tmp, ctx);
    fmpz_mod_bpoly_init(monicA, ctx);
    fmpz_mod_mat_init(N, 0, 0, fmpz_mod_ctx_modulus(ctx));
    CLD = FLINT_ARRAY_ALLOC(Alenx, slong);
    zassenhaus_prune_init(zas);
    fmpz_mod_bpoly_lift_init(L, ctx);

    Aleny = 0;
    for (i = 0; i < Alenx; i++)
    {
        Aleny = FLINT_MAX(Aleny, A->coeffs[i].length);
        CLD[i] = A->coeffs[i].length;
    }
    mpoly_bivar_cld_bounds(CLD, Alenx);

    zassenhaus_prune_set_degree(zas, Alenx - 1);

    fmpz_zero(alpha_tmp);
    fmpz_zero(alpha_best);
    goto got_alpha;

next_alpha:

    fmpz_add_ui(alpha_tmp, alpha_tmp, 1);

    if (!allow_shift || fmpz_cmp(alpha_tmp, fmpz_mod_ctx_modulus(ctx)) >= 0)
    {
        if (local_fac_best->num > 0)
            goto doit;
        success = 0;
        goto cleanup;
    }

got_alpha:

    fmpz_mod_bpoly_eval(Aeval, A, alpha_tmp, ctx);

    /* if killed leading coeff, get new alpha */
    if (Aeval->length != Alenx)
        goto next_alpha;

    /* note the constant term of Aeval can be zero */

    fmpz_mod_poly_factor(local_fac_tmp, Aeval, ctx);
    r = local_fac_tmp->num;

    zassenhaus_prune_start_add_factors(zas);
    for (i = 0; i < r; i++)
        zassenhaus_prune_add_factor(zas, local_fac_tmp->poly[i].length - 1,
                                         local_fac_tmp->exp[i]);
    zassenhaus_prune_end_add_factors(zas);

    if (r < 2 && local_fac_tmp->exp[0] == 1)
        goto irreducible;
    if (zassenhaus_prune_must_be_irreducible(zas))
        goto irreducible;

    /* if multiple factors, get new alpha */
    for (i = 0; i < r; i++)
    {
        FLINT_ASSERT(local_fac_tmp->poly[i].length > 1);
        FLINT_ASSERT(fmpz_is_one(fmpz_mod_poly_lead(local_fac_tmp->poly + i, ctx)));
        if (local_fac_tmp->exp[i] != 1)
            goto next_alpha;
    }

    /* done if A is constant in y */
    if (Aleny < 2)
    {
        fmpz_mod_tpoly_fit_length(F, r, ctx);
        F->length = r;
        for (i = 0; i < r; i++)
            fmpz_mod_bpoly_set_poly_gen0(F->coeffs + i,
                                         local_fac_tmp->poly + i, ctx);
        success = 1;
        goto cleanup;
    }

    /* alpha_tmp & local_fac_tmp are good; update best */
    if (local_fac_best->num < 1 || local_fac_best->num > local_fac_tmp->num)
    {
        fmpz_set(alpha_best, alpha_tmp);
        fmpz_mod_poly_factor_swap(local_fac_best, local_fac_tmp, ctx);
    }

    if (++local_fac_tries < 2)
        goto next_alpha;

doit:

    fmpz_mod_bpoly_taylor_shift_gen1(A, A, alpha_best, ctx);

    /* local_fac_best is a factorization mod (y - alpha_best) */
    r = local_fac_best->num;

    /* precision for constructing true factors */
    final_order = Aleny;

    /* precision for lifted local factors */
    lift_order = Aleny;
    for (i = 0; i < Alenx - 1; i++)
        if (CLD[i] > 0 && lift_order > CLD[i])
            lift_order = CLD[i];
    lift_order = lift_order + 4;

    /* lift up to y^lift_order */
    fmpz_mod_bpoly_make_monic_series(monicA, A, lift_order, ctx);
    fmpz_mod_bpoly_lift_start(L, local_fac_best->poly, r, monicA, ctx);
    fmpz_mod_bpoly_lift_continue(L, monicA, lift_order, ctx);

    /* the rows of N give the combinations of local factors */
    fmpz_mod_mat_clear(N);
    fmpz_mod_mat_init(N, r, r, fmpz_mod_ctx_modulus(ctx));
    for (i = 0; i < r; i++)
        fmpz_one(fmpz_mod_mat_entry(N, i, i));

    /* size limit on subsets in zassenhaus combination */
    zas_limit = 1;

    lattice_order = 0;
    _lattice(N, L->lifted_fac, L->r, lift_order, CLD, &lattice_order, A, ctx);
    if (fmpz_mod_mat_nrows(N) < 2)
        goto irreducible_shift;
    if (!fmpz_mod_mat_is_reduced(N))
        goto increase;
    if (fmpz_mod_mat_nrows(N) < fmpz_mod_mat_ncols(N)/4*3)
        fmpz_mod_bpoly_lift_combine(L, N, monicA, ctx);

try_zas:

    /* zassenhaus only make sense if N is a nice 0-1 mat */
    FLINT_ASSERT(fmpz_mod_mat_is_reduced(N));

    while (fmpz_mod_mat_nrows(N) > 2 && 2*L->fac_lift_order < final_order)
    {
        lift_order = 2*L->fac_lift_order;
        fmpz_mod_bpoly_make_monic_series(monicA, A, lift_order, ctx);
        fmpz_mod_bpoly_lift_continue(L, monicA, lift_order, ctx);
        _lattice(N, L->lifted_fac, L->r, lift_order, CLD, &lattice_order, A, ctx);
        if (fmpz_mod_mat_nrows(N) < 2)
            goto irreducible_shift;
        if (!fmpz_mod_mat_is_reduced(N))
            goto increase;
        if (fmpz_mod_mat_nrows(N) < fmpz_mod_mat_ncols(N)/4*3)
            fmpz_mod_bpoly_lift_combine(L, N, monicA, ctx);
    }

    if (L->fac_lift_order < final_order)
    {
        lift_order = final_order;
        fmpz_mod_bpoly_make_monic_series(monicA, A, lift_order, ctx);
        fmpz_mod_bpoly_lift_continue(L, monicA, lift_order, ctx);
    }

    /* combine local factors according the rows of N, then by subsets */
    F->length = 0;
    fmpz_mod_neg(malpha_best, alpha_best, ctx);
    success = _zassenhaus(zas, zas_limit, F, malpha_best, N,
                                     L->lifted_fac, L->r, final_order, A, ctx);
    if (success)
        goto cleanup;

    /* first attempt failed, try subsets of size 1 or 2 from now on */
    zas_limit = 2;

more:

    /* increase precision until N is a nice 0-1 mat */
    _lattice(N, L->lifted_fac, L->r, lift_order, CLD, &lattice_order, A, ctx);
    if (fmpz_mod_mat_nrows(N) < 2)
        goto irreducible_shift;
    if (!fmpz_mod_mat_is_reduced(N))
        goto increase;
    if (fmpz_mod_mat_nrows(N) < fmpz_mod_mat_ncols(N)/4*3)
        fmpz_mod_bpoly_lift_combine(L, N, monicA, ctx);
    goto try_zas;

increase:

    if (lift_order < final_order)
        lift_order += 4 + lift_order/2;
    else
        lift_order += 1 + lift_order/8;

    fmpz_mod_bpoly_make_monic_series(monicA, A, lift_order, ctx);
    fmpz_mod_bpoly_lift_continue(L, monicA, lift_order, ctx);
    goto more;

cleanup:

    fmpz_mod_bpoly_lift_clear(L, ctx);

    flint_free(CLD);

    fmpz_clear(alpha_tmp);
    fmpz_clear(alpha_best);
    fmpz_clear(malpha_best);

    fmpz_mod_mat_clear(N);
    fmpz_mod_poly_clear(Aeval, ctx);
    fmpz_mod_poly_factor_clear(local_fac_best, ctx);
    fmpz_mod_poly_factor_clear(local_fac_tmp, ctx);
    fmpz_mod_bpoly_clear(monicA, ctx);

    zassenhaus_prune_clear(zas);

    return success;

irreducible_shift:

    fmpz_mod_neg(malpha_best, alpha_best, ctx);
    fmpz_mod_bpoly_taylor_shift_gen1(A, A, malpha_best, ctx);

irreducible:

    fmpz_mod_tpoly_fit_length(F, 1, ctx);
    F->length = 1;
    fmpz_mod_bpoly_swap(F->coeffs + 0, A, ctx);
    success = 1;
    goto cleanup;
}
