/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "fmpq.h"
#include "fmpz_poly.h"
#include "arb_fmpz_poly.h"
#include "fmpq_poly.h"
#include "fmpz_poly_factor.h"
#include "qqbar.h"
#include "gr_vec.h"
#include "gr_mat.h"
#include "fmpz_mat.h"
#include "ulong_extras.h"
#include "nmod_mat.h"
#include "fmpz_mpoly_q.h"
#include "gr_tower.h"
#include "gr_tower/impl.h"

/* [F_k : F_0] */
#define _degree_at(T, k) gr_tower_degree_at(T, k)

/* Coordinates of x in F_k with respect to the monomial basis over F_0,
   ordered with the lowest level varying fastest. res must have room for
   [F_k : F_0] entries. */
int
gr_tower_get_coeffs_at(gr_ptr res, gr_srcptr x, slong k, gr_tower_t T)
{
    int status = GR_SUCCESS;

    if (k == 0)
    {
        return gr_set(res, x, T->base);
    }
    else
    {
        gr_ctx_struct * ctx = gr_tower_field_at(T, k);
        gr_ctx_struct * below = gr_tower_field_at(T, k - 1);
        slong i, d = gr_tower_step_degree(T, k), D = _degree_at(T, k - 1);
        slong sz = T->base->sizeof_elem;
        gr_poly_t xr;

        gr_poly_init(xr, below);
        status |= gr_poly_quotient_get_poly(xr, x, ctx);

        for (i = 0; i < d && status == GR_SUCCESS; i++)
        {
            if (i < xr->length)
                status |= gr_tower_get_coeffs_at(GR_ENTRY(res, i * D, sz), gr_poly_coeff_srcptr(xr, i, below), k - 1, T);
            else
                status |= _gr_vec_zero(GR_ENTRY(res, i * D, sz), D, T->base);
        }

        gr_poly_clear(xr, below);
    }

    return status;
}

int
gr_tower_set_coeffs_at(gr_ptr res, gr_srcptr coeffs, slong k, gr_tower_t T)
{
    int status = GR_SUCCESS;

    if (k == 0)
    {
        return gr_set(res, coeffs, T->base);
    }
    else
    {
        gr_ctx_struct * below = gr_tower_field_at(T, k - 1);
        slong i, d = gr_tower_step_degree(T, k), D = _degree_at(T, k - 1);
        slong sz = T->base->sizeof_elem;
        gr_poly_struct * poly = res;

        gr_poly_fit_length(poly, d, below);
        for (i = 0; i < d && status == GR_SUCCESS; i++)
            status |= gr_tower_set_coeffs_at(gr_poly_coeff_ptr(poly, i, below), GR_ENTRY(coeffs, i * D, sz), k - 1, T);
        _gr_poly_set_length(poly, d, below);
        _gr_poly_normalise(poly, below);
    }

    return status;
}

int
gr_tower_multiplication_matrix(gr_mat_t res, gr_srcptr x, gr_tower_t T)
{
    gr_ctx_struct * top = gr_tower_field(T);
    slong D = gr_tower_degree(T);
    slong i, j;
    gr_ptr v, w, e;
    int status = GR_SUCCESS;

    if (gr_mat_nrows(res, T->base) != D || gr_mat_ncols(res, T->base) != D)
        return GR_DOMAIN;

    v = gr_heap_init_vec(D, T->base);   /* unit vector */
    w = gr_heap_init_vec(D, T->base);
    GR_TMP_INIT(e, top);

    for (j = 0; j < D && status == GR_SUCCESS; j++)
    {
        status |= _gr_vec_zero(v, D, T->base);
        status |= gr_one(GR_ENTRY(v, j, T->base->sizeof_elem), T->base);
        status |= gr_tower_set_coeffs_at(e, v, T->length, T);
        status |= gr_mul(e, e, x, top);
        status |= gr_tower_get_coeffs_at(w, e, T->length, T);

        /* column j = coordinates of x * e_j */
        for (i = 0; i < D; i++)
            status |= gr_set(gr_mat_entry_ptr(res, i, j, T->base), GR_ENTRY(w, i, T->base->sizeof_elem), T->base);
    }

    gr_heap_clear_vec(v, D, T->base);
    gr_heap_clear_vec(w, D, T->base);
    GR_TMP_CLEAR(e, top);

    return status;
}

int
gr_tower_charpoly(gr_poly_t res, gr_srcptr x, gr_tower_t T)
{
    gr_mat_t M;
    slong D = gr_tower_degree(T);
    int status;

    gr_mat_init(M, D, D, T->base);
    status = gr_tower_multiplication_matrix(M, x, T);
    if (status == GR_SUCCESS)
        status = gr_mat_charpoly(res, M, T->base);
    gr_mat_clear(M, T->base);
    return status;
}

/*
    Monic polynomial of least degree over F_0 annihilating x in the ring
    F_n, computed by incremental Gaussian elimination on the coordinate
    vectors of 1, x, x^2, ...

    Rows of E hold reduced vectors with pivots; rows of C hold the
    corresponding combinations of powers of x.
*/
/*
    The same over the rationals, through the flat representation: the
    coordinate vectors of the powers of x are the rows of an integer
    matrix (each scaled by its denominator), the degree of the
    annihilating polynomial is its rank, and the dependency is a kernel
    vector. Returns GR_UNABLE if the tower has transcendental generators
    or is too large.
*/

static int
_gr_tower_annihilating_poly_flat(gr_poly_t res, gr_srcptr x, gr_tower_t T)
{
    gr_tower_flat_struct * F = &T->flat;
    slong D = gr_tower_degree(T);
    slong n, r, k, i, nvars;
    fmpz_mat_t A, B, N;
    fmpz * scale;
    fmpz_mpoly_q_t xf, pw;
    ulong * exp;
    slong * stride, * var_of;
    int status = GR_SUCCESS;

    if (!GR_TOWER_BASE_IS_CONSTS(T) || !GR_TOWER_HAS_CAP(T, GR_TOWER_CAP_RATIONAL) || D > GR_TOWER_OPTION(T, GR_TOWER_OPT_ANNIHILATING_FLAT_LIMIT))
        return GR_UNABLE;

    gr_tower_flat_ensure(F);
    nvars = F->cap;
    fmpz_mpoly_q_init(xf, F->mctx);
    fmpz_mpoly_q_init(pw, F->mctx);

    if (gr_tower_flat_set_nested_at(xf, x, T->length, F) != GR_SUCCESS)
    {
        fmpz_mpoly_q_clear(xf, F->mctx);
        fmpz_mpoly_q_clear(pw, F->mctx);
        return GR_UNABLE;
    }

    /* index of a monomial in the basis: the lowest level varies fastest */
    stride = flint_malloc(sizeof(slong) * (T->length + 1));
    var_of = flint_malloc(sizeof(slong) * (T->length + 1));
    stride[1] = 1;
    for (k = 1; k <= T->length; k++)
    {
        var_of[k] = GR_TOWER_FLAT_VAR(F, k);
        if (k < T->length)
            stride[k + 1] = stride[k] * gr_tower_step_degree(T, k);
    }

    fmpz_mat_init(A, D + 1, D);
    scale = _fmpz_vec_init(D + 1);
    exp = flint_malloc(sizeof(ulong) * nvars);

    r = -1;
    fmpz_mpoly_q_one(pw, F->mctx);
    for (n = 0; n <= D && status == GR_SUCCESS; n++)
    {
        const fmpz_mpoly_struct * num;

        /* the degree of the minimal polynomial divides D: at those n,
           check whether x^n already depends on the lower powers (the
           powers up to the degree suffice, and their coordinates are
           smaller than those of the higher powers) */
        if (n > 0 && n < D && D % n == 0)
        {
            /* cheap filter: the rank modulo a prime */
            fmpz_mat_t W;
            nmod_mat_t Wp;
            slong rp;
            fmpz_mat_window_init(W, A, 0, 0, n, D);
            /* (a prime near 2^(FLINT_BITS - 2): 2^62 + 135 on 64-bit machines) */
            nmod_mat_init(Wp, n, D, n_nextprime(UWORD(1) << (FLINT_BITS - 2), 1));
            fmpz_mat_get_nmod_mat(Wp, W);
            rp = nmod_mat_rank(Wp);
            nmod_mat_clear(Wp);
            if (rp < n)
            {
                /* exact check: a kernel vector of the transpose */
                fmpz_mat_t Wt, Nt;
                fmpz_mat_init(Wt, D, n);
                fmpz_mat_init(Nt, n, n);
                fmpz_mat_transpose(Wt, W);
                if (fmpz_mat_nullspace(Nt, Wt) >= 1)
                    r = n - 1;
                fmpz_mat_clear(Wt);
                fmpz_mat_clear(Nt);
            }
            fmpz_mat_window_clear(W);
            if (r >= 0)
                break;
        }

        status |= gr_tower_flat_reduce(pw, F);
        if (status != GR_SUCCESS)
            break;
        num = fmpz_mpoly_q_numref(pw);
        if (!fmpz_mpoly_is_fmpz(fmpz_mpoly_q_denref(pw), F->mctx))
        {
            status = GR_UNABLE;
            break;
        }
        fmpz_mpoly_get_fmpz(scale + n, fmpz_mpoly_q_denref(pw), F->mctx);

        for (i = 0; i < num->length; i++)
        {
            slong idx = 0;
            fmpz_mpoly_get_term_exp_ui(exp, num, i, F->mctx);
            for (k = 1; k <= T->length; k++)
                idx += exp[var_of[k]] * stride[k];
            if (idx >= D)
            {
                status = GR_UNABLE;   /* not reduced: should not happen */
                break;
            }
            fmpz_mpoly_get_term_coeff_fmpz(fmpz_mat_entry(A, n, idx), num, i, F->mctx);
        }

        if (n < D)
            fmpz_mpoly_q_mul(pw, pw, xf, F->mctx);
    }

    if (status == GR_SUCCESS)
    {
        /* the powers 1, ..., x^(r-1) are independent and x^r depends on them */
        if (r < 0)
            r = fmpz_mat_rank(A);

        fmpz_mat_init(B, D, r + 1);
        fmpz_mat_init(N, r + 1, r + 1);
        for (n = 0; n <= r; n++)
            for (i = 0; i < D; i++)
                fmpz_set(fmpz_mat_entry(B, i, n), fmpz_mat_entry(A, n, i));

        if (fmpz_mat_nullspace(N, B) != 1)
        {
            status = GR_UNABLE;
        }
        else
        {
            fmpq_t c;
            fmpq_init(c);
            gr_poly_fit_length(res, r + 1, T->base);
            for (n = 0; n <= r && status == GR_SUCCESS; n++)
            {
                /* the rows are scale[n] times the coordinates of x^n */
                fmpz_mul(fmpq_numref(c), fmpz_mat_entry(N, n, 0), scale + n);
                fmpz_one(fmpq_denref(c));
                status |= gr_set_fmpq(gr_poly_coeff_ptr(res, n, T->base), c, T->base);
            }
            _gr_poly_set_length(res, r + 1, T->base);
            _gr_poly_normalise(res, T->base);
            status |= gr_poly_make_monic(res, res, T->base);
            fmpq_clear(c);
        }

        fmpz_mat_clear(B);
        fmpz_mat_clear(N);
    }

    fmpz_mat_clear(A);
    _fmpz_vec_clear(scale, D + 1);
    flint_free(exp);
    flint_free(stride);
    flint_free(var_of);
    fmpz_mpoly_q_clear(xf, F->mctx);
    fmpz_mpoly_q_clear(pw, F->mctx);

    return status;
}

int
_gr_tower_annihilating_poly(gr_poly_t res, gr_srcptr x, gr_tower_t T)
{
    gr_ctx_struct * top = gr_tower_field(T);
    gr_ctx_struct * F = T->base;
    slong D = gr_tower_degree(T);
    slong sz = F->sizeof_elem;
    slong r, i, j, n;
    gr_mat_t E, C;
    slong * pivot;
    gr_ptr v, c, xpow, s;
    int status = GR_SUCCESS;
    int found = 0;

    status = _gr_tower_annihilating_poly_flat(res, x, T);
    if (status == GR_SUCCESS)
        return status;
    status = GR_SUCCESS;

    /* (the coordinates of the powers of x fill a D x D matrix with
       entries whose size grows with D: not attempted for large D) */
    if (D > GR_TOWER_OPTION(T, GR_TOWER_OPT_ANNIHILATING_FLAT_LIMIT))
        return GR_UNABLE;

    gr_mat_init(E, D + 1, D, F);
    gr_mat_init(C, D + 1, D + 1, F);
    pivot = flint_malloc(sizeof(slong) * (D + 1));
    v = gr_heap_init_vec(D, F);
    c = gr_heap_init_vec(D + 1, F);
    GR_TMP_INIT(xpow, top);
    GR_TMP_INIT(s, F);

    status |= gr_one(xpow, top);

    for (n = 0; n <= D && status == GR_SUCCESS; n++)
    {
        /* v = coordinates of x^n, c = e_n */
        status |= gr_tower_get_coeffs_at(v, xpow, T->length, T);
        status |= _gr_vec_zero(c, D + 1, F);
        status |= gr_one(GR_ENTRY(c, n, sz), F);

        /* reduce v against the existing rows */
        for (r = 0; r < n && status == GR_SUCCESS; r++)
        {
            gr_srcptr vp = GR_ENTRY(v, pivot[r], sz);

            if (gr_is_zero(vp, F) != T_TRUE)
            {
                /* copy the scalar: vp aliases an entry of v */
                status |= gr_set(s, vp, F);
                status |= _gr_vec_submul_scalar(v, gr_mat_entry_srcptr(E, r, 0, F), D, s, F);
                status |= _gr_vec_submul_scalar(c, gr_mat_entry_srcptr(C, r, 0, F), D + 1, s, F);
            }
        }

        if (status != GR_SUCCESS)
            break;

        /* find a pivot */
        for (j = 0; j < D; j++)
            if (gr_is_zero(GR_ENTRY(v, j, sz), F) != T_TRUE)
                break;

        if (j == D)
        {
            /* dependency: sum c_i x^i = 0 */
            gr_poly_fit_length(res, n + 1, F);
            for (i = 0; i <= n; i++)
                status |= gr_set(gr_poly_coeff_ptr(res, i, F), GR_ENTRY(c, i, sz), F);
            _gr_poly_set_length(res, n + 1, F);
            _gr_poly_normalise(res, F);
            status |= gr_poly_make_monic(res, res, F);
            found = 1;
            break;
        }

        /* normalize so that the pivot is 1, and store */
        pivot[n] = j;
        status |= gr_inv(s, GR_ENTRY(v, j, sz), F);
        status |= _gr_vec_mul_scalar(v, v, D, s, F);
        status |= _gr_vec_mul_scalar(c, c, D + 1, s, F);
        status |= _gr_vec_set(gr_mat_entry_ptr(E, n, 0, F), v, D, F);
        status |= _gr_vec_set(gr_mat_entry_ptr(C, n, 0, F), c, D + 1, F);

        /* keep the pivot column clear in the earlier rows (not strictly
           necessary for correctness since we reduce sequentially) */

        status |= gr_mul(xpow, xpow, x, top);
    }

    if (status == GR_SUCCESS && !found)
        status = GR_UNABLE;

    gr_mat_clear(E, F);
    gr_mat_clear(C, F);
    flint_free(pivot);
    gr_heap_clear_vec(v, D, F);
    gr_heap_clear_vec(c, D + 1, F);
    GR_TMP_CLEAR(xpow, top);
    GR_TMP_CLEAR(s, F);

    return status;
}

/* (the enclosure of a top-field element, for the factor selection) */
typedef struct { gr_srcptr x; gr_tower_struct * T; } _get_z_elem_arg;

static int
_get_z_elem(acb_t z, slong prec, void * arg)
{
    _get_z_elem_arg * a = (_get_z_elem_arg *) arg;
    return gr_tower_get_acb(z, a->x, prec, a->T);
}

int
gr_tower_get_fmpz_poly_minpoly(fmpz_poly_t res, gr_srcptr x, gr_tower_t T)
{
    gr_poly_t p;
    fmpq_poly_t pq;
    fmpz_poly_t pz;
    fmpz_poly_factor_t fac;
    slong i;
    int status = GR_SUCCESS;

    if (!GR_TOWER_BASE_IS_CONSTS(T) || !GR_TOWER_HAS_CAP(T, GR_TOWER_CAP_RATIONAL))
        return GR_UNABLE;

    gr_poly_init(p, T->base);
    fmpq_poly_init(pq);
    fmpz_poly_init(pz);
    fmpz_poly_factor_init(fac);

    status = _gr_tower_annihilating_poly(p, x, T);

    if (status == GR_SUCCESS)
    {
        fmpq_poly_fit_length(pq, p->length);
        for (i = 0; i < p->length; i++)
            fmpq_poly_set_coeff_fmpq(pq, i, gr_poly_coeff_srcptr(p, i, T->base));
        fmpq_poly_get_numerator(pz, pq);
        fmpz_poly_factor(fac, pz);

        /* select the irreducible factor vanishing at the enclosure of x */
        {
            _get_z_elem_arg a = { x, T };
            slong which = _gr_tower_select_fmpz_factor(fac, _get_z_elem, &a, T);
            if (which >= 0)
                fmpz_poly_set(res, fac->p + which);
            else
                status = GR_UNABLE;
        }
    }

    gr_poly_clear(p, T->base);
    fmpq_poly_clear(pq);
    fmpz_poly_clear(pz);
    fmpz_poly_factor_clear(fac);

    return status;
}

int
gr_tower_get_qqbar(qqbar_t res, gr_srcptr x, gr_tower_t T)
{
    fmpz_poly_t p;
    acb_t z, z2;
    mag_t rad;
    slong prec;
    int status;
    int found = 0;

    fmpz_poly_init(p);
    acb_init(z);
    acb_init(z2);
    mag_init(rad);

    status = gr_tower_get_fmpz_poly_minpoly(p, x, T);

    if (status == GR_SUCCESS)
    {
        if (fmpz_poly_degree(p) == 1)
        {
            fmpq_t q;
            fmpq_init(q);
            fmpz_neg(fmpq_numref(q), p->coeffs);
            fmpz_set(fmpq_denref(q), p->coeffs + 1);
            fmpq_canonicalise(q);
            qqbar_set_fmpq(res, q);
            fmpq_clear(q);
            found = 1;
        }

        for (prec = GR_TOWER_DEFAULT_PREC; !found && prec <= 1000000 && status == GR_SUCCESS; prec *= 2)
        {
            slong prec2;

            status |= gr_tower_get_acb(z, x, prec, T);

            for (prec2 = prec / 2; !found && prec2 < 2 * prec; prec2 *= 2)
            {
                acb_set(z2, z);
                acb_get_mag(rad, z);
                mag_mul_2exp_si(rad, rad, -prec2);
                acb_add_error_mag(z2, rad);

                if (qqbar_set_fmpz_poly_root(res, p, z2, 2 * prec2))
                    found = 1;
            }
        }

        if (!found && status == GR_SUCCESS)
            status = GR_UNABLE;
    }

    fmpz_poly_clear(p);
    acb_clear(z);
    acb_clear(z2);
    mag_clear(rad);

    return status;
}
