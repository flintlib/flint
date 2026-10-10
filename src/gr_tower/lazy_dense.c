/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Lazy fields: dense products in number fields and dense forms of elements. */

/* (for recursive mutexes: pthread_mutexattr_settype) */
#define _GNU_SOURCE

#include "gr_tower/lazy_impl.h"

/* -------------------------------------------------------------------- */
/* dense polynomial and matrix products (number fields)                  */
/* -------------------------------------------------------------------- */

/* the common tower of n elements (merging as needed; the trivial tower
   if they are all rational), or NULL */
static gr_tower_flat_struct *
_gr_tower_lazy_common_vec(gr_tower_lazy_elem_struct * const * e, slong n, gr_ctx_t ctx)
{
    gr_tower_flat_struct * U = NULL;
    slong i, lu = 0;

    for (i = 0; i < n; i++)
    {
        const gr_tower_lazy_elem_struct * ei = _gr_tower_lazy_flat_view(e[i]);
        if (ei->level == 0)
            continue;
        if (U == NULL)
        {
            U = ei->F;
            lu = ei->level;
        }
        else if (ei->F == U)
            lu = FLINT_MAX(lu, ei->level);
        else
        {
            U = _gr_tower_lazy_common_tower(U, lu, ei->F, ei->level, ctx);
            if (U == NULL)
                return NULL;
            lu = U->T->num_gens;
        }
    }

    return (U == NULL) ? LAZY(ctx)->trivial : U;
}

/* t (initialized in the current context of U) = e mapped into U */
static int
_gr_tower_lazy_get_in(fmpz_mpoly_q_t t, gr_tower_lazy_elem_t e, gr_tower_flat_struct * U, gr_ctx_t ctx)
{
    e = _gr_tower_lazy_flat_view(e);
    if (e->F == U)
    {
        fmpz_mpoly_q_set(t, &e->elem.flat.data, U->mctx);
        return GR_SUCCESS;
    }
    if (e->level == 0)
    {
        fmpz_t c;
        fmpz_init(c);
        fmpz_mpoly_get_fmpz(c, fmpz_mpoly_q_numref(&e->elem.flat.data), e->elem.flat.mctx);
        fmpz_mpoly_set_fmpz(fmpz_mpoly_q_numref(t), c, U->mctx);
        fmpz_mpoly_get_fmpz(c, fmpz_mpoly_q_denref(&e->elem.flat.data), e->elem.flat.mctx);
        fmpz_mpoly_set_fmpz(fmpz_mpoly_q_denref(t), c, U->mctx);
        fmpz_clear(c);
        return GR_SUCCESS;
    }
    return _gr_tower_lazy_map_element(t, &e->elem.flat.data, e->elem.flat.mctx, e->F, e->level, U, ctx);
}

/* the elements in a common tower U: tt (allocated, initialized in U's
   context); returns U, or NULL if not applicable (nothing allocated) */
static gr_tower_flat_struct *
_gr_tower_lazy_vec_in_common(fmpz_mpoly_q_struct ** tt, gr_tower_lazy_elem_struct * const * e, slong n, gr_ctx_t ctx)
{
    gr_tower_flat_struct * U = _gr_tower_lazy_common_vec(e, n, ctx);
    fmpz_mpoly_q_struct * t;
    slong i;
    int status = GR_SUCCESS;

    if (U == NULL || U == LAZY(ctx)->trivial || !_gr_tower_flat_dense_applicable(U))
        return NULL;

    gr_tower_flat_ensure(U);
    t = flint_malloc(sizeof(fmpz_mpoly_q_struct) * n);
    for (i = 0; i < n; i++)
        fmpz_mpoly_q_init(t + i, U->mctx);
    for (i = 0; i < n && status == GR_SUCCESS; i++)
        status = _gr_tower_lazy_get_in(t + i, e[i], U, ctx);

    if (status != GR_SUCCESS)
    {
        for (i = 0; i < n; i++)
            fmpz_mpoly_q_clear(t + i, U->mctx);
        flint_free(t);
        return NULL;
    }

    *tt = t;
    return U;
}

#define LAZY_DENSE_POLY_MIN_LEN 4

int
_gr_tower_lazy_poly_mullow(gr_tower_lazy_elem_struct * res, const gr_tower_lazy_elem_struct * poly1, slong len1,
    const gr_tower_lazy_elem_struct * poly2, slong len2, slong n, gr_ctx_t ctx)
{
    gr_tower_lazy_elem_struct ** e;
    fmpz_mpoly_q_struct * t, * r;
    fmpz_mpoly_q_struct ** px, ** py, ** pr;
    gr_tower_flat_struct * U;
    slong i;
    int ok = 0;

    if (FLINT_MIN(len1, len2) >= LAZY_DENSE_POLY_MIN_LEN)
    {
        e = flint_malloc(sizeof(gr_tower_lazy_elem_struct *) * (len1 + len2));
        for (i = 0; i < len1; i++)
            e[i] = (gr_tower_lazy_elem_struct *) (poly1 + i);
        for (i = 0; i < len2; i++)
            e[len1 + i] = (gr_tower_lazy_elem_struct *) (poly2 + i);

        U = _gr_tower_lazy_vec_in_common(&t, e, len1 + len2, ctx);
        if (U != NULL)
        {
            r = flint_malloc(sizeof(fmpz_mpoly_q_struct) * n);
            px = flint_malloc(sizeof(fmpz_mpoly_q_struct *) * (len1 + len2 + n));
            py = px + len1;
            pr = px + len1 + len2;
            for (i = 0; i < len1 + len2; i++)
                px[i] = t + i;
            for (i = 0; i < n; i++)
            {
                fmpz_mpoly_q_init(r + i, U->mctx);
                pr[i] = r + i;
            }

            ok = _gr_tower_flat_poly_mullow_dense(pr, px, len1, py, len2, n, U);

            if (ok)
                for (i = 0; i < n; i++)
                    _gr_tower_lazy_install(res + i, U, r + i, ctx);

            for (i = 0; i < n; i++)
                fmpz_mpoly_q_clear(r + i, U->mctx);
            for (i = 0; i < len1 + len2; i++)
                fmpz_mpoly_q_clear(t + i, U->mctx);
            flint_free(r);
            flint_free(t);
            flint_free(px);
        }
        flint_free(e);
    }

    if (ok)
        return GR_SUCCESS;

    return _gr_poly_mullow_generic(res, poly1, len1, poly2, len2, n, ctx);
}

#define LAZY_DENSE_MAT_MIN 3

int
_gr_tower_lazy_mat_mul(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    slong r = gr_mat_nrows(A, ctx), s = gr_mat_ncols(A, ctx), c = gr_mat_ncols(B, ctx);
    slong i, j;
    int ok = 0;

    if (gr_mat_nrows(B, ctx) != s || gr_mat_nrows(C, ctx) != r || gr_mat_ncols(C, ctx) != c)
        return GR_DOMAIN;

    if (r >= LAZY_DENSE_MAT_MIN && s >= LAZY_DENSE_MAT_MIN && c >= LAZY_DENSE_MAT_MIN)
    {
        gr_tower_lazy_elem_struct ** e = flint_malloc(sizeof(gr_tower_lazy_elem_struct *) * (r * s + s * c));
        fmpz_mpoly_q_struct * t, * w;
        fmpz_mpoly_q_struct ** pa, ** pb, ** pc;
        gr_tower_flat_struct * U;

        for (i = 0; i < r; i++)
            for (j = 0; j < s; j++)
                e[i * s + j] = (gr_tower_lazy_elem_struct *) gr_mat_entry_srcptr(A, i, j, ctx);
        for (i = 0; i < s; i++)
            for (j = 0; j < c; j++)
                e[r * s + i * c + j] = (gr_tower_lazy_elem_struct *) gr_mat_entry_srcptr(B, i, j, ctx);

        U = _gr_tower_lazy_vec_in_common(&t, e, r * s + s * c, ctx);
        if (U != NULL)
        {
            w = flint_malloc(sizeof(fmpz_mpoly_q_struct) * r * c);
            pa = flint_malloc(sizeof(fmpz_mpoly_q_struct *) * (r * s + s * c + r * c));
            pb = pa + r * s;
            pc = pb + s * c;
            for (i = 0; i < r * s + s * c; i++)
                pa[i] = t + i;
            for (i = 0; i < r * c; i++)
            {
                fmpz_mpoly_q_init(w + i, U->mctx);
                pc[i] = w + i;
            }

            ok = _gr_tower_flat_mat_mul_dense(pc, pa, pb, r, s, c, U);

            if (ok)
                for (i = 0; i < r; i++)
                    for (j = 0; j < c; j++)
                        _gr_tower_lazy_install(gr_mat_entry_ptr(C, i, j, ctx), U, w + i * c + j, ctx);

            for (i = 0; i < r * c; i++)
                fmpz_mpoly_q_clear(w + i, U->mctx);
            for (i = 0; i < r * s + s * c; i++)
                fmpz_mpoly_q_clear(t + i, U->mctx);
            flint_free(w);
            flint_free(t);
            flint_free(pa);
        }
        flint_free(e);
    }

    if (ok)
        return GR_SUCCESS;

    return gr_mat_mul_generic(C, A, B, ctx);
}

/*
    Algorithm selection for determinants and linear systems. Over number
    fields with univariate moduli (the dense arithmetic applies), the
    fraction arithmetic of Gaussian elimination produces elements of
    large height when the degree D of the field is large compared to the
    size n of the matrix, and the division-free Berkowitz algorithm
    (determinants) or fraction-free elimination (systems) wins; for small
    degrees, elimination (recursive LU, using the fast matrix products)
    is fastest. Rational matrices go to fmpq_mat. Measured on random
    matrices with small entries in Q(zeta_N) and Q(sqrt(-163)): see the
    design notes.
*/
static slong
_gr_tower_lazy_mat_field_degree(const gr_mat_t A, gr_ctx_t ctx, int * rational)
{
    slong r = gr_mat_nrows(A, ctx), c = gr_mat_ncols(A, ctx), i, j;
    gr_tower_lazy_elem_struct ** e;
    gr_tower_flat_struct * U;
    slong D = -1;

    *rational = 0;
    if (r == 0 || c == 0)
        return -1;

    e = flint_malloc(sizeof(gr_tower_lazy_elem_struct *) * r * c);
    for (i = 0; i < r; i++)
        for (j = 0; j < c; j++)
            e[i * c + j] = (gr_tower_lazy_elem_struct *) gr_mat_entry_srcptr(A, i, j, ctx);
    U = _gr_tower_lazy_common_vec(e, r * c, ctx);
    flint_free(e);

    if (U == LAZY(ctx)->trivial)
        *rational = 1;
    else if (U != NULL && _gr_tower_flat_dense_applicable(U))
        D = gr_tower_degree(U->T);

    return D;
}

int
_gr_tower_lazy_mat_det(gr_tower_lazy_elem_t res, const gr_mat_t A, gr_ctx_t ctx)
{
    slong n = gr_mat_nrows(A, ctx), D, i, j;
    int rational;

    if (n != gr_mat_ncols(A, ctx))
        return GR_DOMAIN;

    if (n <= 4)
        return gr_mat_det_cofactor(res, A, ctx);

    D = _gr_tower_lazy_mat_field_degree(A, ctx, &rational);

    if (rational)
    {
        fmpq_mat_t Q;
        fmpq_t d;
        fmpq_mat_init(Q, n, n);
        fmpq_init(d);
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
            {
                const gr_tower_lazy_elem_struct * x = gr_mat_entry_srcptr(A, i, j, ctx);
                if (LAZY_REPR(x) == LAZY_REPR_RATIONAL)
                    fmpq_set(fmpq_mat_entry(Q, i, j), &x->elem.q);
                else
                {
                    x = _gr_tower_lazy_flat_view(x);
                    (void) fmpz_mpoly_q_get_fmpq(fmpq_mat_entry(Q, i, j), &x->elem.flat.data, x->elem.flat.mctx);
                }
            }
        fmpq_mat_det(d, Q);
        _gr_tower_lazy_set_fmpq(res, d, ctx);
        fmpq_mat_clear(Q);
        fmpq_clear(d);
        return GR_SUCCESS;
    }

    if (D >= 1 && D < n)
        return gr_mat_det_lu(res, A, ctx);

    /* (otherwise division-free: with transcendental generators, LU
       accumulates denominators) */
    return gr_mat_det_berkowitz(res, A, ctx);
}

int
_gr_tower_lazy_mat_nonsingular_solve(gr_mat_t X, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    slong n = gr_mat_nrows(A, ctx), D;
    int rational;

    D = _gr_tower_lazy_mat_field_degree(A, ctx, &rational);

    if (D >= 1 && 2 * D >= n)
        return gr_mat_nonsingular_solve_fflu(X, A, B, ctx);

    return gr_mat_nonsingular_solve_lu(X, A, B, ctx);
}

/* -------------------------------------------------------------------- */
/* dense forms: elements of number fields of one generator               */
/* -------------------------------------------------------------------- */

/*
    An element of a tower F whose first generator a (definition order 0)
    is algebraic and proven, of degree 2 <= d <= GR_TOWER_OPT_DENSE_FORM_DEGREE_LIMIT,
    with a monic integral modulus M, and which involves only a (level 1)
    with an integer denominator, may be kept as a canonical fmpq_poly in
    a of length < d, referring to an immutable descriptor (a, M) of F.
    Arithmetic between such elements of the same descriptor and tower,
    and with rational elements, needs neither the towers nor the lock when
    the result already lives in that tower (no reference count changes);
    other operations restore the flat form (_gr_tower_lazy_update, under the lock).
    The dense form is made by the operations on elements at the outermost
    level only (the implementations, under the lock, read the flat data
    directly), and only for elements that are short or dense enough
    (DENSE_SHORT_LENGTH, GR_TOWER_OPT_DENSE_FORM_SPARSITY): a sparse
    element such as a power of a root of unity in a field of large
    degree stays flat, where it has a few terms rather than a
    coefficient for every power of a below its degree. The result of an
    operation with a flat operand that stays flat is not made dense
    either (see LOCKED_DENSE): a running sum of sparse terms would
    otherwise change its form at every step.
*/

/* makes x (not read concurrently: a result) dense, keeping the storage
   of an existing dense form; the value is then to be set */
static void
_dense_fit(gr_tower_lazy_elem_t x)
{
    if (x->repr != LAZY_REPR_DENSE)
    {
        _gr_tower_lazy_clear_data(x);
        x->repr = LAZY_REPR_DENSE;
        fmpq_poly_init(&x->elem.dense.poly);
        x->elem.dense.nf = NULL;
        LAZY_POISON_FLAT(x);
    }
}

static const _gr_tower_dense_field_struct *
_dense_field(gr_tower_flat_struct * F, slong gid, const fmpz_poly_struct * M)
{
    _gr_tower_dense_field_struct * D;
    slong i;

    for (D = F->dense_fields; D != NULL; D = D->next)
        if (D->gid == gid && D->src == M && D->version == F->ideal_version)
            return D;
    for (D = F->dense_fields; D != NULL; D = D->next)
    {
        if (D->gid == gid && fmpz_poly_equal(&D->modulus, M))
        {
            D->src = M;
            D->version = F->ideal_version;
            return D;
        }
    }

    D = flint_malloc(sizeof(_gr_tower_dense_field_struct));
    D->gid = gid;
    D->degree = fmpz_poly_degree(M);
    fmpz_poly_init(&D->modulus);
    fmpz_poly_set(&D->modulus, M);
    D->ones = 1;
    for (i = 0; i < M->length; i++)
        if (!fmpz_is_one(M->coeffs + i))
            D->ones = 0;
    D->src = M;
    D->version = F->ideal_version;
    D->next = F->dense_fields;
    F->dense_fields = D;
    return D;
}

/* The first generator of the tower of x, when the elements of the
   subfield it generates can have a dense form: its modulus (monic,
   integral, univariate) and generator record, else NULL. */
static const fmpz_poly_struct *
_dense_generator(const gr_tower_gen_struct ** gp, const gr_tower_flat_struct * F)
{
    const gr_tower_gen_struct * g;

    if (F->T->num_gens == 0)
        return NULL;
    g = GR_TOWER_GEN(F->T, 0);
    if (g->kind != GR_TOWER_ALGEBRAIC || g->index != 1 || g->status != GR_TOWER_STATUS_PROVEN)
        return NULL;
    if (F->ideal_len < 1 || F->ideal_univar == NULL || F->ideal_univar[0] == NULL)
        return NULL;
    *gp = g;
    return F->ideal_univar[0];
}

/* Polynomials of at most this length are dense forms whatever their
   sparsity: an operation on them costs little in either form, and
   mixing forms costs conversions (measured: the DFT benchmark in
   Q(zeta_256), of degree 128, is fastest with every element dense). */
#define DENSE_SHORT_LENGTH 256

/* 1 if the numerator num (lex order) of a flat element is a polynomial
   in the variable v longer than DENSE_SHORT_LENGTH with fewer than its
   length / max_terms_ratio terms (the leading term has the highest power
   of v, which decides this before the terms are read) */
static int
_dense_too_sparse(const fmpz_mpoly_struct * num, const fmpz_mpoly_ctx_t mctx, slong v, slong max_terms_ratio)
{
    slong off, shift, len;
    ulong mask;

    if (max_terms_ratio <= 0 || num->length == 0 || num->bits > FLINT_BITS || mctx->minfo->ord != ORD_LEX)
        return 0;
    mpoly_gen_offset_shift_sp(&off, &shift, v, num->bits, mctx->minfo);
    mask = (-UWORD(1)) >> (FLINT_BITS - num->bits);
    len = ((num->exps[off] >> shift) & mask) + 1;
    return len > DENSE_SHORT_LENGTH && num->length * max_terms_ratio < len;
}

/* The flat element x (level <= 1, reduced, integer denominator) as a
   canonical polynomial in the first generator, given that it involves
   no other generator and has degree < d in it; returns 0 if x does not
   have this form, or (with max_terms_ratio > 0) if the polynomial is
   longer than DENSE_SHORT_LENGTH and has fewer than (its length) /
   max_terms_ratio nonzero coefficients. */
static int
_dense_poly_from_flat(fmpq_poly_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t mctx, slong v, slong d, slong max_terms_ratio)
{
    const fmpz_mpoly_struct * num = fmpz_mpoly_q_numref(x);
    const fmpz_mpoly_struct * den = fmpz_mpoly_q_denref(x);
    slong N, off, shift, t, len;
    ulong mask, e;

    if (!fmpz_mpoly_is_fmpz(den, mctx) || num->bits > FLINT_BITS || mctx->minfo->ord != ORD_LEX)
        return 0;

    N = mpoly_words_per_exp_sp(num->bits, mctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, v, num->bits, mctx->minfo);
    mask = (-UWORD(1)) >> (FLINT_BITS - num->bits);

    if (_dense_too_sparse(num, mctx, v, max_terms_ratio))
        return 0;

    /* every term must be a power of a alone (checked through the total
       of the other fields being zero) */
    len = 0;
    for (t = 0; t < num->length; t++)
    {
        slong k;
        e = (num->exps[N * t + off] >> shift) & mask;
        if (e >= (ulong) d)
            return 0;
        for (k = 0; k < N; k++)
        {
            ulong w = num->exps[N * t + k];
            if (k == off)
                w &= ~(mask << shift);
            if (w != 0)
                return 0;
        }
        len = FLINT_MAX(len, (slong) e + 1);
    }

    fmpq_poly_fit_length(res, len);
    _fmpz_vec_zero(res->coeffs, len);
    for (t = 0; t < num->length; t++)
    {
        e = (num->exps[N * t + off] >> shift) & mask;
        fmpz_set(res->coeffs + e, num->coeffs + t);
    }
    _fmpq_poly_set_length(res, len);
    fmpz_mpoly_get_fmpz(fmpq_poly_denref(res), den, mctx);
    if (fmpz_sgn(fmpq_poly_denref(res)) < 0)
    {
        fmpz_neg(fmpq_poly_denref(res), fmpq_poly_denref(res));
        _fmpz_vec_neg(res->coeffs, res->coeffs, len);
    }
    /* (the flat form is canonical: the content of the numerator is
       coprime to the integer denominator) */
    return 1;
}

/* 0 if x (flat) cannot take a dense form, by tests that do not bring x
   up to date or read all of its terms; 1 if it may */
static int
_gr_tower_lazy_dense_precheck(const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    gr_tower_flat_struct * F = x->F;

    if (!_gr_tower_lazy_dense_candidate(x) || F == LAZY(ctx)->trivial || _gr_tower_lazy_tls_depth != 1)
        return 0;
    if (x->structure_version == F->T->structure_version && x->elem.flat.mctx == F->mctx &&
        x->reduced_version == F->ideal_version && x->level == 1 && F->T->length >= 1 &&
        _dense_too_sparse(fmpz_mpoly_q_numref(&x->elem.flat.data), x->elem.flat.mctx, GR_TOWER_FLAT_VAR(F, 1),
            LAZY(ctx)->options[GR_TOWER_OPT_DENSE_FORM_SPARSITY]))
        return 0;
    return 1;
}

/* The dense form of x (flat, not shallow) in P, when it applies (under
   the lock, at the outermost level), with its descriptor; 0 if it does
   not apply. */
static int
_gr_tower_lazy_dense_extract(fmpq_poly_t P, const _gr_tower_dense_field_struct ** D, gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    gr_tower_flat_struct * F = x->F;
    const gr_tower_gen_struct * g;
    const fmpz_poly_struct * M;
    slong d;

    if (x->repr != LAZY_REPR_FLAT || x->shallow || F == LAZY(ctx)->trivial || _gr_tower_lazy_tls_depth != 1)
        return 0;

    /* (quick rejections before bringing x up to date) */
    if (F->T->num_gens == 0 || GR_TOWER_GEN(F->T, 0)->kind != GR_TOWER_ALGEBRAIC ||
        GR_TOWER_GEN(F->T, 0)->status != GR_TOWER_STATUS_PROVEN)
        return 0;
    if (x->structure_version == F->T->structure_version && x->elem.flat.mctx == F->mctx && x->level != 1)
        return 0;

    _gr_tower_lazy_update(x);
    if (x->level != 1 || x->reduced_version != F->ideal_version)
        return 0;
    M = _dense_generator(&g, F);
    if (M == NULL)
        return 0;
    d = fmpz_poly_degree(M);
    if (d < 2 || d > LAZY(ctx)->options[GR_TOWER_OPT_DENSE_FORM_DEGREE_LIMIT])
        return 0;

    if (!_dense_poly_from_flat(P, &x->elem.flat.data, x->elem.flat.mctx, GR_TOWER_FLAT_VAR(F, 1), d,
            LAZY(ctx)->options[GR_TOWER_OPT_DENSE_FORM_SPARSITY]))
        return 0;

    *D = _dense_field(F, g->gid, M);
    return 1;
}

/* x takes the dense form P (consumed); x may be an operand read
   concurrently by other threads: its flat data is not read by them --
   they take the lock for a flat operand -- and the dense data is
   written before the representation is published */
static void
_gr_tower_lazy_dense_commit(gr_tower_lazy_elem_t x, fmpq_poly_t P, const _gr_tower_dense_field_struct * D)
{
    fmpz_mpoly_q_clear(&x->elem.flat.data, x->elem.flat.mctx);
    LAZY_POISON_FLAT(x);
    x->elem.dense.poly = *P;
    x->elem.dense.nf = D;
    LAZY_REPR_SET(x, LAZY_REPR_DENSE);
}

/* the dense representation of x, when it applies */
void
_gr_tower_lazy_try_dense(gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    fmpq_poly_t P;
    const _gr_tower_dense_field_struct * D;

    fmpq_poly_init(P);
    if (_gr_tower_lazy_dense_extract(P, &D, x, ctx))
        _gr_tower_lazy_dense_commit(x, P, D);
    else
        fmpq_poly_clear(P);
}

/* After a successful locked operation res = op(x, y) at the outermost
   level (x, y may be NULL or alias res): the operands and the result
   take dense forms when every flat operand can; if one stays flat,
   nothing is converted (a running sum of sparse terms would otherwise
   be converted to the dense form, and read back through a flat view,
   at every step). */
void
_gr_tower_lazy_dense_after(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx)
{
    fmpq_poly_t P[2];
    const _gr_tower_dense_field_struct * D[2];
    gr_tower_lazy_elem_struct * op[2];
    slong i, n = 0;
    int ok = 1;

    if (x != NULL && x != res)
        op[n++] = (gr_tower_lazy_elem_struct *) x;
    if (y != NULL && y != res && y != x)
        op[n++] = (gr_tower_lazy_elem_struct *) y;

    /* (cheap rejections of all operands before any is converted) */
    for (i = 0; i < n; i++)
        if (LAZY_IS_FLAT(op[i]) && !_gr_tower_lazy_dense_precheck(op[i], ctx))
            return;
    if (n == 0 && !_gr_tower_lazy_dense_candidate(res))
        return;

    for (i = 0; i < n; i++)
        fmpq_poly_init(P[i]);

    for (i = 0; i < n && ok; i++)
    {
        if (!LAZY_IS_FLAT(op[i]))
            D[i] = NULL;     /* rational or dense already */
        else if (!_gr_tower_lazy_dense_extract(P[i], D + i, op[i], ctx))
            ok = 0;
    }

    for (i = 0; i < n; i++)
    {
        if (ok && LAZY_IS_FLAT(op[i]))
            _gr_tower_lazy_dense_commit(op[i], P[i], D[i]);
        else
            fmpq_poly_clear(P[i]);
    }

    if (ok && _gr_tower_lazy_dense_candidate(res))
        _gr_tower_lazy_try_dense(res, ctx);
}

/* the flat data of a dense x, in the current context of F */
void
_gr_tower_lazy_dense_get_flat(fmpz_mpoly_q_t res, const gr_tower_lazy_elem_t x, gr_tower_flat_struct * F)
{
    const _gr_tower_dense_field_struct * D = x->elem.dense.nf;
    const gr_tower_gen_struct * g;
    const fmpq_poly_struct * P = &x->elem.dense.poly;
    fmpz_mpoly_struct * num;
    slong d = D->degree, d0, v, N, off, shift, e, k, nnz;
    flint_bitcnt_t bits;

    gr_tower_flat_ensure(F);
    d0 = gr_tower_gid_order(F->T, D->gid);
    if (d0 < 0)
        flint_throw(FLINT_ERROR, "(%s): generator not found\n", __func__);
    v = GR_TOWER_FLAT_VAR_D(F, d0);

    num = fmpz_mpoly_q_numref(res);
    for (e = nnz = 0; e < P->length; e++)
        nnz += !fmpz_is_zero(P->coeffs + e);
    bits = mpoly_fix_bits(FLINT_MAX(MPOLY_MIN_BITS, 1 + FLINT_BIT_COUNT(d)), F->mctx->minfo);
    fmpz_mpoly_fit_length_reset_bits(num, nnz, bits, F->mctx);
    N = mpoly_words_per_exp_sp(bits, F->mctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, v, bits, F->mctx->minfo);
    for (e = P->length - 1, k = 0; e >= 0; e--)
    {
        if (fmpz_is_zero(P->coeffs + e))
            continue;
        mpoly_monomial_zero(num->exps + N * k, N);
        num->exps[N * k + off] = ((ulong) e) << shift;
        fmpz_set(num->coeffs + k, P->coeffs + e);
        k++;
    }
    _fmpz_mpoly_set_length(num, k, F->mctx);
    fmpz_mpoly_set_fmpz(fmpz_mpoly_q_denref(res), fmpq_poly_denref(P), F->mctx);

    /* (still reduced unless the generator has become of lower degree) */
    g = GR_TOWER_GEN(F->T, d0);
    if (!(g->kind == GR_TOWER_ALGEBRAIC && gr_tower_step_degree(F->T, g->index) == d))
        GR_MUST_SUCCEED(gr_tower_flat_reduce(res, F));
}

void
_gr_tower_lazy_dense_view_dense(gr_tower_lazy_dense_view_struct * V, const gr_tower_lazy_elem_t x)
{
    V->p = &x->elem.dense.poly;
}

/* (a shallow view of the rational number q, valid while q is unchanged) */
void
_gr_tower_lazy_dense_view_fmpq(gr_tower_lazy_dense_view_struct * V, const fmpq_t q)
{
    V->c0 = *fmpq_numref(q);
    V->tmp.coeffs = &V->c0;
    V->tmp.alloc = 1;
    V->tmp.length = !fmpz_is_zero(&V->c0);
    V->tmp.den[0] = *fmpq_denref(q);
    V->p = &V->tmp;
}

/* whether x has a view over the descriptor D of F */
int
_gr_tower_lazy_dense_view(gr_tower_lazy_dense_view_struct * V, const gr_tower_lazy_elem_t x, const _gr_tower_dense_field_struct * D, const gr_tower_flat_struct * F, gr_ctx_t ctx)
{
    if (LAZY_HAS_DENSE(x))
    {
        if (x->elem.dense.nf != D || x->F != F)
            return 0;
        _gr_tower_lazy_dense_view_dense(V, x);
        return 1;
    }
    if (LAZY_IS_RATIONAL(x, LAZY(ctx)))
    {
        _gr_tower_lazy_dense_view_fmpq(V, &x->elem.q);
        return 1;
    }
    return 0;
}

/* The numerator C (length n, with room for n coefficients) reduced
   modulo the monic modulus of D: the length becomes < d (C is not
   normalized). */
static void
_dense_reduce(fmpz * C, slong n, const _gr_tower_dense_field_struct * D)
{
    slong d = D->degree, i, j;
    const fmpz * M = D->modulus.coeffs;

    if (n <= d)
        return;

    if (D->ones)
    {
        /* x^(d+1) = 1, then x^d = -(1 + ... + x^(d-1)) */
        for (i = n - 1; i >= d + 1; i--)
        {
            fmpz_add(C + i - d - 1, C + i - d - 1, C + i);
            fmpz_zero(C + i);
        }
        if (!fmpz_is_zero(C + d))
        {
            for (i = 0; i < d; i++)
                fmpz_sub(C + i, C + i, C + d);
            fmpz_zero(C + d);
        }
        return;
    }

    if (d >= 64)
    {
        /* a large modulus: fast division unless it is sparse (as in
           _mul_dense_univar) */
        slong nnz = 0;
        for (j = 0; j < d; j++)
            nnz += !fmpz_is_zero(M + j);
        if (nnz > 16)
        {
            fmpz * R = _fmpz_vec_init(n);
            _fmpz_poly_rem(R, C, n, M, d + 1);
            _fmpz_vec_swap(C, R, d);
            _fmpz_vec_zero(C + d, n - d);
            _fmpz_vec_clear(R, n);
            return;
        }
    }

    for (i = n - 1; i >= d; i--)
    {
        if (fmpz_is_zero(C + i))
            continue;
        for (j = 0; j < d; j++)
            if (!fmpz_is_zero(M + j))
                fmpz_submul(C + i - d + j, C + i, M + j);
        fmpz_zero(C + i);
    }
}

/* the exponent of the single term of P, or -1 */
static slong
_dense_monomial(const fmpq_poly_struct * P)
{
    slong i;
    for (i = 0; i < P->length - 1; i++)
        if (!fmpz_is_zero(P->coeffs + i))
            return -1;
    return P->length - 1;
}

/* R = X * Y mod M (R canonical; R may alias neither X nor Y) */
static void
_dense_mul(fmpq_poly_t R, const fmpq_poly_struct * X, const fmpq_poly_struct * Y, const _gr_tower_dense_field_struct * D)
{
    slong e, n;

    if (X->length == 0 || Y->length == 0)
    {
        fmpq_poly_zero(R);
        return;
    }

    if (X->length < Y->length)
        FLINT_SWAP(const fmpq_poly_struct *, X, Y);

    if (Y->length == 1)
    {
        fmpq_t c;
        *fmpq_numref(c) = Y->coeffs[0];
        *fmpq_denref(c) = Y->den[0];
        fmpq_poly_scalar_mul_fmpq(R, X, c);
        return;
    }

    n = X->length + Y->length - 1;
    fmpq_poly_fit_length(R, n);

    if ((e = _dense_monomial(Y)) >= 0 || (e = _dense_monomial(X)) >= 0)
    {
        /* a monomial factor (a power of a root of unity, say): a shift */
        const fmpq_poly_struct * P = (_dense_monomial(Y) == e) ? X : Y;
        const fmpq_poly_struct * Q = (P == X) ? Y : X;
        _fmpz_vec_zero(R->coeffs, e);
        _fmpz_vec_scalar_mul_fmpz(R->coeffs + e, P->coeffs, P->length, Q->coeffs + e);
        _fmpz_vec_zero(R->coeffs + e + P->length, n - e - P->length);
    }
    else
    {
        _fmpz_poly_mul(R->coeffs, X->coeffs, X->length, Y->coeffs, Y->length);
    }

    fmpz_mul(fmpq_poly_denref(R), X->den, Y->den);
    _dense_reduce(R->coeffs, n, D);
    _fmpq_poly_set_length(R, FLINT_MIN(n, D->degree));
    _fmpq_poly_normalise(R);
    fmpq_poly_canonicalise(R);
}

/* R = 1 / Y, Y nonzero; returns 0 if the modulus turns out to be
   reducible (not for proven steps) */
static int
_dense_inv(fmpq_poly_t R, const fmpq_poly_struct * Y, const _gr_tower_dense_field_struct * D)
{
    slong d = D->degree;

    if (Y->length == 1)
    {
        fmpq_poly_fit_length(R, 1);
        fmpz_set(R->coeffs, Y->den);
        fmpz_set(fmpq_poly_denref(R), Y->coeffs);
        _fmpq_poly_set_length(R, 1);
        fmpq_poly_canonicalise(R);
        return 1;
    }

    if (d == 2)
    {
        /* 1/(u + v x) = (u - p v - v x) / (u^2 - p u v + q v^2) for the
           modulus x^2 + p x + q; times the denominator of Y */
        const fmpz * u = Y->coeffs, * w = Y->coeffs + 1, * q = D->modulus.coeffs, * p = D->modulus.coeffs + 1;
        fmpz_t a, b;
        fmpz_init(a);
        fmpz_init(b);
        fmpz_mul(a, u, u);
        fmpz_mul(b, u, w);
        fmpz_submul(a, b, p);
        fmpz_mul(b, w, w);
        fmpz_addmul(a, b, q);              /* the norm N */
        fmpq_poly_fit_length(R, 2);
        fmpz_mul(R->coeffs + 0, p, w);
        fmpz_sub(R->coeffs + 0, u, R->coeffs + 0);
        fmpz_neg(R->coeffs + 1, w);
        _fmpz_vec_scalar_mul_fmpz(R->coeffs, R->coeffs, 2, Y->den);
        fmpz_swap(fmpq_poly_denref(R), a);
        _fmpq_poly_set_length(R, 2);
        _fmpq_poly_normalise(R);
        fmpq_poly_canonicalise(R);
        fmpz_clear(a);
        fmpz_clear(b);
        return 1;
    }
    else
    {
        /* 1/g = t/r for s M + t g = r (the resultant), g = Y / content */
        fmpz * g, * sv;
        fmpz_t r, cont;
        int ok;
        g = _fmpz_vec_init(Y->length);
        sv = _fmpz_vec_init(Y->length);
        fmpz_init(r);
        fmpz_init(cont);
        fmpq_poly_fit_length(R, d + 1);
        _fmpz_vec_zero(R->coeffs, d + 1);
        _fmpz_vec_content(cont, Y->coeffs, Y->length);
        _fmpz_vec_scalar_divexact_fmpz(g, Y->coeffs, Y->length, cont);
        _fmpz_poly_xgcd(r, sv, R->coeffs, D->modulus.coeffs, d + 1, g, Y->length);
        ok = !fmpz_is_zero(r);
        if (ok)
        {
            _fmpz_vec_scalar_mul_fmpz(R->coeffs, R->coeffs, d, Y->den);
            fmpz_mul(fmpq_poly_denref(R), r, cont);
            _fmpq_poly_set_length(R, d);
            _fmpq_poly_normalise(R);
            fmpq_poly_canonicalise(R);
        }
        _fmpz_vec_clear(g, Y->length);
        _fmpz_vec_clear(sv, Y->length);
        fmpz_clear(r);
        fmpz_clear(cont);
        return ok;
    }
}

/* A per-thread scratch polynomial for results: it is swapped with the
   storage of the result, so that operations in steady state allocate
   nothing (freed by flint_cleanup in the thread). */
static FLINT_TLS_PREFIX fmpq_poly_struct _dense_scratch;
static FLINT_TLS_PREFIX int _dense_scratch_ready = 0;

static void
_dense_scratch_cleanup(void)
{
    if (_dense_scratch_ready)
        fmpq_poly_clear(&_dense_scratch);
    _dense_scratch_ready = 0;
}

static fmpq_poly_struct *
_dense_scratch_get(void)
{
    if (!_dense_scratch_ready)
    {
        fmpq_poly_init(&_dense_scratch);
        flint_register_cleanup_function(_dense_scratch_cleanup);
        _dense_scratch_ready = 1;
    }
    return &_dense_scratch;
}

/* res = X op Y over D (res lives in the tower of D, not shallow);
   returns 0 if not handled (division by zero) */
int
_gr_tower_lazy_dense_op(gr_tower_lazy_elem_t res, const gr_tower_lazy_dense_view_struct * X, const gr_tower_lazy_dense_view_struct * Y, int op, const _gr_tower_dense_field_struct * D)
{
    fmpq_poly_struct * R;
    fmpq_poly_t I;
    const fmpq_poly_struct * x = X->p, * y = Y->p;

    if (op == DENSE_DIV && y->length == 0)
        return 0;

    if ((op == DENSE_ADD || op == DENSE_SUB) && res->repr == LAZY_REPR_DENSE)
    {
        /* in place in the storage of res, which (being dense) is an
           operand only through its own poly, an alias fmpq_poly allows */
        if (op == DENSE_ADD)
            fmpq_poly_add(&res->elem.dense.poly, x, y);
        else
            fmpq_poly_sub(&res->elem.dense.poly, x, y);
        res->elem.dense.nf = D;
        res->level = (res->elem.dense.poly.length > 1);
        return 1;
    }

    R = _dense_scratch_get();
    if (op == DENSE_ADD)
        fmpq_poly_add(R, x, y);
    else if (op == DENSE_SUB)
        fmpq_poly_sub(R, x, y);
    else if (op == DENSE_MUL)
        _dense_mul(R, x, y, D);
    else
    {
        fmpq_poly_init(I);
        if (!_dense_inv(I, y, D))
        {
            fmpq_poly_clear(I);
            flint_throw(FLINT_ERROR, "(%s): reducible modulus\n", __func__);
        }
        _dense_mul(R, x, I, D);
        fmpq_poly_clear(I);
    }

    _dense_fit(res);
    fmpq_poly_swap(&res->elem.dense.poly, R);
    res->elem.dense.nf = D;
    /* (the prefix: the generator a, unless the value is a constant; the
       versions are brought up to date when the element is next read
       under the lock) */
    res->level = (res->elem.dense.poly.length > 1);
    return 1;
}

/* lock-free binary operation with a dense operand: returns 0 if it does
   not apply */
int
_gr_tower_lazy_dense_binary(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, int op, gr_ctx_t ctx)
{
    const gr_tower_lazy_elem_struct * D = LAZY_HAS_DENSE(x) ? x : (LAZY_HAS_DENSE(y) ? y : NULL);
    gr_tower_lazy_dense_view_struct X, Y;

    if (D == NULL || _gr_tower_lazy_tls_depth != 0 || res->shallow || res->F != D->F)
        return 0;
    if (!_gr_tower_lazy_dense_view(&X, x, D->elem.dense.nf, D->F, ctx) || !_gr_tower_lazy_dense_view(&Y, y, D->elem.dense.nf, D->F, ctx))
        return 0;
    return _gr_tower_lazy_dense_op(res, &X, &Y, op, D->elem.dense.nf);
}

/* with a rational scalar */
int
_gr_tower_lazy_dense_scalar(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const fmpq_t c, int op, gr_ctx_t ctx)
{
    gr_tower_lazy_dense_view_struct X, Y;

    if (!LAZY_HAS_DENSE(x) || _gr_tower_lazy_tls_depth != 0 || res->shallow || res->F != x->F)
        return 0;
    _gr_tower_lazy_dense_view_dense(&X, x);
    _gr_tower_lazy_dense_view_fmpq(&Y, c);
    return _gr_tower_lazy_dense_op(res, &X, &Y, op, x->elem.dense.nf);
}

/* (a cheap test before _gr_tower_lazy_try_dense: flat, first generator a proven
   algebraic one, and not known to be at another level) */
int
_gr_tower_lazy_dense_candidate(const gr_tower_lazy_elem_t x)
{
    const gr_tower_struct * T;
    if (x == NULL || x->repr != LAZY_REPR_FLAT || x->shallow)
        return 0;
    T = x->F->T;
    return T->num_gens != 0 && T->gens[0].kind == GR_TOWER_ALGEBRAIC &&
           T->gens[0].status == GR_TOWER_STATUS_PROVEN &&
           (x->level == 1 || x->structure_version != T->structure_version || x->elem.flat.mctx != x->F->mctx);
}

/* x as a polynomial in the first generator a of its tower, with the
   modulus of a (see gr_tower_lazy_get_fmpq_poly); under the lock, x
   being owned by the caller (it may be reduced in place) */
int
_gr_tower_lazy_get_fmpq_poly_impl(fmpq_poly_t res, fmpz_poly_t modulus, gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    gr_tower_flat_struct * F;
    const gr_tower_gen_struct * g;
    const fmpz_poly_struct * M;

    if (LAZY_REPR(x) == LAZY_REPR_RATIONAL)
    {
        fmpq_poly_set_fmpq(res, &x->elem.q);
        if (modulus != NULL)
        {
            fmpz_poly_zero(modulus);
            fmpz_poly_set_coeff_ui(modulus, 1, 1);
        }
        return GR_SUCCESS;
    }

    if (LAZY_REPR(x) == LAZY_REPR_DENSE)
    {
        fmpq_poly_set(res, &x->elem.dense.poly);
        if (modulus != NULL)
            fmpz_poly_set(modulus, &x->elem.dense.nf->modulus);
        return GR_SUCCESS;
    }

    _gr_tower_lazy_update(x);
    F = x->F;
    if (x->reduced_version != F->ideal_version)
    {
        GR_MUST_SUCCEED(gr_tower_flat_reduce(&x->elem.flat.data, F));
        x->reduced_version = F->ideal_version;
        _gr_tower_lazy_shrink(x);
    }

    if (x->level == 0)
    {
        fmpq_t c;
        int status;
        fmpq_init(c);
        status = fmpz_mpoly_q_get_fmpq(c, &x->elem.flat.data, x->elem.flat.mctx) ? GR_SUCCESS : GR_UNABLE;
        if (status == GR_SUCCESS)
        {
            fmpq_poly_set_fmpq(res, c);
            if (modulus != NULL)
            {
                fmpz_poly_zero(modulus);
                fmpz_poly_set_coeff_ui(modulus, 1, 1);
            }
        }
        fmpq_clear(c);
        return status;
    }

    /* (x involves another generator, or the first is not algebraic) */
    if (x->level != 1 || F->T->num_gens == 0)
        return GR_DOMAIN;
    g = GR_TOWER_GEN(F->T, 0);
    if (g->kind != GR_TOWER_ALGEBRAIC || g->index != 1)
        return GR_DOMAIN;
    /* (a modulus which is not monic integral, or an algebraic denominator) */
    if (F->ideal_len < 1 || F->ideal_univar == NULL || F->ideal_univar[0] == NULL)
        return GR_UNABLE;
    M = F->ideal_univar[0];

    if (!_dense_poly_from_flat(res, &x->elem.flat.data, x->elem.flat.mctx, GR_TOWER_FLAT_VAR(F, 1), fmpz_poly_degree(M), 0))
        return GR_UNABLE;
    if (modulus != NULL)
        fmpz_poly_set(modulus, M);
    return GR_SUCCESS;
}
