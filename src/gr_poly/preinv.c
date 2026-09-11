/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_generic.h"
#include "gr_poly.h"

#if FLINT_USES_PTHREAD
#include <pthread.h>
#endif

/*
    A gr_poly_preinv_t holds precomputed data for reduction modulo a fixed
    polynomial f, in one of several representations:

    - GR_POLY_PREINV_NEWTON: the inverse of the reversal of f as a power
      series to length len(f) (Newton division: two multiplications).
    - GR_POLY_PREINV_SPARSE: the exponents and (negated, normalised)
      coefficients of the few nonzero terms of f below the leading term
      (schoolbook reduction costing nz operations per reduced coefficient).
    - GR_POLY_PREINV_TRANSFORMED (todo): Newton inverse with precomputed
      transforms of f and its inverse for long dense moduli, building on
      the gr_transformed_poly representations.

    The division method _gr_poly_preinv_divrem is ring-overloadable
    (GR_METHOD_POLY_DIVREM_PREINV) so that rings can supply fast inner
    loops for the sparse representation; the generic implementation
    handles all representations.
*/

void
gr_poly_preinv_init(gr_poly_preinv_t P, gr_ctx_t FLINT_UNUSED(ctx))
{
    P->kind = GR_POLY_PREINV_NEWTON;
    P->owned = 0;
    P->monic = 0;
    P->lenf = 0;
    P->f = NULL;
    P->finv = NULL;
    P->lenfinv = 0;
    P->nz = 0;
    P->exps = NULL;
    P->coeffs = NULL;
    P->lcinv = NULL;
    P->tctx = NULL;
    P->tfinv = NULL;
    P->cctx = NULL;
    P->L = 0;
    P->tf = NULL;
    P->scratch = NULL;
}

/*
    Scratch storage for the transformed representation. Transformed
    elements are large (a full set of transforms) and allocating them
    for each reduction costs more than the transforms themselves, so
    they are kept in the object. To keep the object usable from several
    threads at once (for example inside multithreaded factorization),
    the scratch is a pool of entries, each used by one thread at a time
    and handed out under a lock; the pool grows when all entries are in
    use. The pool is logically part of the object's mutable cache.
*/
typedef struct
{
    gr_ptr t[2];        /* elements of tctx */
    gr_ptr c[2];        /* elements of cctx */
    gr_ptr coeffs;      /* 2 lenf base ring elements */
    int in_use;
}
preinv_scratch_entry;

typedef struct
{
    preinv_scratch_entry ** entries;    /* separately allocated, so that
                                           growing the pool does not move
                                           entries in use by other threads */
    slong num;
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
#endif
}
preinv_scratch_pool;

static preinv_scratch_pool *
_preinv_pool(const gr_poly_preinv_t P)
{
    return (preinv_scratch_pool *) P->scratch;
}

static void
_preinv_pool_init(gr_poly_preinv_t P)
{
    preinv_scratch_pool * pool = flint_malloc(sizeof(preinv_scratch_pool));
    pool->entries = NULL;
    pool->num = 0;
#if FLINT_USES_PTHREAD
    pthread_mutex_init(&pool->mutex, NULL);
#endif
    P->scratch = pool;
}

static void
_preinv_pool_clear(gr_poly_preinv_t P, gr_ctx_t ctx)
{
    preinv_scratch_pool * pool = _preinv_pool(P);
    slong i;

    if (pool == NULL)
        return;

    for (i = 0; i < pool->num; i++)
    {
        gr_heap_clear(pool->entries[i]->t[0], P->tctx);
        gr_heap_clear(pool->entries[i]->t[1], P->tctx);
        gr_heap_clear(pool->entries[i]->c[0], P->cctx);
        gr_heap_clear(pool->entries[i]->c[1], P->cctx);
        gr_heap_clear_vec(pool->entries[i]->coeffs, 2 * P->lenf, ctx);
        flint_free(pool->entries[i]);
    }

    flint_free(pool->entries);
#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&pool->mutex);
#endif
    flint_free(pool);
    P->scratch = NULL;
}

/* Returns an entry for exclusive use by the caller. */
static preinv_scratch_entry *
_preinv_pool_acquire(const gr_poly_preinv_t P, gr_ctx_t ctx)
{
    preinv_scratch_pool * pool = _preinv_pool(P);
    preinv_scratch_entry * e;
    slong i;

#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&pool->mutex);
#endif

    for (i = 0; i < pool->num; i++)
        if (!pool->entries[i]->in_use)
            break;

    if (i == pool->num)
    {
        e = flint_malloc(sizeof(preinv_scratch_entry));
        e->t[0] = gr_heap_init(P->tctx);
        e->t[1] = gr_heap_init(P->tctx);
        e->c[0] = gr_heap_init(P->cctx);
        e->c[1] = gr_heap_init(P->cctx);
        e->coeffs = gr_heap_init_vec(2 * P->lenf, ctx);
        e->in_use = 0;
        pool->entries = flint_realloc(pool->entries, (pool->num + 1) * sizeof(preinv_scratch_entry *));
        pool->entries[pool->num] = e;
        pool->num++;
    }

    e = pool->entries[i];
    e->in_use = 1;

#if FLINT_USES_PTHREAD
    pthread_mutex_unlock(&pool->mutex);
#endif

    return e;
}

static void
_preinv_pool_release(const gr_poly_preinv_t P, preinv_scratch_entry * e)
{
#if FLINT_USES_PTHREAD
    preinv_scratch_pool * pool = _preinv_pool(P);
    pthread_mutex_lock(&pool->mutex);
    e->in_use = 0;
    pthread_mutex_unlock(&pool->mutex);
#else
    e->in_use = 0;
#endif
}

static void
_gr_poly_preinv_free(gr_poly_preinv_t P, gr_ctx_t ctx)
{
    if (P->owned)
    {
        if (P->f != NULL)
            gr_heap_clear_vec((gr_ptr) P->f, P->lenf, ctx);
        if (P->finv != NULL)
            gr_heap_clear_vec(P->finv, P->lenfinv, ctx);
    }

    if (P->exps != NULL)
        flint_free(P->exps);
    if (P->coeffs != NULL)
        gr_heap_clear_vec(P->coeffs, P->nz, ctx);
    if (P->lcinv != NULL)
        gr_heap_clear(P->lcinv, ctx);

    _preinv_pool_clear(P, ctx);

    if (P->tctx != NULL)
    {
        if (P->tfinv != NULL)
            gr_heap_clear(P->tfinv, P->tctx);
        gr_ctx_clear(P->tctx);
        flint_free(P->tctx);
    }

    if (P->cctx != NULL)
    {
        if (P->tf != NULL)
            gr_heap_clear(P->tf, P->cctx);
        gr_ctx_clear(P->cctx);
        flint_free(P->cctx);
    }

    gr_poly_preinv_init(P, ctx);
}

void
gr_poly_preinv_clear(gr_poly_preinv_t P, gr_ctx_t ctx)
{
    _gr_poly_preinv_free(P, ctx);
}

void
_gr_poly_preinv_init_newton_shallow(gr_poly_preinv_t P, gr_srcptr f, slong lenf,
    gr_srcptr finv, slong lenfinv, gr_ctx_t ctx)
{
    gr_poly_preinv_init(P, ctx);
    P->kind = GR_POLY_PREINV_NEWTON;
    P->owned = 0;
    P->lenf = lenf;
    P->f = f;
    P->finv = (gr_ptr) finv;
    P->lenfinv = lenfinv;
    P->monic = (lenf >= 1 && gr_is_one(GR_ENTRY(f, lenf - 1, ctx->sizeof_elem), ctx) == T_TRUE);
}

static int
_gr_poly_preinv_set_f(gr_poly_preinv_t P, gr_srcptr f, slong lenf, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    _gr_poly_preinv_free(P, ctx);

    if (lenf == 0)
        return GR_DOMAIN;

    P->owned = 1;
    P->lenf = lenf;
    P->f = gr_heap_init_vec(lenf, ctx);
    status |= _gr_vec_set((gr_ptr) P->f, f, lenf, ctx);
    P->monic = (gr_is_one(GR_ENTRY(f, lenf - 1, ctx->sizeof_elem), ctx) == T_TRUE);

    return status;
}

int
_gr_poly_preinv_set_newton(gr_poly_preinv_t P, gr_srcptr f, slong lenf, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    gr_ptr t;

    status |= _gr_poly_preinv_set_f(P, f, lenf, ctx);
    if (status != GR_SUCCESS)
        return status;

    P->kind = GR_POLY_PREINV_NEWTON;

    P->finv = gr_heap_init_vec(lenf, ctx);
    P->lenfinv = lenf;

    GR_TMP_INIT_VEC(t, lenf, ctx);
    status |= _gr_poly_reverse(t, f, lenf, lenf, ctx);
    status |= _gr_poly_inv_series(P->finv, t, lenf, lenf, ctx);
    GR_TMP_CLEAR_VEC(t, lenf, ctx);

    return status;
}

/* The Newton representation without an inverse: ordinary division is
   used. Preferable for short moduli over rings with cheap elements. */
int
_gr_poly_preinv_set_plain(gr_poly_preinv_t P, gr_srcptr f, slong lenf, gr_ctx_t ctx)
{
    int status = _gr_poly_preinv_set_f(P, f, lenf, ctx);
    P->kind = GR_POLY_PREINV_NEWTON;
    return status;
}

int
_gr_poly_preinv_set_sparse(gr_poly_preinv_t P, gr_srcptr f, slong lenf, gr_ctx_t ctx)
{
    slong i, nz, sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    truth_t z;

    if (lenf < 2)
        return GR_DOMAIN;

    /* count the nonzero terms; entries that cannot be verified zero are
       treated as nonzero */
    nz = 0;
    for (i = 0; i < lenf - 1; i++)
        nz += (gr_is_zero(GR_ENTRY(f, i, sz), ctx) != T_TRUE);

    status |= _gr_poly_preinv_set_f(P, f, lenf, ctx);
    if (status != GR_SUCCESS)
        return status;

    P->kind = GR_POLY_PREINV_SPARSE;
    P->nz = nz;
    P->exps = flint_malloc(sizeof(slong) * FLINT_MAX(nz, 1));
    P->coeffs = gr_heap_init_vec(FLINT_MAX(nz, 1), ctx);

    if (!P->monic)
    {
        P->lcinv = gr_heap_init(ctx);
        status |= gr_inv(P->lcinv, GR_ENTRY(f, lenf - 1, sz), ctx);
        if (status != GR_SUCCESS)
            return status;
    }

    nz = 0;
    for (i = 0; i < lenf - 1; i++)
    {
        z = gr_is_zero(GR_ENTRY(f, i, sz), ctx);
        if (z != T_TRUE)
        {
            P->exps[nz] = i;
            status |= gr_neg(GR_ENTRY(P->coeffs, nz, sz), GR_ENTRY(f, i, sz), ctx);
            if (!P->monic)
                status |= gr_mul(GR_ENTRY(P->coeffs, nz, sz), GR_ENTRY(P->coeffs, nz, sz), P->lcinv, ctx);
            nz++;
        }
    }

    return status;
}

/*
    Transformed representation for long dense moduli (as in NTL):

    - the quotient Q is computed as the window [lenf - 1, lenA) of
      A_hi * rev(finv) (a middle product) in a linear ring of transformed
      polynomials, with the transform of rev(finv) precomputed (one
      forward and one inverse transform of length ~ 2n instead of a full
      multiplication);

    - the remainder R = A - Q f is computed modulo x^L - 1 with L >= n,
      which loses nothing since deg R < n: the product Q f is formed in
      a cyclic ring of transformed polynomials of length L ~ n with the
      transform of f (mod x^L - 1) precomputed, so that the second
      multiplication costs one forward and one inverse transform of
      half the length.

    Altogether a reduction costs about one multiplication instead of
    the two of Newton division. The Newton inverse is also kept, for
    dividends longer than the precomputed capacity.
*/
int
_gr_poly_preinv_set_transformed(gr_poly_preinv_t P, gr_srcptr f, slong lenf, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    gr_transformed_poly_workload_t wl = {{ 4, 2, 2, 6, 0, 1 }};
    gr_ptr ff;
    slong i, n = lenf - 1, sz = ctx->sizeof_elem;

    status |= _gr_poly_preinv_set_newton(P, f, lenf, ctx);
    if (status != GR_SUCCESS)
        return status;

    if (P->lenfinv == 0)
        return GR_UNABLE;

    /* linear ring: products of length up to 2 lenf - 2 = 2n (the
       quotient product has length at most lenQ + lenf - 1 <= 2n, a
       product of two reduced polynomials at most 2n - 1) with at most
       lenf accumulated elementary products per coefficient; the bound
       is kept tight so that transform sizes are not rounded up
       needlessly (e.g. for n a power of two) */
    P->tctx = flint_malloc(sizeof(gr_ctx_struct));
    if (gr_ctx_init_gr_poly_transformed_repr(P->tctx, ctx, 2 * lenf - 2, lenf, wl) != GR_SUCCESS)
    {
        flint_free(P->tctx);
        P->tctx = NULL;
        return GR_UNABLE;
    }

    /* cyclic ring of length L >= n */
    P->L = n;
    P->cctx = flint_malloc(sizeof(gr_ctx_struct));
    if (gr_ctx_init_gr_poly_transformed_cyclic_repr(P->cctx, ctx, &P->L, lenf, wl) != GR_SUCCESS)
    {
        flint_free(P->cctx);
        P->cctx = NULL;
        gr_ctx_clear(P->tctx);
        flint_free(P->tctx);
        P->tctx = NULL;
        return GR_UNABLE;
    }

    /* transform of the reversed inverse */
    P->tfinv = gr_heap_init(P->tctx);
    GR_TMP_INIT_VEC(ff, P->lenfinv, ctx);
    status |= _gr_poly_reverse(ff, P->finv, P->lenfinv, P->lenfinv, ctx);
    status |= _gr_set_gr_poly(P->tfinv, ff, P->lenfinv, ctx, P->tctx);
    GR_TMP_CLEAR_VEC(ff, P->lenfinv, ctx);


    _preinv_pool_init(P);

    /* f mod x^L - 1 */
    GR_TMP_INIT_VEC(ff, P->L, ctx);
    status |= _gr_vec_set(ff, f, FLINT_MIN(lenf, P->L), ctx);
    for (i = P->L; i < lenf; i++)
        status |= gr_add(GR_ENTRY(ff, i - P->L, sz), GR_ENTRY(ff, i - P->L, sz), GR_ENTRY(f, i, sz), ctx);
    P->tf = gr_heap_init(P->cctx);
    status |= _gr_set_gr_poly(P->tf, ff, P->L, ctx, P->cctx);
    GR_TMP_CLEAR_VEC(ff, P->L, ctx);

    if (status == GR_SUCCESS)
        P->kind = GR_POLY_PREINV_TRANSFORMED;
    else
        status = GR_UNABLE;

    return status;
}

/*
    Division of A of length lenA <= 2 lenf - 2 using the precomputed
    transforms; R gets lenf - 1 coefficients, Q (if not NULL) the
    lenA - lenf + 1 quotient coefficients.
*/
static int
_gr_poly_preinv_divrem_transformed(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA,
    const gr_poly_preinv_t P, gr_ctx_t ctx)
{
    slong lenf = P->lenf, n = lenf - 1, lenQ = lenA - lenf + 1, L = P->L;
    slong i, sz = ctx->sizeof_elem;
    gr_ptr x, y, xc, yc, q, u;
    preinv_scratch_entry * e;
    int status = GR_SUCCESS;

    FLINT_ASSERT(lenA <= 2 * lenf - 2);

    e = _preinv_pool_acquire(P, ctx);
    x = e->t[0];
    y = e->t[1];
    xc = e->c[0];
    yc = e->c[1];
    q = e->coeffs;
    u = GR_ENTRY(q, lenf, sz);

    /* Q = window [lenf - 1, lenA) of A_hi * rev(finv) (middle product) */
    status |= _gr_set_gr_poly(x, GR_ENTRY(A, lenf - 1, sz), lenQ, ctx, P->tctx);
    status |= gr_mul(y, x, P->tfinv, P->tctx);
    status |= _gr_get_gr_poly_window_destructive(q, y, lenf - 1, lenA, ctx, P->tctx);

    /* R = (A - Q f) mod (x^L - 1), of which only the first n coefficients
       are nonzero */
    status |= _gr_set_gr_poly(xc, q, lenQ, ctx, P->cctx);
    status |= gr_mul(yc, xc, P->tf, P->cctx);
    status |= _gr_get_gr_poly_window_destructive(u, yc, 0, n, ctx, P->cctx);

    /* A mod x^L - 1, first n coefficients */
    status |= _gr_vec_sub(R, A, u, n, ctx);
    for (i = L; i < lenA; i += L)
        status |= _gr_vec_add(R, R, GR_ENTRY(A, i, sz), FLINT_MIN(lenA - i, n), ctx);

    if (Q != NULL)
        status |= _gr_vec_set(Q, q, lenQ, ctx);

    _preinv_pool_release(P, e);

    return status;
}

/* todo: tuning */
#define GR_POLY_PREINV_PLAIN_CUTOFF 8

/*
    Default selection of the representation. Rings can overload this
    (GR_METHOD_POLY_PREINV_SET) to take the degree, the sparsity and the
    size of the coefficients into account. The default uses the sparse
    representation when f has at most 8 nonzero terms below the leading
    term (so that a reduction costs about 8 multiplications per
    coefficient, cheaper than Newton division for any multiplication
    algorithm of practical interest), ordinary division for moduli of
    length at most 8, and otherwise the Newton representation.
*/
int
gr_generic_poly_preinv_set(gr_poly_preinv_struct * P, gr_srcptr f, slong lenf, gr_ctx_t ctx)
{
    slong i, nz, sz = ctx->sizeof_elem;
    const slong cutoff = 8;

    /* For short moduli, ordinary division is at least as fast as the
       two multiplications of Newton division and the inverse need not be
       computed; rings with cheap elements (nmod, mpn_mod) use larger
       cutoffs in their own selection. */
    if (lenf <= GR_POLY_PREINV_PLAIN_CUTOFF)
        return _gr_poly_preinv_set_plain(P, f, lenf, ctx);

    nz = 0;
    for (i = 0; i < lenf - 1 && nz <= cutoff; i++)
        nz += (gr_is_zero(GR_ENTRY(f, i, sz), ctx) != T_TRUE);

    if (nz <= cutoff)
        return _gr_poly_preinv_set_sparse(P, f, lenf, ctx);
    else
        return _gr_poly_preinv_set_newton(P, f, lenf, ctx);
}

int
gr_poly_preinv_set_newton(gr_poly_preinv_t P, const gr_poly_t f, gr_ctx_t ctx)
{
    return _gr_poly_preinv_set_newton(P, f->coeffs, f->length, ctx);
}

int
gr_poly_preinv_set_sparse(gr_poly_preinv_t P, const gr_poly_t f, gr_ctx_t ctx)
{
    return _gr_poly_preinv_set_sparse(P, f->coeffs, f->length, ctx);
}

int
gr_poly_preinv_set_transformed(gr_poly_preinv_t P, const gr_poly_t f, gr_ctx_t ctx)
{
    return _gr_poly_preinv_set_transformed(P, f->coeffs, f->length, ctx);
}

int
gr_poly_preinv_set(gr_poly_preinv_t P, const gr_poly_t f, gr_ctx_t ctx)
{
    return _gr_poly_preinv_set(P, f->coeffs, f->length, ctx);
}

/*
    Generic schoolbook reduction by a sparse modulus: with
    f = lc x^n + sum_k f_{e_k} x^{e_k}, we have x^n = sum_k c_k x^{e_k}
    where c_k = -f_{e_k} / lc, so each leading coefficient r_i of the
    remainder is eliminated by adding r_i c_k to the coefficient of
    x^{e_k + i - n} for each k (the quotient coefficient being r_i / lc).
    Costs nz multiplications per eliminated coefficient. Q may be NULL.
*/
static int
_gr_poly_preinv_divrem_sparse_generic(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA,
    const gr_poly_preinv_t P, gr_ctx_t ctx)
{
    slong i, k, n = P->lenf - 1, nz = P->nz, sz = ctx->sizeof_elem;
    gr_ptr r, c;
    int status = GR_SUCCESS;

    GR_TMP_INIT_VEC(r, lenA, ctx);
    GR_TMP_INIT(c, ctx);

    status |= _gr_vec_set(r, A, lenA, ctx);

    for (i = lenA - 1; i >= n && status == GR_SUCCESS; i--)
    {
        gr_ptr ri = GR_ENTRY(r, i, sz);

        if (gr_is_zero(ri, ctx) == T_TRUE)
        {
            if (Q != NULL)
                status |= gr_zero(GR_ENTRY(Q, i - n, sz), ctx);
            continue;
        }

        /* quotient coefficient r_i / lc; the update uses the coefficients
           c_k = -f_k / lc, so r_i itself multiplies them */
        if (Q != NULL)
        {
            if (P->monic)
                status |= gr_set(GR_ENTRY(Q, i - n, sz), ri, ctx);
            else
                status |= gr_mul(GR_ENTRY(Q, i - n, sz), ri, P->lcinv, ctx);
        }

        status |= gr_set(c, ri, ctx);

        for (k = 0; k < nz; k++)
            status |= gr_addmul(GR_ENTRY(r, P->exps[k] + i - n, sz), c, GR_ENTRY(P->coeffs, k, sz), ctx);
    }

    status |= _gr_vec_set(R, r, n, ctx);

    GR_TMP_CLEAR_VEC(r, lenA, ctx);
    GR_TMP_CLEAR(c, ctx);

    return status;
}

/*
    Generic division by a precomputed modulus: Q of length lenA - lenf + 1
    and R of length lenf - 1 (Q may be NULL for the sparse representation).
    Requires lenA >= lenf.
*/
int
gr_generic_poly_divrem_preinv(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA,
    const gr_poly_preinv_struct * P, gr_ctx_t ctx)
{
    switch (P->kind)
    {
        case GR_POLY_PREINV_SPARSE:
            return _gr_poly_preinv_divrem_sparse_generic(Q, R, A, lenA, P, ctx);

        case GR_POLY_PREINV_TRANSFORMED:
            if (lenA <= 2 * P->lenf - 2)
                return _gr_poly_preinv_divrem_transformed(Q, R, A, lenA, P, ctx);
            /* fall through to Newton division for long dividends */
            /* fallthrough */

        case GR_POLY_PREINV_NEWTON:
        default:
            {
                gr_ptr QQ = Q;
                int status;

                if (Q == NULL)
                    GR_TMP_INIT_VEC(QQ, lenA - P->lenf + 1, ctx);

                /* no inverse: ordinary division */
                if (P->lenfinv == 0)
                {
                    status = _gr_poly_divrem(QQ, R, A, lenA, P->f, P->lenf, ctx);
                }
                else if (lenA <= 2 * P->lenf - 1)
                {
                    status = _gr_poly_divrem_newton_n_preinv(QQ, R, A, lenA, P->f, P->lenf, P->finv, P->lenfinv, ctx);
                }
                else
                {
                    /* The inverse only supports quotients of length at
                       most lenf; reduce long dividends blockwise from
                       the top. */
                    slong lenf = P->lenf, k, m, sz = ctx->sizeof_elem;
                    gr_ptr W, Rt;

                    GR_TMP_INIT_VEC(W, lenA + lenf - 1, ctx);
                    Rt = GR_ENTRY(W, lenA, sz);
                    status = _gr_vec_set(W, A, lenA, ctx);
                    m = lenA;

                    while (m > 2 * lenf - 1 && status == GR_SUCCESS)
                    {
                        k = m - (2 * lenf - 1);
                        status |= _gr_poly_divrem_newton_n_preinv(GR_ENTRY(QQ, k, sz), Rt,
                            GR_ENTRY(W, k, sz), 2 * lenf - 1, P->f, lenf, P->finv, P->lenfinv, ctx);
                        status |= _gr_vec_set(GR_ENTRY(W, k, sz), Rt, lenf - 1, ctx);
                        m = k + lenf - 1;
                    }

                    status |= _gr_poly_divrem_newton_n_preinv(QQ, R, W, m, P->f, lenf, P->finv, P->lenfinv, ctx);
                    GR_TMP_CLEAR_VEC(W, lenA + lenf - 1, ctx);
                }

                if (Q == NULL)
                    GR_TMP_CLEAR_VEC(QQ, lenA - P->lenf + 1, ctx);

                return status;
            }
    }
}

int
_gr_poly_preinv_rem(gr_ptr R, gr_srcptr A, slong lenA, const gr_poly_preinv_t P, gr_ctx_t ctx)
{
    return _gr_poly_preinv_divrem(NULL, R, A, lenA, P, ctx);
}

int
gr_poly_preinv_divrem(gr_poly_t Q, gr_poly_t R, const gr_poly_t A, const gr_poly_preinv_t P, gr_ctx_t ctx)
{
    slong lenA = A->length, lenf = P->lenf;
    int status = GR_SUCCESS;

    if (lenf == 0)
        return GR_DOMAIN;

    if (lenA < lenf)
    {
        status |= gr_poly_set(R, A, ctx);
        status |= gr_poly_zero(Q, ctx);
        return status;
    }

    if (Q == A || R == A || Q == R)
    {
        gr_poly_t t, u;
        gr_poly_init(t, ctx);
        gr_poly_init(u, ctx);
        status |= gr_poly_preinv_divrem(t, u, A, P, ctx);
        gr_poly_swap(Q, t, ctx);
        gr_poly_swap(R, u, ctx);
        gr_poly_clear(t, ctx);
        gr_poly_clear(u, ctx);
        return status;
    }

    gr_poly_fit_length(Q, lenA - lenf + 1, ctx);
    gr_poly_fit_length(R, lenf - 1, ctx);
    status |= _gr_poly_preinv_divrem(Q->coeffs, R->coeffs, A->coeffs, lenA, P, ctx);
    _gr_poly_set_length_normalise(Q, lenA - lenf + 1, ctx);
    _gr_poly_set_length_normalise(R, lenf - 1, ctx);

    return status;
}

int
gr_poly_preinv_rem(gr_poly_t R, const gr_poly_t A, const gr_poly_preinv_t P, gr_ctx_t ctx)
{
    slong lenA = A->length, lenf = P->lenf;
    int status = GR_SUCCESS;

    if (lenf == 0)
        return GR_DOMAIN;

    if (lenA < lenf)
        return gr_poly_set(R, A, ctx);

    if (R == A)
    {
        gr_poly_t t;
        gr_poly_init(t, ctx);
        status |= gr_poly_preinv_rem(t, A, P, ctx);
        gr_poly_swap(R, t, ctx);
        gr_poly_clear(t, ctx);
        return status;
    }

    gr_poly_fit_length(R, lenf - 1, ctx);
    status |= _gr_poly_preinv_rem(R->coeffs, A->coeffs, lenA, P, ctx);
    _gr_poly_set_length_normalise(R, lenf - 1, ctx);

    return status;
}

/* res = poly1 * poly2 mod f. Requires len1, len2 < lenf; writes
   lenf - 1 coefficients. Aliasing of res with the inputs is allowed
   (the product is formed in scratch space). */
int
_gr_poly_preinv_mulmod(gr_ptr res, gr_srcptr poly1, slong len1,
    gr_srcptr poly2, slong len2, const gr_poly_preinv_t P, gr_ctx_t ctx)
{
    slong lenf = P->lenf, lenT;
    gr_ptr T;
    int status = GR_SUCCESS;

    if (len1 == 0 || len2 == 0)
        return _gr_vec_zero(res, lenf - 1, ctx);

    lenT = len1 + len2 - 1;

    if (lenT < lenf)
    {
        /* no reduction necessary */
        if (res == poly1 || res == poly2)
        {
            GR_TMP_INIT_VEC(T, lenT, ctx);
            status |= _gr_poly_mul(T, poly1, len1, poly2, len2, ctx);
            status |= _gr_vec_set(res, T, lenT, ctx);
            GR_TMP_CLEAR_VEC(T, lenT, ctx);
        }
        else
        {
            status |= _gr_poly_mul(res, poly1, len1, poly2, len2, ctx);
        }
        status |= _gr_vec_zero(GR_ENTRY(res, lenT, ctx->sizeof_elem), lenf - 1 - lenT, ctx);
        return status;
    }

    GR_TMP_INIT_VEC(T, lenT, ctx);

    if (P->kind == GR_POLY_PREINV_TRANSFORMED)
    {
        /* form the product in the transformed ring */
        gr_ctx_struct * tctx = P->tctx;
        gr_ptr x, y;
        slong lenp;

        preinv_scratch_entry * e = _preinv_pool_acquire(P, ctx);
        x = e->t[0];
        y = e->t[1];
        status |= _gr_set_gr_poly(x, poly1, len1, ctx, tctx);
        status |= _gr_set_gr_poly(y, poly2, len2, ctx, tctx);
        status |= gr_mul(x, x, y, tctx);
        status |= _gr_get_gr_poly_destructive(T, &lenp, x, ctx, tctx);
        status |= _gr_vec_zero(GR_ENTRY(T, lenp, ctx->sizeof_elem), lenT - lenp, ctx);
        _preinv_pool_release(P, e);
    }
    else
    {
        status |= _gr_poly_mul(T, poly1, len1, poly2, len2, ctx);
    }

    status |= _gr_poly_preinv_rem(res, T, lenT, P, ctx);
    GR_TMP_CLEAR_VEC(T, lenT, ctx);

    return status;
}

int
gr_poly_preinv_mulmod(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2,
    const gr_poly_preinv_t P, gr_ctx_t ctx)
{
    slong lenf = P->lenf;
    int status = GR_SUCCESS;

    if (lenf == 0)
        return GR_DOMAIN;

    if (lenf == 1 || poly1->length == 0 || poly2->length == 0)
        return gr_poly_zero(res, ctx);

    if (poly1->length >= lenf || poly2->length >= lenf)
    {
        gr_poly_t t1, t2;
        gr_poly_init(t1, ctx);
        gr_poly_init(t2, ctx);
        status |= gr_poly_preinv_rem(t1, poly1, P, ctx);
        status |= gr_poly_preinv_rem(t2, poly2, P, ctx);
        if (status == GR_SUCCESS)
            status |= gr_poly_preinv_mulmod(res, t1, t2, P, ctx);
        gr_poly_clear(t1, ctx);
        gr_poly_clear(t2, ctx);
        return status;
    }

    gr_poly_fit_length(res, lenf - 1, ctx);
    status |= _gr_poly_preinv_mulmod(res->coeffs, poly1->coeffs, poly1->length,
        poly2->coeffs, poly2->length, P, ctx);
    _gr_poly_set_length_normalise(res, lenf - 1, ctx);

    return status;
}
