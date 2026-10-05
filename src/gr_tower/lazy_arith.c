/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Lazy fields: arithmetic, comparisons and inverses. */

/* (for recursive mutexes: pthread_mutexattr_settype) */
#define _GNU_SOURCE

#include "gr_tower/lazy_impl.h"

/* Installs a computed element r (in U's current context) into res. */
void
_gr_tower_lazy_install(gr_tower_lazy_elem_t res, gr_tower_flat_struct * U, fmpz_mpoly_q_t r, gr_ctx_t ctx)
{
    if (res->repr == LAZY_REPR_FLAT && res->elem.flat.mctx == U->mctx)
    {
        /* (the old data goes to r, which the caller clears in U's
           context: the same) */
        fmpz_mpoly_q_swap(&res->elem.flat.data, r, U->mctx);
    }
    else
    {
        _gr_tower_lazy_clear_data(res);
        res->repr = LAZY_REPR_FLAT;
        LAZY_POISON_DENSE(res);
        LAZY_POISON_RATIONAL(res);
        res->elem.flat.mctx = U->mctx;
        fmpz_mpoly_q_init(&res->elem.flat.data, res->elem.flat.mctx);
        fmpz_mpoly_q_swap(&res->elem.flat.data, r, res->elem.flat.mctx);
    }
    _gr_tower_lazy_ref(U);
    _gr_tower_lazy_unref(LAZY(ctx), res->F);
    res->F = U;
    res->shallow = 0;
    res->reduced_version = U->ideal_version;
    _gr_tower_lazy_shrink(res);
    _gr_tower_lazy_normalize_rational(res, ctx);
}

/* res = x, moved to a fresh copy of the prefix of its tower in which it
   lives (registered in the context) */
void
_gr_tower_lazy_prefix_copy(gr_tower_lazy_elem_t res, gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    gr_tower_flat_struct * G;
    fmpz_mpoly_q_t u;

    x = _gr_tower_lazy_flat_view(x);
    G = _gr_tower_lazy_new_tower(ctx);
    gr_tower_set_prefix(G->T, x->F->T, x->level);
    gr_tower_flat_ensure(G);
    fmpz_mpoly_q_init(u, G->mctx);
    _gr_tower_flat_transport(u, &x->elem.flat.data, x->F, G);
    _gr_tower_lazy_fresh(res, G, G->T->num_gens, ctx);
    _gr_tower_lazy_install(res, G, u, ctx);
    fmpz_mpoly_q_clear(u, G->mctx);
}

/* Whether a polynomial involves only transcendental variables of F. */
static int
_flat_poly_is_trans(const fmpz_mpoly_t p, gr_tower_flat_struct * F)
{
    int * used;
    slong d;
    int ok = 1;

    if (fmpz_mpoly_is_fmpz(p, F->mctx))
        return 1;

    used = flint_malloc(sizeof(int) * F->cap);
    fmpz_mpoly_used_vars(used, p, F->mctx);
    for (d = 0; d < F->T->num_gens && ok; d++)
        if (used[GR_TOWER_FLAT_VAR_D(F, d)] && F->T->gens[d].kind == GR_TOWER_ALGEBRAIC)
            ok = 0;
    /* variables of no generator (spare capacity) cannot occur */
    flint_free(used);
    return ok;
}

/*
    Whether the sum or difference of x and y (both in U) is reduced
    without further work: both are reduced with respect to the current
    ideal and have denominators in the transcendental variables only
    (multiplying a reduced numerator by such a denominator, and adding,
    keeps the degrees in the algebraic variables below those of the
    triangular set, whose leading coefficients are in the transcendental
    variables too).
*/
static int
_gr_tower_lazy_sum_is_reduced(const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_tower_flat_struct * U)
{
    return x->F == U && y->F == U &&
           x->reduced_version == U->ideal_version && y->reduced_version == U->ideal_version &&
           _flat_poly_is_trans(fmpz_mpoly_q_denref(&x->elem.flat.data), U) &&
           _flat_poly_is_trans(fmpz_mpoly_q_denref(&y->elem.flat.data), U);
}

/*
    The operands of a binary operation in a common tower U, as pointers:
    to the data of the elements themselves when both already live in U
    (no copies), otherwise to the mapped copies xx, yy (initialized; the
    caller clears them when *copied is set).
*/
static int
_gr_tower_lazy_operands(gr_tower_flat_struct ** U, const fmpz_mpoly_q_struct ** px, const fmpz_mpoly_q_struct ** py,
    fmpz_mpoly_q_t xx, fmpz_mpoly_q_t yy, int * copied, gr_tower_lazy_elem_t x, gr_tower_lazy_elem_t y, gr_ctx_t ctx)
{
    int status;

    x = _gr_tower_lazy_flat_view(x);
    y = _gr_tower_lazy_flat_view(y);

    if (x->F == y->F)
    {
        *U = x->F;
        *px = &x->elem.flat.data;
        *py = &y->elem.flat.data;
        *copied = 0;
        return GR_SUCCESS;
    }

    /* (a rational operand: a constant in the tower of the other one) */
    if (x->level == 0 || y->level == 0)
    {
        gr_tower_lazy_elem_struct * r = (x->level == 0) ? x : y;
        gr_tower_lazy_elem_struct * o = (r == x) ? y : x;

        if (fmpz_mpoly_q_is_fmpq(&r->elem.flat.data, r->elem.flat.mctx))
        {
            fmpz_t c;
            fmpz_mpoly_q_struct * rr = (r == x) ? xx : yy;
            fmpz_init(c);
            *U = o->F;
            fmpz_mpoly_q_init(xx, o->elem.flat.mctx);
            fmpz_mpoly_q_init(yy, o->elem.flat.mctx);
            fmpz_mpoly_get_fmpz(c, fmpz_mpoly_q_numref(&r->elem.flat.data), r->elem.flat.mctx);
            fmpz_mpoly_set_fmpz(fmpz_mpoly_q_numref(rr), c, o->elem.flat.mctx);
            fmpz_mpoly_get_fmpz(c, fmpz_mpoly_q_denref(&r->elem.flat.data), r->elem.flat.mctx);
            fmpz_mpoly_set_fmpz(fmpz_mpoly_q_denref(rr), c, o->elem.flat.mctx);
            fmpz_clear(c);
            *px = (r == x) ? xx : &x->elem.flat.data;
            *py = (r == y) ? yy : &y->elem.flat.data;
            *copied = 1;
            return GR_SUCCESS;
        }
    }

    *U = NULL;
    status = _gr_tower_lazy_common(U, xx, yy, x, y, ctx);
    *px = xx;
    *py = yy;
    /* (nothing is allocated when no common tower was found) */
    *copied = (*U != NULL);
    return status;
}

/* res = x + y or x - y for integer denominators (Henrici: the gcd of the
   numerator and the denominator divides gcd(dx, dy)); returns 0 if a
   denominator is not an integer */
/* per-thread scratch polynomials for sums and products (their storage is
   exchanged with the results: no allocation in the common case); an
   fmpz_mpoly's exponent storage is sized for its context, so a slot is
   reset when used in another one */
static FLINT_TLS_PREFIX fmpz_mpoly_struct _gr_tower_lazy_scratch[3];
static FLINT_TLS_PREFIX const fmpz_mpoly_ctx_struct * _gr_tower_lazy_scratch_mctx[3];
static FLINT_TLS_PREFIX slong _gr_tower_lazy_scratch_nvars[3];   /* (against a context reallocated at the same address) */
static FLINT_TLS_PREFIX int _gr_tower_lazy_scratch_ready = 0;

static void
_gr_tower_lazy_scratch_clear(slong i)
{
    _fmpz_vec_clear(_gr_tower_lazy_scratch[i].coeffs, _gr_tower_lazy_scratch[i].alloc);
    flint_free(_gr_tower_lazy_scratch[i].exps);
    memset(_gr_tower_lazy_scratch + i, 0, sizeof(fmpz_mpoly_struct));
    _gr_tower_lazy_scratch[i].bits = MPOLY_MIN_BITS;
    _gr_tower_lazy_scratch_mctx[i] = NULL;
}

static void
_gr_tower_lazy_scratch_cleanup(void)
{
    slong i;
    if (!_gr_tower_lazy_scratch_ready)
        return;
    for (i = 0; i < 3; i++)
        _gr_tower_lazy_scratch_clear(i);
    _gr_tower_lazy_scratch_ready = 0;
}

static fmpz_mpoly_struct *
_gr_tower_lazy_scratch_get(slong i, const fmpz_mpoly_ctx_t mctx)
{
    if (!_gr_tower_lazy_scratch_ready)
    {
        slong j;
        memset(_gr_tower_lazy_scratch, 0, sizeof(_gr_tower_lazy_scratch));
        for (j = 0; j < 3; j++)
        {
            _gr_tower_lazy_scratch[j].bits = MPOLY_MIN_BITS;
            _gr_tower_lazy_scratch_mctx[j] = NULL;
        }
        _gr_tower_lazy_scratch_ready = 1;
        flint_register_cleanup_function(_gr_tower_lazy_scratch_cleanup);
    }
    if (_gr_tower_lazy_scratch_mctx[i] != mctx || _gr_tower_lazy_scratch_nvars[i] != mctx->minfo->nvars)
    {
        _gr_tower_lazy_scratch_clear(i);
        _gr_tower_lazy_scratch_mctx[i] = mctx;
        _gr_tower_lazy_scratch_nvars[i] = mctx->minfo->nvars;
    }
    return _gr_tower_lazy_scratch + i;
}

static int
_flat_add_int_den(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, int subtract, const fmpz_mpoly_ctx_t mctx)
{
    const fmpz_mpoly_struct * dxp = fmpz_mpoly_q_denref(x), * dyp = fmpz_mpoly_q_denref(y);
    fmpz_t dx, dy, g, a, b;
    fmpz_mpoly_struct * num = fmpz_mpoly_q_numref(res);
    slong i;

    if (!fmpz_mpoly_is_fmpz(dxp, mctx) || !fmpz_mpoly_is_fmpz(dyp, mctx))
        return 0;

    fmpz_init(dx);
    fmpz_init(dy);
    fmpz_init(g);
    fmpz_mpoly_get_fmpz(dx, dxp, mctx);
    fmpz_mpoly_get_fmpz(dy, dyp, mctx);

    if (fmpz_equal(dx, dy))
    {
        /* (into the scratch polynomial, whose storage the result then
           exchanges: no temporary for in-place sums) */
        fmpz_mpoly_struct * S = _gr_tower_lazy_scratch_get(0, mctx);
        if (subtract)
            fmpz_mpoly_sub(S, fmpz_mpoly_q_numref(x), fmpz_mpoly_q_numref(y), mctx);
        else
            fmpz_mpoly_add(S, fmpz_mpoly_q_numref(x), fmpz_mpoly_q_numref(y), mctx);
        fmpz_mpoly_swap(num, S, mctx);
        fmpz_set(g, dx);
    }
    else
    {
        fmpz_mpoly_struct * S = _gr_tower_lazy_scratch_get(0, mctx), * t = _gr_tower_lazy_scratch_get(1, mctx);
        fmpz_init(a);
        fmpz_init(b);
        fmpz_gcd(g, dx, dy);
        fmpz_divexact(a, dy, g);
        fmpz_divexact(b, dx, g);
        fmpz_mpoly_scalar_mul_fmpz(t, fmpz_mpoly_q_numref(y), b, mctx);
        if (fmpz_is_one(a))
        {
            if (subtract)
                fmpz_mpoly_sub(S, fmpz_mpoly_q_numref(x), t, mctx);
            else
                fmpz_mpoly_add(S, fmpz_mpoly_q_numref(x), t, mctx);
        }
        else
        {
            fmpz_mpoly_struct * R = _gr_tower_lazy_scratch_get(2, mctx);
            fmpz_mpoly_scalar_mul_fmpz(R, fmpz_mpoly_q_numref(x), a, mctx);
            if (subtract)
                fmpz_mpoly_sub(S, R, t, mctx);
            else
                fmpz_mpoly_add(S, R, t, mctx);
        }
        fmpz_mpoly_swap(num, S, mctx);
        fmpz_mul(dx, dx, a);
        fmpz_clear(a);
        fmpz_clear(b);
    }

    if (fmpz_mpoly_is_zero(num, mctx))
        fmpz_one(dx);
    else if (!fmpz_is_one(g))
    {
        /* gcd of g and the coefficients, stopping at 1 */
        fmpz_abs(g, g);
        for (i = 0; i < num->length && !fmpz_is_one(g); i++)
            fmpz_gcd(g, g, num->coeffs + i);
        if (!fmpz_is_one(g))
        {
            fmpz_mpoly_scalar_divexact_fmpz(num, num, g, mctx);
            fmpz_divexact(dx, dx, g);
        }
    }

    fmpz_mpoly_set_fmpz(fmpz_mpoly_q_denref(res), dx, mctx);
    fmpz_clear(dx);
    fmpz_clear(dy);
    fmpz_clear(g);
    return 1;
}

/* (whether the result of an operation in U can be written directly to
   res: the operations below allow aliasing) */
#define LAZY_IN_PLACE(res, U) ((res)->F == (U) && (res)->repr == LAZY_REPR_FLAT && (res)->elem.flat.mctx == (U)->mctx && !(res)->shallow)

/* after an in-place operation on res in its own tower */
static void
_gr_tower_lazy_install_in_place(gr_tower_lazy_elem_t res, int status)
{
    if (status == GR_SUCCESS)
    {
        res->reduced_version = res->F->ideal_version;
        _gr_tower_lazy_shrink(res);
    }
    else
        res->reduced_version = 0;
}

#define LAZY_ADDITIVE_OP(name, opfunc, subtract) \
int \
name(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx) \
{ \
    gr_tower_flat_struct * U; \
    fmpz_mpoly_q_t xx, yy, rt; \
    fmpz_mpoly_q_struct * r; \
    const fmpz_mpoly_q_struct * px, * py; \
    int status, reduced, copied, inplace; \
    x = _gr_tower_lazy_flat_view(x); \
    y = _gr_tower_lazy_flat_view(y); \
    status = _gr_tower_lazy_operands(&U, &px, &py, xx, yy, &copied, (gr_tower_lazy_elem_struct *) x, (gr_tower_lazy_elem_struct *) y, ctx); \
    if (status != GR_SUCCESS) \
    { \
        if (copied) { fmpz_mpoly_q_clear(xx, U->mctx); fmpz_mpoly_q_clear(yy, U->mctx); } \
        return status; \
    } \
    reduced = _gr_tower_lazy_sum_is_reduced(x, y, U); \
    inplace = LAZY_IN_PLACE(res, U); \
    if (inplace) \
        r = &res->elem.flat.data; \
    else \
    { \
        fmpz_mpoly_q_init(rt, U->mctx); \
        r = rt; \
    } \
    if (!_flat_add_int_den(r, px, py, subtract, U->mctx)) \
        opfunc(r, px, py, U->mctx); \
    if (!reduced) \
        status = gr_tower_flat_reduce(r, U); \
    if (inplace) \
    { _gr_tower_lazy_install_in_place(res, status); _gr_tower_lazy_normalize_rational(res, ctx); } \
    else \
    { \
        if (status == GR_SUCCESS) \
            _gr_tower_lazy_install(res, U, r, ctx); \
        fmpz_mpoly_q_clear(rt, U->mctx); \
    } \
    if (copied) { fmpz_mpoly_q_clear(xx, U->mctx); fmpz_mpoly_q_clear(yy, U->mctx); } \
    return status; \
}

LAZY_ADDITIVE_OP(_gr_tower_lazy_add, fmpz_mpoly_q_add, 0)
LAZY_ADDITIVE_OP(_gr_tower_lazy_sub, fmpz_mpoly_q_sub, 1)
static int _flat_is_rational(const fmpz_mpoly_q_t x, gr_tower_flat_struct * F);

/* res = x * y for integer denominators (canonical by Gauss's lemma:
   cancelling gcd(content(nx), dy) and gcd(content(ny), dx) suffices, no
   polynomial gcd); returns 0 if a denominator is not an integer */
static int
_flat_mul_int_den(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, const fmpz_mpoly_ctx_t mctx)
{
    const fmpz_mpoly_struct * nx = fmpz_mpoly_q_numref(x), * ny = fmpz_mpoly_q_numref(y);
    const fmpz_mpoly_struct * dxp = fmpz_mpoly_q_denref(x), * dyp = fmpz_mpoly_q_denref(y);
    fmpz_mpoly_struct * S;
    fmpz_t dx, dy, g1, g2;

    if (!fmpz_mpoly_is_fmpz(dxp, mctx) || !fmpz_mpoly_is_fmpz(dyp, mctx))
        return 0;

    if (nx->length == 0 || ny->length == 0)
    {
        fmpz_mpoly_q_zero(res, mctx);
        return 1;
    }

    fmpz_init(dx);
    fmpz_init(dy);
    fmpz_init(g1);
    fmpz_init(g2);
    fmpz_mpoly_get_fmpz(dx, dxp, mctx);
    fmpz_mpoly_get_fmpz(dy, dyp, mctx);

    S = _gr_tower_lazy_scratch_get(0, mctx);
    fmpz_mpoly_mul(S, nx, ny, mctx);

    if (!fmpz_is_one(dy))
    {
        _fmpz_vec_content(g1, nx->coeffs, nx->length);
        fmpz_gcd(g1, g1, dy);
    }
    else
        fmpz_one(g1);
    if (!fmpz_is_one(dx))
    {
        _fmpz_vec_content(g2, ny->coeffs, ny->length);
        fmpz_gcd(g2, g2, dx);
    }
    else
        fmpz_one(g2);

    fmpz_mul(g1, g1, g2);           /* (the common factor) */
    fmpz_mul(dx, dx, dy);
    if (!fmpz_is_one(g1))
    {
        _fmpz_vec_scalar_divexact_fmpz(S->coeffs, S->coeffs, S->length, g1);
        fmpz_divexact(dx, dx, g1);
    }

    fmpz_mpoly_swap(fmpz_mpoly_q_numref(res), S, mctx);
    fmpz_mpoly_set_fmpz(fmpz_mpoly_q_denref(res), dx, mctx);

    fmpz_clear(dx);
    fmpz_clear(dy);
    fmpz_clear(g1);
    fmpz_clear(g2);
    return 1;
}

int
_gr_tower_lazy_mul(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx)
{
    gr_tower_flat_struct * U;
    fmpz_mpoly_q_t xx, yy, rt;
    fmpz_mpoly_q_struct * r;
    const fmpz_mpoly_q_struct * px, * py;
    int status, copied, inplace;

    x = _gr_tower_lazy_flat_view(x);
    y = _gr_tower_lazy_flat_view(y);
    status = _gr_tower_lazy_operands(&U, &px, &py, xx, yy, &copied, (gr_tower_lazy_elem_struct *) x, (gr_tower_lazy_elem_struct *) y, ctx);
    if (status != GR_SUCCESS)
    {
        if (copied) { fmpz_mpoly_q_clear(xx, U->mctx); fmpz_mpoly_q_clear(yy, U->mctx); }
        return status;
    }

    inplace = LAZY_IN_PLACE(res, U);
    if (inplace)
        r = &res->elem.flat.data;
    else
    {
        fmpz_mpoly_q_init(rt, U->mctx);
        r = rt;
    }

    /* a rational factor: a scalar multiplication (no reduction) */
    if (_flat_is_rational(py, U) || _flat_is_rational(px, U))
    {
        if (!_flat_mul_int_den(r, px, py, U->mctx))
            fmpz_mpoly_q_mul(r, px, py, U->mctx);
    }
    /* number fields with univariate moduli: dense products (reduced) */
    else if (!(_gr_tower_flat_dense_applicable(U) && !_gr_tower_flat_mul_prefers_sparse(px, py, U) &&
               _gr_tower_flat_mul_dense(r, px, py, U)))
    {
        if (!_flat_mul_int_den(r, px, py, U->mctx))
            fmpz_mpoly_q_mul(r, px, py, U->mctx);
        status = gr_tower_flat_reduce(r, U);
    }

    if (inplace)
    {
        _gr_tower_lazy_install_in_place(res, status);
        _gr_tower_lazy_normalize_rational(res, ctx);
    }
    else
    {
        if (status == GR_SUCCESS)
            _gr_tower_lazy_install(res, U, r, ctx);
        fmpz_mpoly_q_clear(rt, U->mctx);
    }
    if (copied) { fmpz_mpoly_q_clear(xx, U->mctx); fmpz_mpoly_q_clear(yy, U->mctx); }
    return status;
}

/* arithmetic with a rational scalar, in the tower of x (the result
   stays reduced; no merging of towers, and no temporary element) */
#define LAZY_SCALAR_ADD 0
#define LAZY_SCALAR_SUB 1
#define LAZY_SCALAR_MUL 2
#define LAZY_SCALAR_DIV 3

/* Prepares res for an operation modifying its flat data in place: the
   data is brought to the current context (a shallow copy is replaced by
   a private one) and the dense form, which would no longer agree with
   it, is dropped. */
static void
_gr_tower_lazy_prepare_in_place(gr_tower_lazy_elem_t res, gr_ctx_t ctx)
{
    if (res->shallow)
    {
        gr_tower_lazy_elem_struct t;
        _gr_tower_lazy_init(&t, ctx);
        GR_MUST_SUCCEED(_gr_tower_lazy_set(&t, res, ctx));
        *res = t;
    }
    _gr_tower_lazy_update(res);
}

static int
_gr_tower_lazy_scalar_fmpz(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const fmpz_t c, int op, gr_ctx_t ctx)
{
    if (op >= LAZY_SCALAR_MUL && fmpz_is_zero(c))
        return (op == LAZY_SCALAR_DIV) ? GR_DOMAIN : _gr_tower_lazy_zero(res, ctx);
    if (op <= LAZY_SCALAR_SUB && fmpz_is_zero(c))
        return _gr_tower_lazy_set(res, x, ctx);
    if (op >= LAZY_SCALAR_MUL && fmpz_is_one(c))
        return _gr_tower_lazy_set(res, x, ctx);

    if (res != x)
        _gr_tower_lazy_set(res, x, ctx);
    _gr_tower_lazy_prepare_in_place(res, ctx);
    if (op <= LAZY_SCALAR_SUB && fmpz_mpoly_is_fmpz(fmpz_mpoly_q_denref(&res->elem.flat.data), res->elem.flat.mctx))
    {
        /* an integer denominator d: num + c d (still coprime to d) */
        fmpz_t t;
        fmpz_init(t);
        fmpz_mpoly_get_fmpz(t, fmpz_mpoly_q_denref(&res->elem.flat.data), res->elem.flat.mctx);
        fmpz_mul(t, t, c);
        if (op == LAZY_SCALAR_SUB)
            fmpz_neg(t, t);
        fmpz_mpoly_add_fmpz(fmpz_mpoly_q_numref(&res->elem.flat.data), fmpz_mpoly_q_numref(&res->elem.flat.data), t, res->elem.flat.mctx);
        if (fmpz_mpoly_is_zero(fmpz_mpoly_q_numref(&res->elem.flat.data), res->elem.flat.mctx))
            fmpz_mpoly_one(fmpz_mpoly_q_denref(&res->elem.flat.data), res->elem.flat.mctx);
        fmpz_clear(t);
    }
    else if (op == LAZY_SCALAR_ADD)
        fmpz_mpoly_q_add_fmpz(&res->elem.flat.data, &res->elem.flat.data, c, res->elem.flat.mctx);
    else if (op == LAZY_SCALAR_SUB)
        fmpz_mpoly_q_sub_fmpz(&res->elem.flat.data, &res->elem.flat.data, c, res->elem.flat.mctx);
    else if (op == LAZY_SCALAR_MUL)
        fmpz_mpoly_q_mul_fmpz(&res->elem.flat.data, &res->elem.flat.data, c, res->elem.flat.mctx);
    else
        fmpz_mpoly_q_div_fmpz(&res->elem.flat.data, &res->elem.flat.data, c, res->elem.flat.mctx);
    _gr_tower_lazy_normalize_rational(res, ctx);
    return GR_SUCCESS;
}

static int
_gr_tower_lazy_scalar_fmpq(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const fmpq_t c, int op, gr_ctx_t ctx)
{
    if (fmpz_is_one(fmpq_denref(c)))
        return _gr_tower_lazy_scalar_fmpz(res, x, fmpq_numref(c), op, ctx);
    if (res != x)
        _gr_tower_lazy_set(res, x, ctx);
    _gr_tower_lazy_prepare_in_place(res, ctx);
    if (op == LAZY_SCALAR_ADD)
        fmpz_mpoly_q_add_fmpq(&res->elem.flat.data, &res->elem.flat.data, c, res->elem.flat.mctx);
    else if (op == LAZY_SCALAR_SUB)
        fmpz_mpoly_q_sub_fmpq(&res->elem.flat.data, &res->elem.flat.data, c, res->elem.flat.mctx);
    else if (op == LAZY_SCALAR_MUL)
        fmpz_mpoly_q_mul_fmpq(&res->elem.flat.data, &res->elem.flat.data, c, res->elem.flat.mctx);
    else
        fmpz_mpoly_q_div_fmpq(&res->elem.flat.data, &res->elem.flat.data, c, res->elem.flat.mctx);
    _gr_tower_lazy_normalize_rational(res, ctx);
    return GR_SUCCESS;
}

#define LAZY_SCALAR_METHODS(opname, op) \
int _gr_tower_lazy_##opname##_fmpz(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const fmpz_t c, gr_ctx_t ctx) \
{ return _gr_tower_lazy_scalar_fmpz(res, x, c, op, ctx); } \
int _gr_tower_lazy_##opname##_fmpq(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const fmpq_t c, gr_ctx_t ctx) \
{ return _gr_tower_lazy_scalar_fmpq(res, x, c, op, ctx); } \
int _gr_tower_lazy_##opname##_si(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, slong c, gr_ctx_t ctx) \
{ fmpz_t t; int status; fmpz_init_set_si(t, c); status = _gr_tower_lazy_scalar_fmpz(res, x, t, op, ctx); fmpz_clear(t); return status; } \
int _gr_tower_lazy_##opname##_ui(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, ulong c, gr_ctx_t ctx) \
{ fmpz_t t; int status; fmpz_init_set_ui(t, c); status = _gr_tower_lazy_scalar_fmpz(res, x, t, op, ctx); fmpz_clear(t); return status; }

LAZY_SCALAR_METHODS(add, LAZY_SCALAR_ADD)
LAZY_SCALAR_METHODS(sub, LAZY_SCALAR_SUB)
LAZY_SCALAR_METHODS(mul, LAZY_SCALAR_MUL)
LAZY_SCALAR_METHODS(div, LAZY_SCALAR_DIV)

int
_gr_tower_lazy_neg(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    if (res != x)
    {
        int status = _gr_tower_lazy_set(res, x, ctx);
        if (status != GR_SUCCESS) return status;
    }
    if (res->repr == LAZY_REPR_RATIONAL)
        fmpq_neg(&res->elem.q, &res->elem.q);
    else if (res->repr == LAZY_REPR_DENSE)
        fmpq_poly_neg(&res->elem.dense.poly, &res->elem.dense.poly);
    else
    {
        _gr_tower_lazy_prepare_in_place(res, ctx);
        fmpz_mpoly_q_neg(&res->elem.flat.data, &res->elem.flat.data, res->elem.flat.mctx);
    }
    return GR_SUCCESS;
}

/*
    Equality of algebraic numbers in towers whose merge would be large
    (cos(pi/257) as 2956 nested square roots, in a tower of nominal
    degree 4096, against cos(pi/257) in Q(zeta_257): their merge has
    nominal degree 2^20, in which the zero test is hopeless): with P the
    minimal polynomial of the element s of the smaller tower (found by
    linear algebra in that tower), the other element o equals s if and
    only if P(o) = 0 (a zero test in the tower of o alone) and o is the
    root of P which s is (the roots of P being isolated numerically).
    Returns T_UNKNOWN when this does not apply.
*/
static FLINT_TLS_PREFIX slong _gr_tower_lazy_split_depth = 0;

static slong
_gr_tower_lazy_tower_degree_capped(const gr_tower_struct * T, slong cap)
{
    slong k, d = 1;
    for (k = 1; k <= T->length; k++)
    {
        d *= gr_tower_step_degree(T, k);
        if (d > cap)
            return cap + 1;
    }
    return d;
}

static truth_t
_gr_tower_lazy_equal_minpoly(gr_tower_lazy_elem_t a, gr_tower_lazy_elem_t b, gr_ctx_t ctx)
{
    gr_tower_lazy_elem_struct * s, * o;
    gr_tower_struct * T;
    gr_ctx_struct * top;
    gr_ptr xt;
    fmpz_poly_t P;
    slong da, db, i, prec, deg, lim;
    truth_t res = T_UNKNOWN;
    int status;

    a = _gr_tower_lazy_flat_view(a);
    b = _gr_tower_lazy_flat_view(b);
    lim = LAZY(ctx)->options[GR_TOWER_OPT_MINPOLY_DEGREE_LIMIT];
    da = _gr_tower_lazy_tower_degree_capped(a->F->T, lim);
    db = _gr_tower_lazy_tower_degree_capped(b->F->T, lim);
    s = (da <= db) ? a : b;
    o = (da <= db) ? b : a;
    if (FLINT_MIN(da, db) > lim || s->F->T->num_trans != 0)
        return T_UNKNOWN;

    /* the minimal polynomial of s */
    T = s->F->T;
    top = gr_tower_field(T);
    fmpz_poly_init(P);
    GR_TMP_INIT(xt, top);
    {
        slong k = _gr_tower_lazy_alg_level(s);
        gr_ctx_struct * Fk = gr_tower_field_at(T, k);
        gr_ptr u;
        GR_TMP_INIT(u, Fk);
        status = gr_tower_flat_get_nested_at(u, &s->elem.flat.data, k, s->F);
        if (status == GR_SUCCESS)
            status = gr_tower_promote(xt, u, k, T->length, T);
        GR_TMP_CLEAR(u, Fk);
    }
    if (status == GR_SUCCESS)
        status = gr_tower_get_fmpz_poly_minpoly(P, xt, T);
    GR_TMP_CLEAR(xt, top);

    deg = fmpz_poly_degree(P);
    if (status == GR_SUCCESS && deg >= 1)
    {
        gr_tower_lazy_elem_struct p;
        truth_t t;

        /* P(o) */
        _gr_tower_lazy_init(&p, ctx);
        status = _gr_tower_lazy_zero(&p, ctx);
        for (i = deg; i >= 0 && status == GR_SUCCESS; i--)
        {
            status = _gr_tower_lazy_mul(&p, &p, o, ctx);
            if (status == GR_SUCCESS)
                status = _gr_tower_lazy_add_fmpz(&p, &p, P->coeffs + i, ctx);
        }
        t = (status == GR_SUCCESS) ? _gr_tower_lazy_is_zero(&p, ctx) : T_UNKNOWN;
        _gr_tower_lazy_clear(&p, ctx);

        if (t == T_FALSE)
            res = T_FALSE;
        else if (t == T_TRUE)
        {
            /* o is a root of P, as s is: the same one? */
            acb_ptr roots = _acb_vec_init(deg);
            acb_t zs, zo;
            acb_init(zs);
            acb_init(zo);
            for (prec = 64; prec <= 16 * GR_TOWER_DEFAULT_PREC * 64; prec *= 2)
            {
                slong is = -1, io = -1, ns = 0, no = 0;

                arb_fmpz_poly_complex_roots(roots, P, 0, prec);
                if (_gr_tower_lazy_get_acb_impl(zs, s, prec, ctx) != GR_SUCCESS ||
                    _gr_tower_lazy_get_acb_impl(zo, o, prec, ctx) != GR_SUCCESS)
                    break;
                for (i = 0; i < deg; i++)
                {
                    if (acb_overlaps(roots + i, zs))
                        is = i, ns++;
                    if (acb_overlaps(roots + i, zo))
                        io = i, no++;
                }
                if (ns == 1 && no == 1)
                {
                    res = (is == io) ? T_TRUE : T_FALSE;
                    break;
                }
            }
            acb_clear(zs);
            acb_clear(zo);
            _acb_vec_clear(roots, deg);
        }
    }

    fmpz_poly_clear(P);
    return res;
}

/*
    The zero test of an algebraic x = f + g in a large tower (nominal
    degree above GR_TOWER_OPT_SPLIT_DEGREE_LIMIT, not all proven), f and g involving
    disjoint sets of generators (with their definitions): a difference of
    two numbers built independently, whose towers have been stacked by
    the merge. Decided as f = -g in towers of their own generators
    (above). Returns T_UNKNOWN when this does not apply.
*/
static truth_t
_gr_tower_lazy_flat_split_is_zero(const fmpz_mpoly_q_t x, gr_tower_flat_struct * F, gr_ctx_t ctx)
{
    gr_tower_struct * T = F->T;
    const fmpz_mpoly_struct * num = fmpz_mpoly_q_numref(x);
    slong G = T->num_gens, n = num->length, nvars, i, d, e, rA = -1, roots;
    slong lim = LAZY(ctx)->options[GR_TOWER_OPT_SPLIT_DEGREE_LIMIT];
    slong * parent, * term_rep, * mapA, * mapB;
    int * clos, * used, * markA, * markB;
    ulong * exp;
    truth_t res = T_UNKNOWN;

    if (lim <= 0 || _gr_tower_lazy_split_depth > 2 || T->num_trans != 0 || T->length < 2 || n < 2 ||
        !GR_TOWER_HAS_CAP(T, GR_TOWER_CAP_RATIONAL) ||
        !fmpz_mpoly_is_fmpz(fmpz_mpoly_q_denref(x), F->mctx) ||
        _gr_tower_lazy_tower_degree_capped(T, lim) <= lim ||
        _gr_tower_all_proven(T))
        return T_UNKNOWN;

    nvars = F->mctx->minfo->nvars;
    clos = flint_calloc(G * G, sizeof(int));
    used = flint_calloc(G, sizeof(int));
    parent = flint_malloc(sizeof(slong) * G);
    term_rep = flint_malloc(sizeof(slong) * n);
    exp = flint_malloc(sizeof(ulong) * nvars);

    for (d = 0; d < G; d++)
    {
        clos[d * G + d] = 1;
        _gr_tower_involved_gens_closure(clos + d * G, T);
        parent[d] = d;
    }

#define FIND(r) do { while (parent[r] != r) r = parent[r] = parent[parent[r]]; } while (0)

    /* components: the generators of a term (with their definitions)
       are joined */
    for (i = 0; i < n; i++)
    {
        slong rep = -1;
        fmpz_mpoly_get_term_exp_ui(exp, num, i, F->mctx);
        for (d = 0; d < G; d++)
        {
            if (exp[GR_TOWER_FLAT_VAR_D(F, d)] == 0)
                continue;
            for (e = 0; e < G; e++)
            {
                if (!clos[d * G + e])
                    continue;
                used[e] = 1;
                if (rep < 0)
                    rep = e;
                else
                {
                    slong r1 = rep, r2 = e;
                    FIND(r1);
                    FIND(r2);
                    if (r1 != r2)
                        parent[r2] = r1;
                }
            }
        }
        term_rep[i] = rep;
    }

    roots = 0;
    for (d = 0; d < G; d++)
    {
        slong r = d;
        FIND(r);
        if (used[d] && r == d)
            roots++;
    }

    if (roots == 2)
    {
        fmpz_mpoly_q_t fA, fB, u;
        gr_tower_flat_struct * GA, * GB;
        gr_tower_lazy_elem_struct a, b;
        fmpz_t c;
        slong na = 0, nb = 0;

        for (i = 0; i < n && rA < 0; i++)
            if (term_rep[i] >= 0)
            {
                rA = term_rep[i];
                FIND(rA);
            }

        markA = flint_calloc(G, sizeof(int));
        markB = flint_calloc(G, sizeof(int));
        mapA = flint_malloc(sizeof(slong) * G);
        mapB = flint_malloc(sizeof(slong) * G);
        for (d = 0; d < G; d++)
        {
            slong r = d;
            FIND(r);
            markA[d] = used[d] && r == rA;
            markB[d] = used[d] && r != rA;
            mapA[d] = markA[d] ? na++ : -1;
            mapB[d] = markB[d] ? nb++ : -1;
        }

        /* f = the terms of A (with the constant), g = minus those of B */
        fmpz_mpoly_q_init(fA, F->mctx);
        fmpz_mpoly_q_init(fB, F->mctx);
        fmpz_init(c);
        for (i = 0; i < n; i++)
        {
            slong r = term_rep[i];
            if (r >= 0)
                FIND(r);
            fmpz_mpoly_get_term_exp_ui(exp, num, i, F->mctx);
            fmpz_mpoly_get_term_coeff_fmpz(c, num, i, F->mctx);
            if (r < 0 || r == rA)
                fmpz_mpoly_push_term_fmpz_ui(fmpz_mpoly_q_numref(fA), c, exp, F->mctx);
            else
            {
                fmpz_neg(c, c);
                fmpz_mpoly_push_term_fmpz_ui(fmpz_mpoly_q_numref(fB), c, exp, F->mctx);
            }
        }
        fmpz_mpoly_sort_terms(fmpz_mpoly_q_numref(fA), F->mctx);
        fmpz_mpoly_combine_like_terms(fmpz_mpoly_q_numref(fA), F->mctx);
        fmpz_mpoly_sort_terms(fmpz_mpoly_q_numref(fB), F->mctx);
        fmpz_mpoly_combine_like_terms(fmpz_mpoly_q_numref(fB), F->mctx);
        fmpz_mpoly_set(fmpz_mpoly_q_denref(fA), fmpz_mpoly_q_denref(x), F->mctx);
        fmpz_mpoly_set(fmpz_mpoly_q_denref(fB), fmpz_mpoly_q_denref(x), F->mctx);
        fmpz_mpoly_q_canonicalise(fA, F->mctx);
        fmpz_mpoly_q_canonicalise(fB, F->mctx);
        fmpz_clear(c);

        /* the two numbers in towers of their own generators */
        GA = _gr_tower_lazy_new_tower(ctx);
        gr_tower_set_subset(GA->T, T, markA);
        gr_tower_flat_ensure(GA);
        GB = _gr_tower_lazy_new_tower(ctx);
        gr_tower_set_subset(GB->T, T, markB);
        gr_tower_flat_ensure(GB);

        _gr_tower_lazy_init(&a, ctx);
        _gr_tower_lazy_init(&b, ctx);
        fmpz_mpoly_q_init(u, GA->mctx);
        _gr_tower_flat_transport_map(u, fA, F, GA, mapA);
        _gr_tower_lazy_install(&a, GA, u, ctx);
        fmpz_mpoly_q_clear(u, GA->mctx);
        fmpz_mpoly_q_init(u, GB->mctx);
        _gr_tower_flat_transport_map(u, fB, F, GB, mapB);
        _gr_tower_lazy_install(&b, GB, u, ctx);
        fmpz_mpoly_q_clear(u, GB->mctx);
        fmpz_mpoly_q_clear(fA, F->mctx);
        fmpz_mpoly_q_clear(fB, F->mctx);

        _gr_tower_lazy_split_depth++;
        res = _gr_tower_lazy_equal_minpoly(&a, &b, ctx);
        _gr_tower_lazy_split_depth--;

        _gr_tower_lazy_clear(&a, ctx);
        _gr_tower_lazy_clear(&b, ctx);
        flint_free(markA);
        flint_free(markB);
        flint_free(mapA);
        flint_free(mapB);
    }

#undef FIND

    flint_free(clos);
    flint_free(used);
    flint_free(parent);
    flint_free(term_rep);
    flint_free(exp);
    return res;
}

static truth_t
_gr_tower_lazy_flat_num_is_zero(fmpz_mpoly_q_t x, gr_tower_flat_struct * F, gr_ctx_t ctx)
{
    truth_t res = _gr_tower_lazy_flat_split_is_zero(x, F, ctx);
    if (res != T_UNKNOWN)
        return res;
    res = gr_tower_flat_num_is_zero(x, F);
    if (res == T_UNKNOWN)
    {
        /* no structure theorem decides it: only more precision than the
           separation stage used can still show x != 0 (1 - erf(100), the
           tower's erfc(100) = 6.4e-4346, needs 14500 bits) */
        acb_t v;
        slong prec, lim = GR_TOWER_OPTION(F->T, GR_TOWER_OPT_NUMERIC_PREC_LIMIT);
        acb_init(v);
        for (prec = 2 * GR_TOWER_OPTION(F->T, GR_TOWER_OPT_PREC_LIMIT); prec <= lim; prec *= 2)
        {
            if (gr_tower_flat_get_acb(v, x, prec, F) != GR_SUCCESS)
                break;
            if (!acb_contains_zero(v))
            {
                res = T_FALSE;
                break;
            }
        }
        acb_clear(v);
    }
    return res;
}

truth_t
_gr_tower_lazy_is_zero(const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    fmpz_mpoly_q_t t;
    truth_t res;

    x = _gr_tower_lazy_flat_view(x);
    fmpz_mpoly_q_init(t, x->elem.flat.mctx);
    fmpz_mpoly_q_set(t, &x->elem.flat.data, x->elem.flat.mctx);
    res = _gr_tower_lazy_flat_num_is_zero(t, x->F, ctx);
    fmpz_mpoly_q_clear(t, x->F->mctx);   /* num_is_zero may have grown the context: t is in the current one */
    return res;
}

/* whether the enclosures of x and y (each in its own tower) are disjoint
   at a moderate precision: then x != y, without merging their towers
   (comparing two roots of a degree-20 polynomial would otherwise build
   their compositum) */
static int
_gr_tower_lazy_enclosures_disjoint(const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx)
{
    acb_t a, b;
    slong prec;
    int res = 0;

    if (x->F == y->F)
        return 0;
    acb_init(a);
    acb_init(b);
    for (prec = GR_TOWER_DEFAULT_PREC; prec <= LAZY(ctx)->options[GR_TOWER_OPT_PREC_LIMIT] && !res; prec *= 2)
    {
        if (_gr_tower_lazy_get_acb_impl(a, x, prec, ctx) != GR_SUCCESS || _gr_tower_lazy_get_acb_impl(b, y, prec, ctx) != GR_SUCCESS)
            break;
        res = !acb_overlaps(a, b);
    }
    acb_clear(a);
    acb_clear(b);
    return res;
}

truth_t
_gr_tower_lazy_equal(const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx)
{
    gr_tower_flat_struct * U;
    fmpz_mpoly_q_t xx, yy;
    truth_t res;

    if (_gr_tower_lazy_enclosures_disjoint(x, y, ctx))
        return T_FALSE;

    if (_gr_tower_lazy_common(&U, xx, yy, (gr_tower_lazy_elem_struct *) x, (gr_tower_lazy_elem_struct *) y, ctx) != GR_SUCCESS)
        return T_UNKNOWN;

    if (fmpz_mpoly_q_equal(xx, yy, U->mctx))
        res = T_TRUE;
    else
    {
        fmpz_mpoly_q_sub(xx, xx, yy, U->mctx);
        res = _gr_tower_lazy_flat_num_is_zero(xx, U, ctx);
    }

    fmpz_mpoly_q_clear(xx, U->mctx);
    fmpz_mpoly_q_clear(yy, U->mctx);
    return res;
}

truth_t
_gr_tower_lazy_is_one(const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    fmpz_mpoly_q_t t;
    fmpz_t one;
    truth_t res;

    x = _gr_tower_lazy_flat_view(x);
    if (fmpz_mpoly_q_is_one(&x->elem.flat.data, x->elem.flat.mctx))
        return T_TRUE;

    fmpz_init_set_ui(one, 1);
    fmpz_mpoly_q_init(t, x->elem.flat.mctx);
    fmpz_mpoly_q_sub_fmpz(t, &x->elem.flat.data, one, x->elem.flat.mctx);
    res = gr_tower_flat_num_is_zero(t, x->F);
    fmpz_mpoly_q_clear(t, x->F->mctx);
    fmpz_clear(one);
    return res;
}

/*
    Rationalizes the denominator of x (in the current context of F): if
    it involves algebraic generators, it is inverted in the nested
    representation (whose elements have denominators in the
    transcendental generators only) and multiplied into the numerator.
    Sums of many fractions would otherwise accumulate products of
    algebraic denominators which polynomial gcds cannot cancel.
*/
static int
_gr_tower_lazy_rationalize(fmpz_mpoly_q_t x, gr_tower_flat_struct * F)
{
    gr_tower_struct * T = F->T;
    int * used;
    slong k, level = 0;
    int status = GR_SUCCESS;

    if (T->length == 0 || fmpz_mpoly_is_fmpz(fmpz_mpoly_q_denref(x), F->mctx))
        return GR_SUCCESS;

    /* a monomial denominator is cleared through the moduli, a small
       polynomial one through norms */
    if (gr_tower_flat_rationalize(x, F) == GR_SUCCESS)
        return GR_SUCCESS;

    /* (a large denominator stays as it is: its inverse in the nested
       field would be enormous, and the fraction is a valid element) */
    if (fmpz_mpoly_q_denref(x)->length > GR_TOWER_OPTION(F->T, GR_TOWER_OPT_RATIONALIZE_LIMIT) / 8)
        return GR_SUCCESS;

    used = flint_malloc(sizeof(int) * F->cap);
    fmpz_mpoly_used_vars(used, fmpz_mpoly_q_denref(x), F->mctx);
    for (k = 1; k <= T->length; k++)
        if (used[GR_TOWER_FLAT_VAR(F, k)])
            level = k;
    flint_free(used);

    if (level == 0)
        return GR_SUCCESS;

    {
        gr_ctx_struct * Fk = gr_tower_field_at(T, level);
        gr_ptr d, e;
        fmpz_mpoly_q_t inv;

        GR_TMP_INIT2(d, e, Fk);
        fmpz_mpoly_q_init(inv, F->mctx);

        fmpz_mpoly_q_set(inv, x, F->mctx);
        fmpz_mpoly_one(fmpz_mpoly_q_numref(inv), F->mctx);
        fmpz_mpoly_swap(fmpz_mpoly_q_numref(inv), fmpz_mpoly_q_denref(inv), F->mctx);   /* inv = den */

        status = gr_tower_flat_get_nested_at(d, inv, level, F);
        if (status == GR_SUCCESS)
            status = gr_tower_inv_at(e, d, level, T);
        if (status == GR_SUCCESS)
            status = gr_tower_flat_set_nested_at(inv, e, level, F);
        if (status == GR_SUCCESS)
        {
            fmpz_mpoly_one(fmpz_mpoly_q_denref(x), F->mctx);
            fmpz_mpoly_q_mul(x, x, inv, F->mctx);
            status = gr_tower_flat_reduce(x, F);
        }

        fmpz_mpoly_q_clear(inv, F->mctx);
        GR_TMP_CLEAR2(d, e, Fk);
    }

    /* the raw fraction remains valid if the nested inversion failed */
    return GR_SUCCESS;
}

/*
    The inverse of x (reduced, integer denominator) in a number field
    tower of degree D <= 64 over QQ (no transcendental generators; the
    algebraic steps proven), by linear algebra: the multiplication matrix
    of x on the monomial basis prod a_k^(e_k) (e_k < d_k over the steps
    of degree > 1), its columns x b_j obtained from one another by
    multiplications by single generators, and the solution of M y = e_0.
    Unlike the dense inverses this needs no univariate moduli (triangular
    towers: sqrt(1 + sqrt(2)), roots of polynomials with algebraic
    coefficients). Returns 1 on success.
*/

static int
_flat_inv_linear(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, gr_tower_flat_struct * F)
{
    gr_tower_struct * T = F->T;
    slong m = 0, D = 1, k, i, j;
    slong * var, * deg, * C;
    fmpz_mpoly_q_struct * cols;
    fmpq_mat_t M, Y, B;
    fmpq * v;
    int ok = 1;

    if (!GR_TOWER_BASE_IS_CONSTS(T) || T->length < 2 || !GR_TOWER_HAS_CAP(T, GR_TOWER_CAP_RATIONAL) ||
        !fmpz_mpoly_is_fmpz(fmpz_mpoly_q_denref(x), F->mctx))
        return 0;

    for (k = 1; k <= T->length; k++)
    {
        slong d = gr_tower_step_degree(T, k);
        if (d > 1)
        {
            if (GR_TOWER_STEP(T, k - 1)->status != GR_TOWER_STATUS_PROVEN)
                return 0;
            m++;
            D *= d;
            if (D > GR_TOWER_OPTION(T, GR_TOWER_OPT_INV_LINEAR_DEGREE_LIMIT))
                return 0;
        }
    }
    if (m < 2)
        return 0;

    gr_tower_flat_ensure(F);
    var = flint_malloc(sizeof(slong) * (4 * (m + 1)));
    deg = var + m + 1;
    C = deg + m + 1;
    for (k = 1, i = 0; k <= T->length; k++)
    {
        slong d = gr_tower_step_degree(T, k);
        if (d > 1)
        {
            var[i] = GR_TOWER_FLAT_VAR(F, k);
            deg[i] = d;
            i++;
        }
    }
    C[0] = 1;
    for (i = 0; i < m; i++)
        C[i + 1] = C[i] * deg[i];

    cols = flint_malloc(sizeof(fmpz_mpoly_q_struct) * D);
    for (j = 0; j < D; j++)
        fmpz_mpoly_q_init(cols + j, F->mctx);
    fmpq_mat_init(M, D, D);
    fmpq_mat_init(Y, D, 1);
    fmpq_mat_init(B, D, 1);
    v = _fmpq_vec_init(D);

    /* columns: x b_j, with b_j = a_k b_(j - C_k) for the first k with e_k > 0 */
    fmpz_mpoly_q_set(cols + 0, x, F->mctx);
    ok = (gr_tower_flat_reduce(cols + 0, F) == GR_SUCCESS);
    for (j = 1; j < D && ok; j++)
    {
        fmpz_mpoly_q_t g;
        for (k = 0; k < m; k++)
            if ((j / C[k]) % deg[k] != 0)
                break;
        fmpz_mpoly_q_init(g, F->mctx);
        fmpz_mpoly_q_gen(g, var[k], F->mctx);
        fmpz_mpoly_q_mul(cols + j, cols + j - C[k], g, F->mctx);
        fmpz_mpoly_q_clear(g, F->mctx);
        ok = (gr_tower_flat_reduce(cols + j, F) == GR_SUCCESS);
    }

    for (j = 0; j < D && ok; j++)
    {
        ok = _gr_tower_lazy_prim_coords(v, cols + j, var, C, m, D, F->mctx);
        for (i = 0; i < D && ok; i++)
            fmpq_set(fmpq_mat_entry(M, i, j), v + i);
    }

    if (ok)
    {
        fmpq_one(fmpq_mat_entry(B, 0, 0));
        ok = fmpq_mat_solve(Y, M, B);
    }

    if (ok)
    {
        /* res = sum y_j b_j */
        fmpz_mpoly_q_t t, g;
        fmpz_mpoly_q_init(t, F->mctx);
        fmpz_mpoly_q_init(g, F->mctx);
        fmpz_mpoly_q_zero(res, F->mctx);
        for (j = 0; j < D; j++)
        {
            if (fmpq_is_zero(fmpq_mat_entry(Y, j, 0)))
                continue;
            fmpz_mpoly_q_one(t, F->mctx);
            for (k = 0; k < m; k++)
            {
                slong ek = (j / C[k]) % deg[k];
                if (ek != 0)
                {
                    fmpz_mpoly_q_gen(g, var[k], F->mctx);
                    fmpz_mpoly_pow_ui(fmpz_mpoly_q_numref(g), fmpz_mpoly_q_numref(g), ek, F->mctx);
                    fmpz_mpoly_q_mul(t, t, g, F->mctx);
                }
            }
            fmpz_mpoly_q_mul_fmpq(t, t, fmpq_mat_entry(Y, j, 0), F->mctx);
            fmpz_mpoly_q_add(res, res, t, F->mctx);
        }
        fmpz_mpoly_q_clear(t, F->mctx);
        fmpz_mpoly_q_clear(g, F->mctx);
    }

    for (j = 0; j < D; j++)
        fmpz_mpoly_q_clear(cols + j, F->mctx);
    flint_free(cols);
    fmpq_mat_clear(M);
    fmpq_mat_clear(Y);
    fmpq_mat_clear(B);
    _fmpq_vec_clear(v, D);
    flint_free(var);
    return ok;
}

/* 1/x for x = c g_1^e_1 ... g_r^e_r / d with roots of unity g_j of
   orders N_j: (d/c) g_1^(N_1 - e_1) ... (exp(-2 pi i k/p) in a cosine,
   which the inversion through the modulus Phi_p of degree p - 1 would
   otherwise compute by a resultant); returns 1 on success */
static int
_flat_inv_root_of_unity_monomial(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, gr_tower_flat_struct * F)
{
    gr_tower_struct * T = F->T;
    const fmpz_mpoly_struct * num = fmpz_mpoly_q_numref(x);
    ulong * exp;
    fmpz_t c;
    slong d, nvars;
    int ok = 1;

    if (num->length != 1 || !fmpz_mpoly_is_fmpz(fmpz_mpoly_q_denref(x), F->mctx) ||
        fmpz_mpoly_is_fmpz(num, F->mctx))
        return 0;

    nvars = F->mctx->minfo->nvars;
    exp = flint_malloc(sizeof(ulong) * nvars);
    fmpz_mpoly_get_term_exp_ui(exp, num, 0, F->mctx);
    {
        slong v, covered = 0, nonzero = 0;
        for (v = 0; v < nvars; v++)
            nonzero += (exp[v] != 0);
        for (d = 0; d < T->num_gens && ok; d++)
        {
            const gr_tower_gen_struct * g = T->gens + d;
            v = GR_TOWER_FLAT_VAR_D(F, d);
            if (exp[v] == 0)
                continue;
            covered++;
            if (g->kind != GR_TOWER_ALGEBRAIC || g->def_kind != GR_TOWER_ROOT_OF_UNITY || g->def_param <= 0)
                ok = 0;
            else
                exp[v] = ((ulong) g->def_param - exp[v] % (ulong) g->def_param) % (ulong) g->def_param;
        }
        ok = ok && (covered == nonzero);
    }

    if (ok)
    {
        fmpz_init(c);
        fmpz_mpoly_get_term_coeff_fmpz(c, num, 0, F->mctx);
        fmpz_mpoly_set(fmpz_mpoly_q_numref(res), fmpz_mpoly_q_denref(x), F->mctx);
        fmpz_mpoly_one(fmpz_mpoly_q_denref(res), F->mctx);
        {
            fmpz_mpoly_t m;
            fmpz_mpoly_init(m, F->mctx);
            fmpz_mpoly_set_coeff_ui_ui(m, 1, exp, F->mctx);
            fmpz_mpoly_mul(fmpz_mpoly_q_numref(res), fmpz_mpoly_q_numref(res), m, F->mctx);
            fmpz_mpoly_clear(m, F->mctx);
        }
        fmpz_mpoly_set_fmpz(fmpz_mpoly_q_denref(res), c, F->mctx);
        fmpz_mpoly_q_canonicalise(res, F->mctx);
        fmpz_clear(c);
        ok = (gr_tower_flat_reduce(res, F) == GR_SUCCESS);
    }

    flint_free(exp);
    return ok;
}

int
_gr_tower_lazy_inv(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    truth_t z = _gr_tower_lazy_is_zero(x, ctx);
    fmpz_mpoly_q_t t;

    if (z == T_TRUE)
        return GR_DOMAIN;
    if (z == T_UNKNOWN)
        return GR_UNABLE;

    x = _gr_tower_lazy_flat_view(x);
    fmpz_mpoly_q_init(t, x->elem.flat.mctx);
    if (!_flat_inv_root_of_unity_monomial(t, &x->elem.flat.data, x->F) &&
        !(_gr_tower_flat_dense_applicable(x->F) && _gr_tower_flat_inv_dense(t, &x->elem.flat.data, x->F)) &&
        !_flat_inv_linear(t, &x->elem.flat.data, x->F))
    {
        fmpz_mpoly_q_inv(t, &x->elem.flat.data, x->elem.flat.mctx);
        _gr_tower_lazy_rationalize(t, x->F);
    }
    _gr_tower_lazy_install(res, x->F, t, ctx);
    fmpz_mpoly_q_clear(t, res->elem.flat.mctx);
    return GR_SUCCESS;
}

/* whether x (in the current context of F) is a rational number */
static int
_flat_is_rational(const fmpz_mpoly_q_t x, gr_tower_flat_struct * F)
{
    return fmpz_mpoly_q_is_fmpq(x, F->mctx);
}

/* xx = xx / yy by the dense arithmetic, if applicable (returns 1); also
   the division by a nonzero rational number */
static int
_gr_tower_lazy_div_dense(fmpz_mpoly_q_t xx, const fmpz_mpoly_q_t yy, gr_tower_flat_struct * U)
{
    fmpz_mpoly_q_t t;
    int ok;

    if (_flat_is_rational(yy, U) && !fmpz_mpoly_q_is_zero(yy, U->mctx))
    {
        fmpz_mpoly_q_div(xx, xx, yy, U->mctx);
        return 1;
    }

    if (!_gr_tower_flat_dense_applicable(U))
        return 0;

    fmpz_mpoly_q_init(t, U->mctx);
    ok = _gr_tower_flat_inv_dense(t, yy, U) && _gr_tower_flat_mul_dense(xx, xx, t, U);
    fmpz_mpoly_q_clear(t, U->mctx);
    return ok;
}

int
_gr_tower_lazy_div(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx)
{
    gr_tower_flat_struct * U;
    fmpz_mpoly_q_t xx, yy;
    truth_t z;
    int status;

    z = _gr_tower_lazy_is_zero(y, ctx);
    if (z == T_TRUE)
        return GR_DOMAIN;
    if (z == T_UNKNOWN)
        return GR_UNABLE;

    status = _gr_tower_lazy_common(&U, xx, yy, (gr_tower_lazy_elem_struct *) x, (gr_tower_lazy_elem_struct *) y, ctx);
    if (status != GR_SUCCESS)
        return status;

    if (_gr_tower_lazy_div_dense(xx, yy, U))
    {
        /* (dense inverse and product in a single generator) */
    }
    else if (gr_tower_flat_reduce(yy, U) == GR_SUCCESS && _flat_inv_linear(yy, yy, U))
    {
        /* (the inverse by linear algebra, then the product) */
        fmpz_mpoly_q_mul(xx, xx, yy, U->mctx);
        status = gr_tower_flat_reduce(xx, U);
    }
    else
    {
        fmpz_mpoly_q_div(xx, xx, yy, U->mctx);
        status = gr_tower_flat_reduce(xx, U);
        if (status == GR_SUCCESS)
            status = _gr_tower_lazy_rationalize(xx, U);
    }
    if (status == GR_SUCCESS)
        _gr_tower_lazy_install(res, U, xx, ctx);
    fmpz_mpoly_q_clear(xx, U->mctx);
    fmpz_mpoly_q_clear(yy, U->mctx);
    return status;
}

/*
    Principal n-th root: first try to express the root in the tower of x;
    otherwise extend that tower in place.
*/
/*
    The principal n-th root of the rational number c. For c = a/b > 0 it
    is root(a b^(n-1)) / b, with the radicand factored into prime powers
    (trial division, then perfect powers of the cofactor, which is
    otherwise kept as an opaque radicand), so that the roots of different
    rational numbers live in the tower generated by the roots of their
    prime factors: sqrt(6) = sqrt(2) sqrt(3). For c < 0, the principal
    root is exp(pi i / n) times the root of |c|.
*/
