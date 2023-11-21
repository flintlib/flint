/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Multivariate polynomials over generic rings */

#include "gr.h"
#include "gr_mpoly.h"

typedef struct
{
    gr_ctx_struct * base_ring;
    mpoly_ctx_t mctx;
}
_gr_gr_mpoly_ctx_t;

#define MPOLYNOMIAL_CTX(ring_ctx) ((_gr_gr_mpoly_ctx_t *)(GR_CTX_DATA_AS_PTR(ring_ctx)))
#define MPOLYNOMIAL_ELEM_CTX(ring_ctx) (MPOLYNOMIAL_CTX(ring_ctx)->base_ring)
#define MPOLYNOMIAL_MCTX(ring_ctx) (MPOLYNOMIAL_CTX(ring_ctx)->mctx)

int _gr_gr_mpoly_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Ring of multivariate polynomials over ");
    gr_ctx_write(out, MPOLYNOMIAL_ELEM_CTX(ctx));
    gr_stream_write(out, " in ");
    gr_stream_write_si(out, MPOLYNOMIAL_MCTX(ctx)->nvars);
    gr_stream_write(out, " variables");
    if (MPOLYNOMIAL_MCTX(ctx)->ord == ORD_LEX)
        gr_stream_write(out, ", lex order");
    else if (MPOLYNOMIAL_MCTX(ctx)->ord == ORD_DEGLEX)
        gr_stream_write(out, ", deglex order");
    else if (MPOLYNOMIAL_MCTX(ctx)->ord == ORD_DEGREVLEX)
        gr_stream_write(out, ", degrevlex order");
    return GR_SUCCESS;
}

void
_gr_gr_mpoly_ctx_clear(gr_ctx_t ctx)
{
    mpoly_ctx_clear(MPOLYNOMIAL_MCTX(ctx));
    flint_free(GR_CTX_DATA_AS_PTR(ctx));
}

truth_t
_gr_gr_mpoly_ctx_is_commutative_ring(gr_ctx_t ctx)
{
    return gr_ctx_is_commutative_ring(MPOLYNOMIAL_ELEM_CTX(ctx));
}

truth_t
_gr_gr_mpoly_ctx_is_integral_domain(gr_ctx_t ctx)
{
    return gr_ctx_is_integral_domain(MPOLYNOMIAL_ELEM_CTX(ctx));
}

truth_t
_gr_gr_mpoly_ctx_is_field(gr_ctx_t ctx)
{
    if (MPOLYNOMIAL_MCTX(ctx)->nvars == 0)
        return gr_ctx_is_field(MPOLYNOMIAL_ELEM_CTX(ctx));

    return T_FALSE;
}

truth_t
_gr_gr_mpoly_ctx_is_threadsafe(gr_ctx_t ctx)
{
    return gr_ctx_is_threadsafe(MPOLYNOMIAL_ELEM_CTX(ctx));
}

void
_gr_gr_mpoly_init(gr_mpoly_t res, gr_ctx_t ctx)
{
    gr_mpoly_init(res, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}

void
_gr_gr_mpoly_clear(gr_mpoly_t res, gr_ctx_t ctx)
{
    gr_mpoly_clear(res, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}

void
_gr_gr_mpoly_swap(gr_mpoly_t poly1, gr_mpoly_t poly2, gr_ctx_t ctx)
{
    gr_mpoly_swap(poly1, poly2, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}

void
_gr_gr_mpoly_set_shallow(gr_mpoly_t res, const gr_mpoly_t poly, gr_ctx_t ctx)
{
    *res = *poly;
}

int
_gr_gr_mpoly_randtest(gr_mpoly_t res, flint_rand_t state, gr_ctx_t ctx)
{
    return gr_mpoly_randtest_bits(res, state, n_randint(state, 5), 1 + n_randint(state, 3), MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}

int
_gr_gr_mpoly_write(gr_stream_t out, gr_mpoly_t poly, gr_ctx_t ctx)
{
    return gr_mpoly_write_pretty(out, poly, NULL, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}

truth_t
_gr_gr_mpoly_equal(const gr_mpoly_t poly1, const gr_mpoly_t poly2, gr_ctx_t ctx)
{
    return gr_mpoly_equal(poly1, poly2, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}

truth_t
_gr_gr_mpoly_is_zero(const gr_mpoly_t poly, gr_ctx_t ctx)
{
    return gr_mpoly_is_zero(poly, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}

truth_t
_gr_gr_mpoly_is_one(const gr_mpoly_t poly, gr_ctx_t ctx)
{
    return gr_mpoly_is_one(poly, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}

int
_gr_gr_mpoly_zero(gr_mpoly_t res, gr_ctx_t ctx)
{
    return gr_mpoly_zero(res, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}

int
_gr_gr_mpoly_one(gr_mpoly_t res, gr_ctx_t ctx)
{
    return gr_mpoly_one(res, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}

int
_gr_gr_mpoly_set(gr_mpoly_t res, const gr_mpoly_t mat, gr_ctx_t ctx)
{
    return gr_mpoly_set(res, mat, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}

int
_gr_gr_mpoly_set_si(gr_mpoly_t res, slong v, gr_ctx_t ctx)
{
    return gr_mpoly_set_si(res, v, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}

int
_gr_gr_mpoly_set_ui(gr_mpoly_t res, ulong v, gr_ctx_t ctx)
{
    return gr_mpoly_set_ui(res, v, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}

int
_gr_gr_mpoly_set_fmpz(gr_mpoly_t res, const fmpz_t v, gr_ctx_t ctx)
{
    return gr_mpoly_set_fmpz(res, v, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}

int
_gr_gr_mpoly_set_fmpq(gr_mpoly_t res, const fmpq_t v, gr_ctx_t ctx)
{
    return gr_mpoly_set_fmpq(res, v, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}

int
_gr_gr_mpoly_neg(gr_mpoly_t res, const gr_mpoly_t mat, gr_ctx_t ctx)
{
    return gr_mpoly_neg(res, mat, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}

int
_gr_gr_mpoly_add(gr_mpoly_t res, const gr_mpoly_t poly1, const gr_mpoly_t poly2, gr_ctx_t ctx)
{
    if (poly1->length + poly2->length > ctx->size_limit)
        return GR_UNABLE | gr_mpoly_zero(res, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));

    return gr_mpoly_add(res, poly1, poly2, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}

int
_gr_gr_mpoly_sub(gr_mpoly_t res, const gr_mpoly_t poly1, const gr_mpoly_t poly2, gr_ctx_t ctx)
{
    if (poly1->length + poly2->length > ctx->size_limit)
        return GR_UNABLE | gr_mpoly_zero(res, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));

    return gr_mpoly_sub(res, poly1, poly2, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}

int
_gr_gr_mpoly_mul(gr_mpoly_t res, const gr_mpoly_t poly1, const gr_mpoly_t poly2, gr_ctx_t ctx)
{
    if (poly1->length * poly2->length > ctx->size_limit)
        return GR_UNABLE | gr_mpoly_zero(res, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));

    return gr_mpoly_mul(res, poly1, poly2, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}


int _gr__gr_gr_mpoly_methods_initialized = 0;

gr_static_method_table _gr__gr_gr_mpoly_methods;

gr_method_tab_input _gr__gr_gr_mpoly_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,   (gr_funcptr) _gr_gr_mpoly_ctx_write},
    {GR_METHOD_CTX_CLEAR,   (gr_funcptr) _gr_gr_mpoly_ctx_clear},
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) gr_generic_ctx_predicate_true},  /* todo: matrices over semirings? */
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) _gr_gr_mpoly_ctx_is_commutative_ring},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) _gr_gr_mpoly_ctx_is_integral_domain},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) _gr_gr_mpoly_ctx_is_field},
    {GR_METHOD_CTX_IS_THREADSAFE,       (gr_funcptr) _gr_gr_mpoly_ctx_is_threadsafe},
    {GR_METHOD_INIT,        (gr_funcptr) _gr_gr_mpoly_init},
    {GR_METHOD_CLEAR,       (gr_funcptr) _gr_gr_mpoly_clear},
    {GR_METHOD_SWAP,        (gr_funcptr) _gr_gr_mpoly_swap},
    {GR_METHOD_SET_SHALLOW, (gr_funcptr) _gr_gr_mpoly_set_shallow},
    {GR_METHOD_RANDTEST,    (gr_funcptr) _gr_gr_mpoly_randtest},
    {GR_METHOD_WRITE,       (gr_funcptr) _gr_gr_mpoly_write},
    {GR_METHOD_ZERO,        (gr_funcptr) _gr_gr_mpoly_zero},
    {GR_METHOD_ONE,         (gr_funcptr) _gr_gr_mpoly_one},
    {GR_METHOD_IS_ZERO,     (gr_funcptr) _gr_gr_mpoly_is_zero},
    {GR_METHOD_IS_ONE,      (gr_funcptr) _gr_gr_mpoly_is_one},
    {GR_METHOD_EQUAL,       (gr_funcptr) _gr_gr_mpoly_equal},
    {GR_METHOD_SET,         (gr_funcptr) _gr_gr_mpoly_set},
    {GR_METHOD_SET_UI,      (gr_funcptr) _gr_gr_mpoly_set_ui},
    {GR_METHOD_SET_SI,      (gr_funcptr) _gr_gr_mpoly_set_si},
    {GR_METHOD_SET_FMPZ,    (gr_funcptr) _gr_gr_mpoly_set_fmpz},
    {GR_METHOD_SET_FMPQ,    (gr_funcptr) _gr_gr_mpoly_set_fmpq},
    {GR_METHOD_NEG,         (gr_funcptr) _gr_gr_mpoly_neg},
    {GR_METHOD_ADD,         (gr_funcptr) _gr_gr_mpoly_add},
    {GR_METHOD_SUB,         (gr_funcptr) _gr_gr_mpoly_sub},
    {GR_METHOD_MUL,         (gr_funcptr) _gr_gr_mpoly_mul},
    {0,                     (gr_funcptr) NULL},
};

void
gr_ctx_init_gr_mpoly(gr_ctx_t ctx, gr_ctx_t base_ring, slong nvars, const ordering_t ord)
{
    ctx->which_ring = GR_CTX_GR_MPOLY;
    ctx->sizeof_elem = sizeof(gr_mpoly_struct);
    GR_CTX_DATA_AS_PTR(ctx) = flint_malloc(sizeof(_gr_gr_mpoly_ctx_t));
    ctx->size_limit = WORD_MAX;

    MPOLYNOMIAL_ELEM_CTX(ctx) = base_ring;
    mpoly_ctx_init(MPOLYNOMIAL_MCTX(ctx), nvars, ord);

    ctx->methods = _gr__gr_gr_mpoly_methods;

    if (!_gr__gr_gr_mpoly_methods_initialized)
    {
        gr_method_tab_init(_gr__gr_gr_mpoly_methods, _gr__gr_gr_mpoly_methods_input);
        _gr__gr_gr_mpoly_methods_initialized = 1;
    }
}
