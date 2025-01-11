/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Multivariate polynomials over generic rings */

#include <string.h>
#include "mpoly.h"
#include "gr.h"
#include "gr_mpoly.h"
#include "gr_generic.h"

typedef struct
{
    gr_ctx_struct * cctx;
    mpoly_ctx_struct * mctx;
    char ** vars;
} _gr_mpoly_ctx_struct;

typedef gr_ctx_struct gr_mpoly_ctx_struct;
typedef gr_mpoly_ctx_struct gr_mpoly_ctx_t[1];

#define GR_MPOLY_MCTX(ctx) (((_gr_mpoly_ctx_struct *) (ctx->data))->mctx)
#define GR_MPOLY_CCTX(ctx) (((_gr_mpoly_ctx_struct *) (ctx->data))->cctx)
#define GR_MPOLY_VARS(ctx) (((_gr_mpoly_ctx_struct *) (ctx->data))->vars)
#define GR_MPOLY_NVARS(ctx) (GR_MPOLY_MCTX(ctx)->nvars)


int gr_mpoly_ctx_write(gr_stream_t out, gr_mpoly_ctx_t ctx)
{
    gr_stream_write(out, "Ring of multivariate polynomials over ");
    gr_ctx_write(out, GR_MPOLY_CCTX(ctx));
    gr_stream_write(out, " in ");
    gr_stream_write_si(out, GR_MPOLY_NVARS(ctx));
    gr_stream_write(out, " variables");
    if (GR_MPOLY_MCTX(ctx)->ord == ORD_LEX)
        gr_stream_write(out, ", lex order");
    else if (GR_MPOLY_MCTX(ctx)->ord == ORD_DEGLEX)
        gr_stream_write(out, ", deglex order");
    else if (GR_MPOLY_MCTX(ctx)->ord == ORD_DEGREVLEX)
        gr_stream_write(out, ", degrevlex order");
    return GR_SUCCESS;
}

void
gr_mpoly_ctx_clear(gr_mpoly_ctx_t ctx)
{
    if (GR_MPOLY_VARS(ctx) != NULL)
    {
        slong i;
        for (i = 0; i < GR_MPOLY_NVARS(ctx); i++)
            flint_free(GR_MPOLY_VARS(ctx)[i]);
        flint_free(GR_MPOLY_VARS(ctx));
    }

    mpoly_ctx_clear(GR_MPOLY_MCTX(ctx));
    flint_free(GR_MPOLY_MCTX(ctx));
}

int
gr_mpoly_ctx_set_gen_names(gr_mpoly_ctx_t ctx, const char ** s)
{
    slong i, nvars, len;

    nvars = GR_MPOLY_NVARS(ctx);

    if (GR_MPOLY_VARS(ctx) == NULL)
    {
        GR_MPOLY_VARS(ctx) = flint_malloc(nvars * sizeof(char *));
        for (i = 0; i < nvars; i++)
            GR_MPOLY_VARS(ctx)[i] = NULL;
    }

    for (i = 0; i < nvars; i++)
    {
        len = strlen(s[i]);
        GR_MPOLY_VARS(ctx)[i] = flint_realloc(GR_MPOLY_VARS(ctx)[i], len + 1);
        memcpy(GR_MPOLY_VARS(ctx)[i], s[i], len + 1);
    }

    return GR_SUCCESS;
}

truth_t
gr_mpoly_ctx_is_ring(gr_mpoly_ctx_t ctx)
{
    return gr_ctx_is_ring(GR_MPOLY_CCTX(ctx));
}

truth_t
gr_mpoly_ctx_is_zero_ring(gr_mpoly_ctx_t ctx)
{
    return gr_ctx_is_zero_ring(GR_MPOLY_CCTX(ctx));
}

truth_t
gr_mpoly_ctx_is_commutative_ring(gr_mpoly_ctx_t ctx)
{
    return gr_ctx_is_commutative_ring(GR_MPOLY_CCTX(ctx));
}

truth_t
gr_mpoly_ctx_is_integral_domain(gr_mpoly_ctx_t ctx)
{
    return gr_ctx_is_integral_domain(GR_MPOLY_CCTX(ctx));
}

truth_t
gr_mpoly_ctx_is_field(gr_mpoly_ctx_t ctx)
{
    if (GR_MPOLY_NVARS(ctx) == 0)
        return gr_ctx_is_field(GR_MPOLY_CCTX(ctx));
    else
        return T_FALSE;
}

truth_t
gr_mpoly_ctx_is_threadsafe(gr_mpoly_ctx_t ctx)
{
    return gr_ctx_is_threadsafe(GR_MPOLY_CCTX(ctx));
}

void
_gr_gr_mpoly_init(gr_mpoly_t res, gr_mpoly_ctx_t ctx)
{
    gr_mpoly_init(res, GR_MPOLY_MCTX(ctx), GR_MPOLY_CCTX(ctx));
}

void
_gr_gr_mpoly_clear(gr_mpoly_t res, gr_mpoly_ctx_t ctx)
{
    gr_mpoly_clear(res, GR_MPOLY_MCTX(ctx), GR_MPOLY_CCTX(ctx));
}

void
_gr_gr_mpoly_swap(gr_mpoly_t poly1, gr_mpoly_t poly2, gr_mpoly_ctx_t ctx)
{
    gr_mpoly_swap(poly1, poly2, GR_MPOLY_MCTX(ctx), GR_MPOLY_CCTX(ctx));
}

void
_gr_gr_mpoly_set_shallow(gr_mpoly_t res, const gr_mpoly_t poly, gr_mpoly_ctx_t ctx)
{
    *res = *poly;
}

int
_gr_gr_mpoly_randtest(gr_mpoly_t res, flint_rand_t state, gr_mpoly_ctx_t ctx)
{
    return gr_mpoly_randtest_bits(res, state, n_randint(state, 5), 1 + n_randint(state, 3), GR_MPOLY_MCTX(ctx), GR_MPOLY_CCTX(ctx));
}

slong
_gr_gr_mpoly_length(const gr_mpoly_t x, gr_mpoly_ctx_t ctx)
{
    return x->length;
}

int
_gr_gr_mpoly_write(gr_stream_t out, gr_mpoly_t poly, gr_mpoly_ctx_t ctx)
{
    return gr_mpoly_write_pretty(out, poly, (const char **) GR_MPOLY_VARS(ctx), GR_MPOLY_MCTX(ctx), GR_MPOLY_CCTX(ctx));
}

truth_t
_gr_gr_mpoly_equal(const gr_mpoly_t poly1, const gr_mpoly_t poly2, gr_mpoly_ctx_t ctx)
{
    return gr_mpoly_equal(poly1, poly2, GR_MPOLY_MCTX(ctx), GR_MPOLY_CCTX(ctx));
}

truth_t
_gr_gr_mpoly_is_zero(const gr_mpoly_t poly, gr_mpoly_ctx_t ctx)
{
    return gr_mpoly_is_zero(poly, GR_MPOLY_MCTX(ctx), GR_MPOLY_CCTX(ctx));
}

truth_t
_gr_gr_mpoly_is_one(const gr_mpoly_t poly, gr_mpoly_ctx_t ctx)
{
    return gr_mpoly_is_one(poly, GR_MPOLY_MCTX(ctx), GR_MPOLY_CCTX(ctx));
}

int
_gr_gr_mpoly_gens(gr_vec_t res, gr_mpoly_ctx_t ctx)
{
    slong i, n;
    int status = GR_SUCCESS;

    n = GR_MPOLY_NVARS(ctx);

    gr_vec_set_length(res, n, ctx);
    for (i = 0; i < n; i++)
        status |= gr_mpoly_gen(((gr_mpoly_struct *) res->entries) + i, i, GR_MPOLY_MCTX(ctx), GR_MPOLY_CCTX(ctx));

    return status;
}

int
_gr_gr_mpoly_gens_recursive(gr_vec_t vec, gr_mpoly_ctx_t ctx)
{
    int status;
    gr_vec_t vec1;
    slong i, n, m;

    /* Get generators of the element ring */
    gr_vec_init(vec1, 0, GR_MPOLY_CCTX(ctx));
    status = gr_gens_recursive(vec1, GR_MPOLY_CCTX(ctx));
    n = vec1->length;

    m = GR_MPOLY_NVARS(ctx);

    gr_vec_set_length(vec, n + m, ctx);

    /* Promote to polynomials */
    for (i = 0; i < n; i++)
    {
        status |= gr_mpoly_set_scalar(gr_vec_entry_ptr(vec, i, ctx),
                gr_vec_entry_srcptr(vec1, i, GR_MPOLY_CCTX(ctx)),
                    GR_MPOLY_MCTX(ctx), GR_MPOLY_CCTX(ctx));
    }

    for (i = 0; i < m; i++)
    {
        status |= gr_mpoly_gen(((gr_mpoly_struct *) vec->entries) + n + i,
            i, GR_MPOLY_MCTX(ctx), GR_MPOLY_CCTX(ctx));
    }

    gr_vec_clear(vec1, GR_MPOLY_CCTX(ctx));

    return status;
}

int
_gr_gr_mpoly_zero(gr_mpoly_t res, gr_mpoly_ctx_t ctx)
{
    return gr_mpoly_zero(res, GR_MPOLY_MCTX(ctx), GR_MPOLY_CCTX(ctx));
}

int
_gr_gr_mpoly_one(gr_mpoly_t res, gr_mpoly_ctx_t ctx)
{
    return gr_mpoly_one(res, GR_MPOLY_MCTX(ctx), GR_MPOLY_CCTX(ctx));
}

int
_gr_gr_mpoly_set(gr_mpoly_t res, const gr_mpoly_t mat, gr_mpoly_ctx_t ctx)
{
    return gr_mpoly_set(res, mat, GR_MPOLY_MCTX(ctx), GR_MPOLY_CCTX(ctx));
}

int
_gr_gr_mpoly_set_si(gr_mpoly_t res, slong v, gr_mpoly_ctx_t ctx)
{
    return gr_mpoly_set_si(res, v, GR_MPOLY_MCTX(ctx), GR_MPOLY_CCTX(ctx));
}

int
_gr_gr_mpoly_set_ui(gr_mpoly_t res, ulong v, gr_mpoly_ctx_t ctx)
{
    return gr_mpoly_set_ui(res, v, GR_MPOLY_MCTX(ctx), GR_MPOLY_CCTX(ctx));
}

int
_gr_gr_mpoly_set_fmpz(gr_mpoly_t res, const fmpz_t v, gr_mpoly_ctx_t ctx)
{
    return gr_mpoly_set_fmpz(res, v, GR_MPOLY_MCTX(ctx), GR_MPOLY_CCTX(ctx));
}

int
_gr_gr_mpoly_set_fmpq(gr_mpoly_t res, const fmpq_t v, gr_mpoly_ctx_t ctx)
{
    return gr_mpoly_set_fmpq(res, v, GR_MPOLY_MCTX(ctx), GR_MPOLY_CCTX(ctx));
}

int
_gr_gr_mpoly_neg(gr_mpoly_t res, const gr_mpoly_t mat, gr_mpoly_ctx_t ctx)
{
    return gr_mpoly_neg(res, mat, GR_MPOLY_MCTX(ctx), GR_MPOLY_CCTX(ctx));
}

int
_gr_gr_mpoly_add(gr_mpoly_t res, const gr_mpoly_t poly1, const gr_mpoly_t poly2, gr_mpoly_ctx_t ctx)
{
    if (poly1->length + poly2->length > ctx->size_limit)
        return GR_UNABLE | gr_mpoly_zero(res, GR_MPOLY_MCTX(ctx), GR_MPOLY_CCTX(ctx));

    return gr_mpoly_add(res, poly1, poly2, GR_MPOLY_MCTX(ctx), GR_MPOLY_CCTX(ctx));
}

int
_gr_gr_mpoly_sub(gr_mpoly_t res, const gr_mpoly_t poly1, const gr_mpoly_t poly2, gr_mpoly_ctx_t ctx)
{
    if (poly1->length + poly2->length > ctx->size_limit)
        return GR_UNABLE | gr_mpoly_zero(res, GR_MPOLY_MCTX(ctx), GR_MPOLY_CCTX(ctx));

    return gr_mpoly_sub(res, poly1, poly2, GR_MPOLY_MCTX(ctx), GR_MPOLY_CCTX(ctx));
}

int
_gr_gr_mpoly_mul(gr_mpoly_t res, const gr_mpoly_t poly1, const gr_mpoly_t poly2, gr_mpoly_ctx_t ctx)
{
    if (poly1->length * poly2->length > ctx->size_limit)
        return GR_UNABLE | gr_mpoly_zero(res, GR_MPOLY_MCTX(ctx), GR_MPOLY_CCTX(ctx));

    return gr_mpoly_mul(res, poly1, poly2, GR_MPOLY_MCTX(ctx), GR_MPOLY_CCTX(ctx));
}


int _gr_mpoly_methods_initialized = 0;

gr_static_method_table _gr_mpoly_methods;

gr_method_tab_input _gr_mpoly_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,   (gr_funcptr) gr_mpoly_ctx_write},
    {GR_METHOD_CTX_CLEAR,   (gr_funcptr) gr_mpoly_ctx_clear},
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) gr_mpoly_ctx_is_ring},
    {GR_METHOD_CTX_IS_ZERO_RING,     (gr_funcptr) gr_mpoly_ctx_is_zero_ring},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) gr_mpoly_ctx_is_commutative_ring},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) gr_mpoly_ctx_is_integral_domain},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) gr_mpoly_ctx_is_field},
    {GR_METHOD_CTX_IS_THREADSAFE,       (gr_funcptr) gr_mpoly_ctx_is_threadsafe},
    {GR_METHOD_CTX_SET_GEN_NAMES,       (gr_funcptr) gr_mpoly_ctx_set_gen_names},
    {GR_METHOD_INIT,        (gr_funcptr) _gr_gr_mpoly_init},
    {GR_METHOD_CLEAR,       (gr_funcptr) _gr_gr_mpoly_clear},
    {GR_METHOD_SWAP,        (gr_funcptr) _gr_gr_mpoly_swap},
    {GR_METHOD_SET_SHALLOW, (gr_funcptr) _gr_gr_mpoly_set_shallow},
    {GR_METHOD_RANDTEST,    (gr_funcptr) _gr_gr_mpoly_randtest},
    {_GR_METHOD_LENGTH,     (gr_funcptr) _gr_gr_mpoly_length},
    {GR_METHOD_WRITE,       (gr_funcptr) _gr_gr_mpoly_write},
    {GR_METHOD_GENS,        (gr_funcptr) _gr_gr_mpoly_gens},
    {GR_METHOD_GENS_RECURSIVE,       (gr_funcptr) _gr_gr_mpoly_gens_recursive},
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
    {GR_METHOD_SET_STR,     (gr_funcptr) gr_generic_set_str_balance_additions},
    {GR_METHOD_NEG,         (gr_funcptr) _gr_gr_mpoly_neg},
    {GR_METHOD_ADD,         (gr_funcptr) _gr_gr_mpoly_add},
    {GR_METHOD_SUB,         (gr_funcptr) _gr_gr_mpoly_sub},
    {GR_METHOD_MUL,         (gr_funcptr) _gr_gr_mpoly_mul},
    {0,                     (gr_funcptr) NULL},
};

/* todo: first arg as gr_mpoly_ctx_t */
void
gr_ctx_init_gr_mpoly(gr_ctx_t ctx, gr_ctx_t base_ring, slong nvars, const ordering_t ord)
{
    ctx->which_ring = GR_CTX_GR_MPOLY;
    ctx->sizeof_elem = sizeof(gr_mpoly_struct);
    ctx->size_limit = WORD_MAX;

    /* by reference */
    GR_MPOLY_CCTX(ctx) = base_ring;

    /* allocated here */
    GR_MPOLY_MCTX(ctx) = flint_malloc(sizeof(mpoly_ctx_struct));
    mpoly_ctx_init(GR_MPOLY_MCTX(ctx), nvars, ord);

    GR_MPOLY_VARS(ctx) = NULL;

    ctx->methods = _gr_mpoly_methods;

    if (!_gr_mpoly_methods_initialized)
    {
        gr_method_tab_init(_gr_mpoly_methods, _gr_mpoly_methods_input);
        _gr_mpoly_methods_initialized = 1;
    }
}
