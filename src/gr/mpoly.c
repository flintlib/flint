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
    gr_ctx_struct * base_ring;
    mpoly_ctx_t mctx;
    char ** vars;
}
_gr_gr_mpoly_ctx_t;

#define MPOLYNOMIAL_CTX(ring_ctx) ((_gr_gr_mpoly_ctx_t *)(GR_CTX_DATA_AS_PTR(ring_ctx)))
#define MPOLYNOMIAL_ELEM_CTX(ring_ctx) (MPOLYNOMIAL_CTX(ring_ctx)->base_ring)
#define MPOLYNOMIAL_MCTX(ring_ctx) (MPOLYNOMIAL_CTX(ring_ctx)->mctx)

static int _gr_gr_mpoly_ctx_write(gr_stream_t out, gr_ctx_t ctx)
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

static void
_gr_gr_mpoly_ctx_clear(gr_ctx_t ctx)
{
    if (MPOLYNOMIAL_CTX(ctx)->vars != NULL)
    {
        slong i;
        for (i = 0; i < MPOLYNOMIAL_MCTX(ctx)->nvars; i++)
            flint_free(MPOLYNOMIAL_CTX(ctx)->vars[i]);
        flint_free(MPOLYNOMIAL_CTX(ctx)->vars);
    }

    mpoly_ctx_clear(MPOLYNOMIAL_MCTX(ctx));
    flint_free(GR_CTX_DATA_AS_PTR(ctx));
}

static int
_gr_gr_mpoly_ctx_set_gen_names(gr_ctx_t ctx, const char ** s)
{
    slong i, nvars, len;

    nvars = MPOLYNOMIAL_MCTX(ctx)->nvars;

    if (MPOLYNOMIAL_CTX(ctx)->vars == NULL)
    {
        MPOLYNOMIAL_CTX(ctx)->vars = flint_malloc(nvars * sizeof(char *));
        for (i = 0; i < nvars; i++)
            MPOLYNOMIAL_CTX(ctx)->vars[i] = NULL;
    }
    else
    {
        for (i = 0; i < nvars; i++)
            flint_free(MPOLYNOMIAL_CTX(ctx)->vars[i]);
    }

    for (i = 0; i < nvars; i++)
    {
        len = strlen(s[i]);
        MPOLYNOMIAL_CTX(ctx)->vars[i] = flint_realloc(MPOLYNOMIAL_CTX(ctx)->vars[i], len + 1);
        memcpy(MPOLYNOMIAL_CTX(ctx)->vars[i], s[i], len + 1);
    }

    return GR_SUCCESS;
}

/* todo: everything is a ring when there are 0 vars? */
static truth_t
_gr_gr_mpoly_ctx_is_ring(gr_ctx_t ctx)
{
    return gr_ctx_is_ring(MPOLYNOMIAL_ELEM_CTX(ctx));
}

static truth_t
_gr_gr_mpoly_ctx_is_commutative_ring(gr_ctx_t ctx)
{
    return gr_ctx_is_commutative_ring(MPOLYNOMIAL_ELEM_CTX(ctx));
}

static truth_t
_gr_gr_mpoly_ctx_is_integral_domain(gr_ctx_t ctx)
{
    return gr_ctx_is_integral_domain(MPOLYNOMIAL_ELEM_CTX(ctx));
}

static truth_t
_gr_gr_mpoly_ctx_is_field(gr_ctx_t ctx)
{
    if (MPOLYNOMIAL_MCTX(ctx)->nvars == 0)
        return gr_ctx_is_field(MPOLYNOMIAL_ELEM_CTX(ctx));

    return T_FALSE;
}

static truth_t
_gr_gr_mpoly_ctx_is_threadsafe(gr_ctx_t ctx)
{
    return gr_ctx_is_threadsafe(MPOLYNOMIAL_ELEM_CTX(ctx));
}

static void
_gr_gr_mpoly_init(gr_mpoly_t res, gr_ctx_t ctx)
{
    gr_mpoly_init(res, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}

static void
_gr_gr_mpoly_clear(gr_mpoly_t res, gr_ctx_t ctx)
{
    gr_mpoly_clear(res, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}

static void
_gr_gr_mpoly_swap(gr_mpoly_t poly1, gr_mpoly_t poly2, gr_ctx_t ctx)
{
    gr_mpoly_swap(poly1, poly2, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}

static void
_gr_gr_mpoly_set_shallow(gr_mpoly_t res, const gr_mpoly_t poly, gr_ctx_t FLINT_UNUSED(ctx))
{
    *res = *poly;
}

static int
_gr_gr_mpoly_randtest(gr_mpoly_t res, flint_rand_t state, gr_ctx_t ctx)
{
    return gr_mpoly_randtest_bits(res, state, n_randint(state, 5), 1 + n_randint(state, 3), MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}

static slong
_gr_gr_mpoly_length(const gr_mpoly_t x, gr_ctx_t FLINT_UNUSED(ctx))
{
    return x->length;
}

static int
_gr_gr_mpoly_write(gr_stream_t out, gr_mpoly_t poly, gr_ctx_t ctx)
{
    return gr_mpoly_write_pretty(out, poly, (const char **) MPOLYNOMIAL_CTX(ctx)->vars, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}

static truth_t
_gr_gr_mpoly_equal(const gr_mpoly_t poly1, const gr_mpoly_t poly2, gr_ctx_t ctx)
{
    return gr_mpoly_equal(poly1, poly2, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}

static truth_t
_gr_gr_mpoly_is_zero(const gr_mpoly_t poly, gr_ctx_t ctx)
{
    return gr_mpoly_is_zero(poly, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}

static truth_t
_gr_gr_mpoly_is_one(const gr_mpoly_t poly, gr_ctx_t ctx)
{
    return gr_mpoly_is_one(poly, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}

static int
_gr_gr_mpoly_gens(gr_vec_t res, gr_ctx_t ctx)
{
    slong i, n;
    int status = GR_SUCCESS;

    n = MPOLYNOMIAL_MCTX(ctx)->nvars;

    gr_vec_set_length(res, n, ctx);
    for (i = 0; i < n; i++)
        status |= gr_mpoly_gen(((gr_mpoly_struct *) res->entries) + i, i, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));

    return status;
}

static int
_gr_gr_mpoly_gens_recursive(gr_vec_t vec, gr_ctx_t ctx)
{
    int status;
    gr_vec_t vec1;
    slong i, n, m;

    /* Get generators of the element ring */
    gr_vec_init(vec1, 0, MPOLYNOMIAL_ELEM_CTX(ctx));
    status = gr_gens_recursive(vec1, MPOLYNOMIAL_ELEM_CTX(ctx));
    n = vec1->length;

    m = MPOLYNOMIAL_MCTX(ctx)->nvars;

    gr_vec_set_length(vec, n + m, ctx);

    /* Promote to polynomials */
    for (i = 0; i < n; i++)
    {
        status |= gr_mpoly_set_scalar(gr_vec_entry_ptr(vec, i, ctx),
                gr_vec_entry_srcptr(vec1, i, MPOLYNOMIAL_ELEM_CTX(ctx)),
                    MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
    }

    for (i = 0; i < m; i++)
    {
        status |= gr_mpoly_gen(((gr_mpoly_struct *) vec->entries) + n + i,
            i, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
    }

    gr_vec_clear(vec1, MPOLYNOMIAL_ELEM_CTX(ctx));

    return status;
}

static int
_gr_gr_mpoly_zero(gr_mpoly_t res, gr_ctx_t ctx)
{
    return gr_mpoly_zero(res, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}

static int
_gr_gr_mpoly_one(gr_mpoly_t res, gr_ctx_t ctx)
{
    return gr_mpoly_one(res, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}

static int
_gr_gr_mpoly_set(gr_mpoly_t res, const gr_mpoly_t mat, gr_ctx_t ctx)
{
    return gr_mpoly_set(res, mat, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}

static int
_gr_gr_mpoly_set_si(gr_mpoly_t res, slong v, gr_ctx_t ctx)
{
    return gr_mpoly_set_si(res, v, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}

static int
_gr_gr_mpoly_set_ui(gr_mpoly_t res, ulong v, gr_ctx_t ctx)
{
    return gr_mpoly_set_ui(res, v, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}

static int
_gr_gr_mpoly_set_fmpz(gr_mpoly_t res, const fmpz_t v, gr_ctx_t ctx)
{
    return gr_mpoly_set_fmpz(res, v, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}

static int
_gr_gr_mpoly_set_fmpq(gr_mpoly_t res, const fmpq_t v, gr_ctx_t ctx)
{
    return gr_mpoly_set_fmpq(res, v, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}

static int
_gr_gr_mpoly_neg(gr_mpoly_t res, const gr_mpoly_t mat, gr_ctx_t ctx)
{
    return gr_mpoly_neg(res, mat, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}

static int
_gr_gr_mpoly_add(gr_mpoly_t res, const gr_mpoly_t poly1, const gr_mpoly_t poly2, gr_ctx_t ctx)
{
    if ((ulong) (poly1->length + poly2->length) > ctx->size_limit)
        return GR_UNABLE | gr_mpoly_zero(res, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));

    return gr_mpoly_add(res, poly1, poly2, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}

static int
_gr_gr_mpoly_sub(gr_mpoly_t res, const gr_mpoly_t poly1, const gr_mpoly_t poly2, gr_ctx_t ctx)
{
    if ((ulong) (poly1->length + poly2->length) > ctx->size_limit)
        return GR_UNABLE | gr_mpoly_zero(res, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));

    return gr_mpoly_sub(res, poly1, poly2, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}

static int
_gr_gr_mpoly_mul(gr_mpoly_t res, const gr_mpoly_t poly1, const gr_mpoly_t poly2, gr_ctx_t ctx)
{
    if ((ulong) (poly1->length * poly2->length) > ctx->size_limit)
        return GR_UNABLE | gr_mpoly_zero(res, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));

    return gr_mpoly_mul(res, poly1, poly2, MPOLYNOMIAL_MCTX(ctx), MPOLYNOMIAL_ELEM_CTX(ctx));
}


int _gr__gr_gr_mpoly_methods_initialized = 0;

gr_static_method_table _gr__gr_gr_mpoly_methods;

gr_method_tab_input _gr__gr_gr_mpoly_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,   (gr_funcptr) _gr_gr_mpoly_ctx_write},
    {GR_METHOD_CTX_CLEAR,   (gr_funcptr) _gr_gr_mpoly_ctx_clear},
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) _gr_gr_mpoly_ctx_is_ring},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) _gr_gr_mpoly_ctx_is_commutative_ring},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) _gr_gr_mpoly_ctx_is_integral_domain},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) _gr_gr_mpoly_ctx_is_field},
    {GR_METHOD_CTX_IS_THREADSAFE,       (gr_funcptr) _gr_gr_mpoly_ctx_is_threadsafe},
    {GR_METHOD_CTX_SET_GEN_NAMES,       (gr_funcptr) _gr_gr_mpoly_ctx_set_gen_names},
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

void
gr_ctx_init_gr_mpoly(gr_ctx_t ctx, gr_ctx_t base_ring, slong nvars, const ordering_t ord)
{
    ctx->which_ring = GR_CTX_GR_MPOLY;
    ctx->sizeof_elem = sizeof(gr_mpoly_struct);
    GR_CTX_DATA_AS_PTR(ctx) = flint_malloc(sizeof(_gr_gr_mpoly_ctx_t));
    ctx->size_limit = WORD_MAX;

    MPOLYNOMIAL_ELEM_CTX(ctx) = base_ring;
    mpoly_ctx_init(MPOLYNOMIAL_MCTX(ctx), nvars, ord);
    MPOLYNOMIAL_CTX(ctx)->vars = NULL;

    ctx->methods = _gr__gr_gr_mpoly_methods;

    if (!_gr__gr_gr_mpoly_methods_initialized)
    {
        gr_method_tab_init(_gr__gr_gr_mpoly_methods, _gr__gr_gr_mpoly_methods_input);
        _gr__gr_gr_mpoly_methods_initialized = 1;
    }
}
