/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "fmpz.h"
#include "fmpz_mod_mpoly_factor.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_generic.h"

typedef struct
{
    fmpz_mod_mpoly_ctx_t mctx;
    char ** vars;
}
_gr_fmpz_mod_mpoly_ctx_t;

#define MPOLYNOMIAL_CTX(ring_ctx) ((_gr_fmpz_mod_mpoly_ctx_t *)(GR_CTX_DATA_AS_PTR(ring_ctx)))
#define MPOLYNOMIAL_MCTX(ring_ctx) (MPOLYNOMIAL_CTX(ring_ctx)->mctx)

int _gr_fmpz_mod_mpoly_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Fraction field of multivariate polynomials over finite field ring (fmpz_mod) mod:");
    gr_stream_write_fmpz(out, MPOLYNOMIAL_MCTX(ctx)->ffinfo->n);
    gr_stream_write(out, " in ");
    gr_stream_write_si(out, MPOLYNOMIAL_MCTX(ctx)->minfo->nvars);
    gr_stream_write(out, " variables");
    if (MPOLYNOMIAL_MCTX(ctx)->minfo->ord == ORD_LEX)
        gr_stream_write(out, ", lex order");
    else if (MPOLYNOMIAL_MCTX(ctx)->minfo->ord == ORD_DEGLEX)
        gr_stream_write(out, ", deglex order");
    else if (MPOLYNOMIAL_MCTX(ctx)->minfo->ord == ORD_DEGREVLEX)
        gr_stream_write(out, ", degrevlex order");
    return GR_SUCCESS;
}

void _gr_fmpz_mod_mpoly_ctx_clear(gr_ctx_t ctx)
{
    if (MPOLYNOMIAL_CTX(ctx)->vars != NULL)
    {
        slong i;
        for (i = 0; i < MPOLYNOMIAL_MCTX(ctx)->minfo->nvars; i++)
            flint_free(MPOLYNOMIAL_CTX(ctx)->vars[i]);
        flint_free(MPOLYNOMIAL_CTX(ctx)->vars);
    }

    fmpz_mod_mpoly_ctx_clear(MPOLYNOMIAL_MCTX(ctx));
    flint_free(GR_CTX_DATA_AS_PTR(ctx));  
}

int
_gr_fmpz_mod_mpoly_ctx_set_gen_names(gr_ctx_t ctx, const char ** s)
{
    slong i, nvars, len;

    nvars = MPOLYNOMIAL_MCTX(ctx)->minfo->nvars;

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

void
_gr_fmpz_mod_mpoly_init(fmpz_mod_mpoly_t res, gr_ctx_t ctx)
{
    fmpz_mod_mpoly_init(res, MPOLYNOMIAL_MCTX(ctx));
}

void
_gr_fmpz_mod_mpoly_clear(fmpz_mod_mpoly_t res, gr_ctx_t ctx)
{
    fmpz_mod_mpoly_clear(res, MPOLYNOMIAL_MCTX(ctx));
}

void
_gr_fmpz_mod_mpoly_swap(fmpz_mod_mpoly_t poly1, fmpz_mod_mpoly_t poly2, gr_ctx_t ctx)
{
    fmpz_mod_mpoly_swap(poly1, poly2, MPOLYNOMIAL_MCTX(ctx));
}

void
_gr_fmpz_mod_mpoly_set_shallow(fmpz_mod_mpoly_t res, const fmpz_mod_mpoly_t poly, gr_ctx_t ctx)
{
    *res = *poly;
}

int
_gr_fmpz_mod_mpoly_randtest(fmpz_mpoly_t res, flint_rand_t state, gr_ctx_t ctx)
{
    fmpz_mod_mpoly_randtest_bits(res, state, n_randint(state, 5), 1 + n_randint(state, 3), MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_mpoly_randtest_small(fmpz_mpoly_t res, flint_rand_t state, gr_ctx_t ctx)
{
    fmpz_mod_mpoly_randtest_bits(res, state, n_randint(state, 3), 1 + n_randint(state, 3), MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

slong
_gr_fmpz_mod_mpoly_length(const fmpz_mod_mpoly_t x, gr_ctx_t ctx)
{
    return x->length;
}

int
_gr_fmpz_mod_mpoly_write(gr_stream_t out, fmpz_mod_mpoly_t poly, gr_ctx_t ctx)
{
    gr_stream_write_free(out, fmpz_mod_mpoly_get_str_pretty(poly, (const char **) MPOLYNOMIAL_CTX(ctx)->vars, MPOLYNOMIAL_MCTX(ctx)));
    return GR_SUCCESS;
}

truth_t
_gr_fmpz_mod_mpoly_equal(const fmpz_mod_mpoly_t poly1, const fmpz_mod_mpoly_t poly2, gr_ctx_t ctx)
{
    return fmpz_mod_mpoly_equal(poly1, poly2, MPOLYNOMIAL_MCTX(ctx)) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fmpz_mod_mpoly_is_zero(const fmpz_mpoly_t poly, gr_ctx_t ctx)
{
    return fmpz_mod_mpoly_is_zero(poly, MPOLYNOMIAL_MCTX(ctx)) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fmpz_mod_mpoly_is_one(const fmpz_mod_mpoly_t poly, gr_ctx_t ctx)
{
    return fmpz_mod_mpoly_is_one(poly, MPOLYNOMIAL_MCTX(ctx)) ? T_TRUE : T_FALSE;
}

int
_gr_fmpz_mod_mpoly_zero(fmpz_mod_mpoly_t res, gr_ctx_t ctx)
{
    fmpz_mod_mpoly_zero(res, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_mpoly_one(fmpz_mod_mpoly_t res, gr_ctx_t ctx)
{
    fmpz_mod_mpoly_one(res, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int _gr_fmpz_mod_mpoly_methods_initialized = 0;

gr_static_method_table _gr_fmpz_mod_mpoly_methods;

gr_method_tab_input _gr_fmpz_mod_mpoly_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,   (gr_funcptr) _gr_fmpz_mod_mpoly_ctx_write},
    {GR_METHOD_CTX_CLEAR,   (gr_funcptr) _gr_fmpz_mod_mpoly_ctx_clear},
    {GR_METHOD_CTX_IS_RING,                     (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING,         (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,          (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN,    (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FIELD,                    (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE,                   (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC,    (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_THREADSAFE,               (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_SET_GEN_NAMES,               (gr_funcptr) _gr_fmpz_mod_mpoly_ctx_set_gen_names},
    {GR_METHOD_INIT,        (gr_funcptr) _gr_fmpz_mod_mpoly_init},
    {GR_METHOD_CLEAR,       (gr_funcptr) _gr_fmpz_mod_mpoly_clear},
    {GR_METHOD_SWAP,        (gr_funcptr) _gr_fmpz_mod_mpoly_swap},
    {GR_METHOD_SET_SHALLOW, (gr_funcptr) _gr_fmpz_mod_mpoly_set_shallow},
    {GR_METHOD_RANDTEST,    (gr_funcptr) _gr_fmpz_mod_mpoly_randtest},
    {GR_METHOD_RANDTEST_SMALL,    (gr_funcptr) _gr_fmpz_mod_mpoly_randtest_small},
    {_GR_METHOD_LENGTH,     (gr_funcptr) _gr_fmpz_mod_mpoly_length},
    {GR_METHOD_WRITE,       (gr_funcptr) _gr_fmpz_mod_mpoly_write},
    {GR_METHOD_ZERO,        (gr_funcptr) _gr_fmpz_mod_mpoly_zero},
    {GR_METHOD_ONE,         (gr_funcptr) _gr_fmpz_mod_mpoly_one},
    {GR_METHOD_IS_ZERO,     (gr_funcptr) _gr_fmpz_mod_mpoly_is_zero},
    {GR_METHOD_IS_ONE,      (gr_funcptr) _gr_fmpz_mod_mpoly_is_one},
    {0,                     (gr_funcptr) NULL},
};



void
gr_ctx_init_fmpz_mpoly(gr_ctx_t ctx, slong nvars, const ordering_t ord, const fmpz *mod)
{
    ctx->which_ring = GR_CTX_FMPZ_MOD_MPOLY_Q;
    ctx->sizeof_elem = sizeof(fmpz_mod_mpoly_struct);
    GR_CTX_DATA_AS_PTR(ctx) = flint_malloc(sizeof(_gr_fmpz_mod_mpoly_ctx_t));
    ctx->size_limit = WORD_MAX;

    fmpz_mod_mpoly_ctx_init(MPOLYNOMIAL_MCTX(ctx), nvars, ord, mod);
    MPOLYNOMIAL_CTX(ctx)->vars = NULL;

    ctx->methods = _gr_fmpz_mod_mpoly_methods;

    if (!_gr_fmpz_mod_mpoly_methods_initialized)
    {
        gr_method_tab_init(_gr_fmpz_mod_mpoly_methods, _gr_fmpz_mod_mpoly_methods_input);
        _gr_fmpz_mod_mpoly_methods_initialized = 1;
    }
}