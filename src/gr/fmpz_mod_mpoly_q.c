/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "fmpz_mod_mpoly_q.h"
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

int _gr_fmpz_mod_mpoly_q_ctx_write(gr_stream_t out, gr_ctx_t ctx)
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

void _gr_fmpz_mod_mpoly_ctx_clear(gr_ctx_t ctx);
int _gr_fmpz_mod_mpoly_ctx_set_gen_names(gr_ctx_t ctx, const char ** s);

#define _gr_fmpz_mod_mpoly_q_ctx_clear _gr_fmpz_mod_mpoly_ctx_clear
#define _gr_fmpz_mod_mpoly_q_ctx_set_gen_names _gr_fmpz_mod_mpoly_ctx_set_gen_names

void
_gr_fmpz_mod_mpoly_q_init(fmpz_mod_mpoly_q_t res, gr_ctx_t ctx)
{
    fmpz_mod_mpoly_q_init(res, MPOLYNOMIAL_MCTX(ctx));
}

void
_gr_fmpz_mod_mpoly_q_clear(fmpz_mod_mpoly_q_t res, gr_ctx_t ctx)
{
    fmpz_mod_mpoly_q_clear(res, MPOLYNOMIAL_MCTX(ctx));
}



int _gr_fmpz_mod_mpoly_q_methods_initialized = 0;

gr_static_method_table _gr_fmpz_mod_mpoly_q_methods;

gr_method_tab_input _gr_fmpz_mod_mpoly_q_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,   (gr_funcptr) _gr_fmpz_mod_mpoly_q_ctx_write},
    {GR_METHOD_CTX_CLEAR,   (gr_funcptr) _gr_fmpz_mod_mpoly_q_ctx_clear},
    {GR_METHOD_CTX_IS_RING,                     (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING,         (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,          (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN,    (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FIELD,                    (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE,                   (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC,    (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_THREADSAFE,               (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_SET_GEN_NAMES,               (gr_funcptr) _gr_fmpz_mod_mpoly_q_ctx_set_gen_names},
    {GR_METHOD_INIT,        (gr_funcptr) _gr_fmpz_mod_mpoly_q_init},
    {GR_METHOD_CLEAR,       (gr_funcptr) _gr_fmpz_mod_mpoly_q_clear},


    {0,                     (gr_funcptr) NULL},
};



void
gr_ctx_init_fmpz_mpoly_q(gr_ctx_t ctx, slong nvars, const ordering_t ord, const fmpz *mod)
{
    ctx->which_ring = GR_CTX_FMPZ_MOD_MPOLY_Q;
    ctx->sizeof_elem = sizeof(fmpz_mod_mpoly_q_struct);
    GR_CTX_DATA_AS_PTR(ctx) = flint_malloc(sizeof(_gr_fmpz_mod_mpoly_ctx_t));
    ctx->size_limit = WORD_MAX;

    fmpz_mod_mpoly_ctx_init(MPOLYNOMIAL_MCTX(ctx), nvars, ord, mod);
    MPOLYNOMIAL_CTX(ctx)->vars = NULL;

    ctx->methods = _gr_fmpz_mod_mpoly_q_methods;

    if (!_gr_fmpz_mod_mpoly_q_methods_initialized)
    {
        gr_method_tab_init(_gr_fmpz_mod_mpoly_q_methods, _gr_fmpz_mod_mpoly_q_methods_input);
        _gr_fmpz_mod_mpoly_q_methods_initialized = 1;
    }
}