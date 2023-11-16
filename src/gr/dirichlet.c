/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "dirichlet.h"
#include "gr.h"

#define DIRICHLET_CTX(ctx) ((dirichlet_group_struct *) (GR_CTX_DATA_AS_PTR(ctx)))

int _gr_dirichlet_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Dirichlet group modulo ");
    gr_stream_write_ui(out, DIRICHLET_CTX(ctx)->q);
    gr_stream_write(out, " (dirichlet_char)");
    return GR_SUCCESS;
}

void
_gr_dirichlet_ctx_clear(gr_ctx_t ctx)
{
    if (DIRICHLET_CTX(ctx) != NULL) /* in case init failed */
    {
        dirichlet_group_clear(DIRICHLET_CTX(ctx));
        flint_free(DIRICHLET_CTX(ctx));
    }
}

int
_gr_dirichlet_init(dirichlet_char_t res, gr_ctx_t ctx)
{
    dirichlet_char_init(res, DIRICHLET_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_dirichlet_clear(dirichlet_char_t res, gr_ctx_t ctx)
{
    dirichlet_char_clear(res);
    return GR_SUCCESS;
}

void
_gr_dirichlet_swap(dirichlet_char_t x, dirichlet_char_t y, gr_ctx_t ctx)
{
    dirichlet_char_struct t = *x;
    *x = *y;
    *y = t;
}

void
_dirichlet_char_print(gr_stream_t out, const dirichlet_group_t G, const dirichlet_char_t x)
{
    gr_stream_write(out, "chi_");
    gr_stream_write_ui(out, G->q);
    gr_stream_write(out, "(");
    gr_stream_write_ui(out, G->q == 1 ? 1 : x->n);
    gr_stream_write(out, ", .)");
/*
    slong k;
    gr_stream_write(out, "Character ");
    gr_stream_write_ui(out, G->q == 1 ? 1 : x->n);
    gr_stream_write(out, " [");
    for (k = 0; k < G->num; k++)
    {
        gr_stream_write_ui(out, x->log[k]);
        if (k + 1 < G->num)
            gr_stream_write(out, ", ");
    }
    gr_stream_write(out, "]");
*/
}

int
_gr_dirichlet_write(gr_stream_t out, dirichlet_char_t x, gr_ctx_t ctx)
{
    _dirichlet_char_print(out, DIRICHLET_CTX(ctx), x);
    return GR_SUCCESS;
}

int
_gr_dirichlet_randtest(dirichlet_char_t res, flint_rand_t state, gr_ctx_t ctx)
{
    ulong phi, i;

    phi = DIRICHLET_CTX(ctx)->phi_q;
    i = n_randint(state, phi);

    dirichlet_char_index(res, DIRICHLET_CTX(ctx), i);
    return GR_SUCCESS;
}

truth_t
_gr_dirichlet_equal(const dirichlet_char_t x, const dirichlet_char_t y, gr_ctx_t ctx)
{
    /* fixme: arb code is inconsistent about reducing mod 1 */
    if (DIRICHLET_CTX(ctx)->q == 1)
        return T_TRUE;

    return dirichlet_char_eq(x, y) ? T_TRUE : T_FALSE;
}

int
_gr_dirichlet_set(dirichlet_char_t res, const dirichlet_char_t x, gr_ctx_t ctx)
{
    dirichlet_char_set(res, DIRICHLET_CTX(ctx), x);
    return GR_SUCCESS;
}

int
_gr_dirichlet_one(dirichlet_char_t res, gr_ctx_t ctx)
{
    dirichlet_char_one(res, DIRICHLET_CTX(ctx));
    return GR_SUCCESS;
}

truth_t
_gr_dirichlet_is_one(const dirichlet_char_t x, gr_ctx_t ctx)
{
    return dirichlet_char_is_principal(DIRICHLET_CTX(ctx), x) ? T_TRUE : T_FALSE;
}

int
_gr_dirichlet_mul(dirichlet_char_t res, const dirichlet_char_t x, const dirichlet_char_t y, gr_ctx_t ctx)
{
    dirichlet_char_mul(res, DIRICHLET_CTX(ctx), x, y);
    return GR_SUCCESS;
}

int
_gr_dirichlet_pow_ui(dirichlet_char_t res, const dirichlet_char_t x, ulong exp, gr_ctx_t ctx)
{
    dirichlet_char_pow(res, DIRICHLET_CTX(ctx), x, exp);
    return GR_SUCCESS;
}

void
_dirichlet_char_pow_fmpz(dirichlet_char_t c, const dirichlet_group_t G, const dirichlet_char_t a, const fmpz_t n)
{
    ulong k;
    ulong nred;

    for (k = 0; k < G->num ; k++)
    {
        nred = fmpz_fdiv_ui(n, G->P[k].phi.n);
        c->log[k] = nmod_mul(a->log[k], nred, G->P[k].phi);
    }

    if (fmpz_sgn(n) >= 0 && 0)  /* todo: which is faster? */
        c->n = nmod_pow_fmpz(a->n, n, G->mod);
    else
        _dirichlet_char_exp(c, G);
}

int
_gr_dirichlet_pow_fmpz(dirichlet_char_t res, const dirichlet_char_t x, const fmpz_t exp, gr_ctx_t ctx)
{
    _dirichlet_char_pow_fmpz(res, DIRICHLET_CTX(ctx), x, exp);
    return GR_SUCCESS;
}

int
_gr_dirichlet_pow_si(dirichlet_char_t res, const dirichlet_char_t x, slong exp, gr_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_set_si(t, exp);
    _dirichlet_char_pow_fmpz(res, DIRICHLET_CTX(ctx), x, t);
    fmpz_clear(t);
    return GR_SUCCESS;
}

/* todo */
int
_gr_dirichlet_inv(dirichlet_char_t res, const dirichlet_char_t x, gr_ctx_t ctx)
{
    return _gr_dirichlet_pow_si(res, x, -1, ctx);
}

/* todo: should be generic. also want left division */
int
_gr_dirichlet_div(dirichlet_char_t res, const dirichlet_char_t x, const dirichlet_char_t y, gr_ctx_t ctx)
{
    dirichlet_char_t t;
    dirichlet_char_init(t, DIRICHLET_CTX(ctx));
    _gr_dirichlet_inv(t, y, ctx);
    dirichlet_char_mul(res, DIRICHLET_CTX(ctx), x, t);
    dirichlet_char_clear(t);
    return GR_SUCCESS;
}


int _dirichlet_methods_initialized = 0;

gr_static_method_table _dirichlet_methods;

gr_method_tab_input _dirichlet_methods_input[] =
{
    {GR_METHOD_CTX_IS_FINITE,
                            (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_MULTIPLICATIVE_GROUP,
                            (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_WRITE,   (gr_funcptr) _gr_dirichlet_ctx_write},
    {GR_METHOD_CTX_CLEAR,   (gr_funcptr) _gr_dirichlet_ctx_clear},
    {GR_METHOD_INIT,        (gr_funcptr) _gr_dirichlet_init},
    {GR_METHOD_CLEAR,       (gr_funcptr) _gr_dirichlet_clear},
    {GR_METHOD_SWAP,        (gr_funcptr) _gr_dirichlet_swap},
    {GR_METHOD_RANDTEST,    (gr_funcptr) _gr_dirichlet_randtest},
    {GR_METHOD_WRITE,       (gr_funcptr) _gr_dirichlet_write},
    {GR_METHOD_ONE,         (gr_funcptr) _gr_dirichlet_one},
    {GR_METHOD_EQUAL,       (gr_funcptr) _gr_dirichlet_equal},
    {GR_METHOD_SET,         (gr_funcptr) _gr_dirichlet_set},
    {GR_METHOD_MUL,         (gr_funcptr) _gr_dirichlet_mul},
    {GR_METHOD_INV,         (gr_funcptr) _gr_dirichlet_inv},
    {GR_METHOD_DIV,         (gr_funcptr) _gr_dirichlet_div},
    {GR_METHOD_POW_UI,      (gr_funcptr) _gr_dirichlet_pow_ui},
    {GR_METHOD_POW_SI,      (gr_funcptr) _gr_dirichlet_pow_si},
    {GR_METHOD_POW_FMPZ,    (gr_funcptr) _gr_dirichlet_pow_fmpz},
    {0,                     (gr_funcptr) NULL},
};

int
gr_ctx_init_dirichlet_group(gr_ctx_t ctx, ulong q)
{
    if (q == 0)
        return GR_DOMAIN;

    ctx->which_ring = GR_CTX_DIRICHLET_GROUP;
    ctx->sizeof_elem = sizeof(dirichlet_char_struct);
    ctx->size_limit = WORD_MAX;

    GR_CTX_DATA_AS_PTR(ctx) = flint_malloc(sizeof(dirichlet_group_struct));

    if (!dirichlet_group_init(DIRICHLET_CTX(ctx), q))
    {
        flint_free(GR_CTX_DATA_AS_PTR(ctx));
        return GR_UNABLE;
    }

    ctx->methods = _dirichlet_methods;

    if (!_dirichlet_methods_initialized)
    {
        gr_method_tab_init(_dirichlet_methods, _dirichlet_methods_input);
        _dirichlet_methods_initialized = 1;
    }

    return GR_SUCCESS;
}
