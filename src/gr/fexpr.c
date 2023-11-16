/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fexpr.h"
#include "fexpr_builtin.h"
#include "gr.h"

int
_gr_fexpr_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Symbolic expressions (fexpr)");
    return GR_SUCCESS;
}

void
_gr_fexpr_init(fexpr_t x, const gr_ctx_t ctx)
{
    fexpr_init(x);
}

void
_gr_fexpr_clear(fexpr_t x, const gr_ctx_t ctx)
{
    fexpr_clear(x);
}

void
_gr_fexpr_swap(fexpr_t x, fexpr_t y, const gr_ctx_t ctx)
{
    fexpr_t t;
    *t = *x;
    *x = *y;
    *y = *t;
}

void
_gr_fexpr_set_shallow(fexpr_t res, const fexpr_t x, const gr_ctx_t ctx)
{
    *res = *x;
}

/* todo */
int
_gr_fexpr_randtest(fexpr_t res, flint_rand_t state, const gr_ctx_t ctx)
{
    fexpr_set_si(res, -5 + (slong) n_randint(state, 10));
    return GR_SUCCESS;
}

int
_gr_fexpr_write(gr_stream_t out, const fexpr_t x, const gr_ctx_t ctx)
{
    gr_stream_write_free(out, fexpr_get_str(x));
    return GR_SUCCESS;
}

int
_gr_fexpr_zero(fexpr_t x, const gr_ctx_t ctx)
{
    fexpr_zero(x);
    return GR_SUCCESS;
}

int
_gr_fexpr_one(fexpr_t x, const gr_ctx_t ctx)
{
    fexpr_set_si(x, 1);
    return GR_SUCCESS;
}

int
_gr_fexpr_set_si(fexpr_t res, slong v, const gr_ctx_t ctx)
{
    fexpr_set_si(res, v);
    return GR_SUCCESS;
}

int
_gr_fexpr_set_ui(fexpr_t res, ulong v, const gr_ctx_t ctx)
{
    fexpr_set_ui(res, v);
    return GR_SUCCESS;
}

int
_gr_fexpr_set_fmpz(fexpr_t res, const fmpz_t v, const gr_ctx_t ctx)
{
    fexpr_set_fmpz(res, v);
    return GR_SUCCESS;
}

int
_gr_fexpr_set_fmpq(fexpr_t res, const fmpq_t v, const gr_ctx_t ctx)
{
    fexpr_set_fmpq(res, v);
    return GR_SUCCESS;
}

int
_gr_fexpr_set_other(fexpr_t res, gr_srcptr x, gr_ctx_t x_ctx, const gr_ctx_t ctx)
{
    return gr_get_fexpr(res, x, x_ctx);
}

int
_gr_fexpr_set(fexpr_t res, const fexpr_t x, const gr_ctx_t ctx)
{
    fexpr_set(res, x);
    return GR_SUCCESS;
}

int
_gr_fexpr_set_d(fexpr_t res, double x, const gr_ctx_t ctx)
{
    fexpr_set_d(res, x);
    return GR_SUCCESS;
}

int
_gr_fexpr_get_fmpz(fmpz_t res, const fexpr_t x, const gr_ctx_t ctx)
{
    if (fexpr_get_fmpz(res, x))
        return GR_SUCCESS;
    else
        return GR_UNABLE;
}

truth_t
_gr_fexpr_is_zero(const fexpr_t x, const gr_ctx_t ctx)
{
    return fexpr_is_zero(x) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fexpr_is_one(const fexpr_t x, const gr_ctx_t ctx)
{
    return fexpr_equal_ui(x, 1) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fexpr_is_neg_one(const fexpr_t x, const gr_ctx_t ctx)
{
    return fexpr_equal_si(x, -1) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fexpr_equal(const fexpr_t x, const fexpr_t y, const gr_ctx_t ctx)
{
    return fexpr_equal(x, y) ? T_TRUE : T_FALSE;
}

int
_gr_fexpr_neg(fexpr_t res, const fexpr_t x, const gr_ctx_t ctx)
{
    fexpr_neg(res, x);
    return GR_SUCCESS;
}

int
_gr_fexpr_add(fexpr_t res, const fexpr_t x, const fexpr_t y, const gr_ctx_t ctx)
{
    fexpr_add(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fexpr_sub(fexpr_t res, const fexpr_t x, const fexpr_t y, const gr_ctx_t ctx)
{
    fexpr_sub(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fexpr_mul(fexpr_t res, const fexpr_t x, const fexpr_t y, const gr_ctx_t ctx)
{
    fexpr_mul(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fexpr_inv(fexpr_t res, const fexpr_t x, const gr_ctx_t ctx)
{
    fexpr_t b;
    fexpr_init(b);
    fexpr_set_si(b, -1);
    fexpr_pow(res, x, b);
    fexpr_clear(b);
    return GR_SUCCESS;
}

int
_gr_fexpr_div(fexpr_t res, const fexpr_t x, const fexpr_t y, const gr_ctx_t ctx)
{
    fexpr_div(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fexpr_pow(fexpr_t res, const fexpr_t x, const fexpr_t y, const gr_ctx_t ctx)
{
    fexpr_pow(res, x, y);
    return GR_SUCCESS;
}

int
_gr_fexpr_pow_ui(fexpr_t res, const fexpr_t x, ulong exp, const gr_ctx_t ctx)
{
    fexpr_t b;
    fexpr_init(b);
    fexpr_set_ui(b, exp);
    fexpr_pow(res, x, b);
    fexpr_clear(b);
    return GR_SUCCESS;
}

int
_gr_fexpr_pow_si(fexpr_t res, const fexpr_t x, slong exp, const gr_ctx_t ctx)
{
    fexpr_t b;
    fexpr_init(b);
    fexpr_set_si(b, exp);
    fexpr_pow(res, x, b);
    fexpr_clear(b);
    return GR_SUCCESS;
}

int
_gr_fexpr_pow_fmpz(fexpr_t res, const fexpr_t x, const fmpz_t exp, const gr_ctx_t ctx)
{
    fexpr_t b;
    fexpr_init(b);
    fexpr_set_fmpz(b, exp);
    fexpr_pow(res, x, b);
    fexpr_clear(b);
    return GR_SUCCESS;
}

int
_gr_fexpr_pow_fmpq(fexpr_t res, const fexpr_t x, const fmpq_t exp, const gr_ctx_t ctx)
{
    fexpr_t b;
    fexpr_init(b);
    fexpr_set_fmpq(b, exp);
    fexpr_pow(res, x, b);
    fexpr_clear(b);
    return GR_SUCCESS;
}

int
_gr_fexpr_pi(fexpr_t res, const gr_ctx_t ctx)
{
    fexpr_set_symbol_builtin(res, FEXPR_Pi);
    return GR_SUCCESS;
}

int
_gr_fexpr_i(fexpr_t res, const gr_ctx_t ctx)
{
    fexpr_set_symbol_builtin(res, FEXPR_NumberI);
    return GR_SUCCESS;
}


int _fexpr_methods_initialized = 0;

gr_static_method_table _fexpr_methods;

gr_method_tab_input _fexpr_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _gr_fexpr_ctx_write},
    {GR_METHOD_INIT,            (gr_funcptr) _gr_fexpr_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) _gr_fexpr_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) _gr_fexpr_swap},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) _gr_fexpr_set_shallow},
    {GR_METHOD_RANDTEST,        (gr_funcptr) _gr_fexpr_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) _gr_fexpr_write},
    {GR_METHOD_ZERO,            (gr_funcptr) _gr_fexpr_zero},
    {GR_METHOD_ONE,             (gr_funcptr) _gr_fexpr_one},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) _gr_fexpr_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) _gr_fexpr_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) _gr_fexpr_is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) _gr_fexpr_equal},
    {GR_METHOD_SET,             (gr_funcptr) _gr_fexpr_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) _gr_fexpr_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) _gr_fexpr_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) _gr_fexpr_set_fmpz},
    {GR_METHOD_SET_FMPQ,        (gr_funcptr) _gr_fexpr_set_fmpq},
    {GR_METHOD_SET_D,           (gr_funcptr) _gr_fexpr_set_d},
    {GR_METHOD_SET_OTHER,       (gr_funcptr) _gr_fexpr_set_other},
    {GR_METHOD_GET_FMPZ,        (gr_funcptr) _gr_fexpr_get_fmpz},
    {GR_METHOD_NEG,             (gr_funcptr) _gr_fexpr_neg},
    {GR_METHOD_ADD,             (gr_funcptr) _gr_fexpr_add},
    {GR_METHOD_SUB,             (gr_funcptr) _gr_fexpr_sub},
    {GR_METHOD_MUL,             (gr_funcptr) _gr_fexpr_mul},
    {GR_METHOD_DIV,             (gr_funcptr) _gr_fexpr_div},
    {GR_METHOD_DIVEXACT,        (gr_funcptr) _gr_fexpr_div},
    {GR_METHOD_INV,             (gr_funcptr) _gr_fexpr_inv},
    {GR_METHOD_POW,             (gr_funcptr) _gr_fexpr_pow},
    {GR_METHOD_POW_UI,          (gr_funcptr) _gr_fexpr_pow_ui},
    {GR_METHOD_POW_SI,          (gr_funcptr) _gr_fexpr_pow_si},
    {GR_METHOD_POW_FMPZ,        (gr_funcptr) _gr_fexpr_pow_fmpz},
    {GR_METHOD_POW_FMPQ,        (gr_funcptr) _gr_fexpr_pow_fmpq},
    {GR_METHOD_I,               (gr_funcptr) _gr_fexpr_i},
    {GR_METHOD_PI,              (gr_funcptr) gr_not_in_domain},
    {0,                         (gr_funcptr) NULL},
};

void
gr_ctx_init_fexpr(gr_ctx_t ctx)
{
    ctx->which_ring = GR_CTX_FEXPR;
    ctx->sizeof_elem = sizeof(fexpr_struct);
    ctx->size_limit = WORD_MAX;
    ctx->methods = _fexpr_methods;

    if (!_fexpr_methods_initialized)
    {
        gr_method_tab_init(_fexpr_methods, _fexpr_methods_input);
        _fexpr_methods_initialized = 1;
    }
}
