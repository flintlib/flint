/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
//#include "fmpq.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "gr.h"
#include "gr_mat.h"

#define NMOD8_CTX_REF(ring_ctx) (((nmod_t *)((ring_ctx))))
#define NMOD8_CTX(ring_ctx) (*NMOD8_CTX_REF(ring_ctx))

typedef unsigned char nmod8_struct;
typedef nmod8_struct nmod8_t[1];

void
nmod8_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Integers mod ");
    gr_stream_write_ui(out, NMOD8_CTX(ctx).n);
    gr_stream_write(out, " (nmod8)");
}

/* todo: n_is_prime is fast, but this should still be cached
   or use a fixed table lookup */
truth_t
nmod8_ctx_is_field(const gr_ctx_t ctx)
{
    return n_is_prime(NMOD8_CTX(ctx).n) ? T_TRUE : T_FALSE;
}

void
nmod8_init(nmod8_t x, const gr_ctx_t ctx)
{
    x[0] = 0;
}

void
nmod8_clear(nmod8_t x, const gr_ctx_t ctx)
{
}

void
nmod8_swap(nmod8_t x, nmod8_t y, const gr_ctx_t ctx)
{
    nmod8_t t;
    *t = *x;
    *x = *y;
    *y = *t;
}

void
nmod8_set_shallow(nmod8_t res, const nmod8_t x, const gr_ctx_t ctx)
{
    *res = *x;
}

int
nmod8_randtest(nmod8_t res, flint_rand_t state, const gr_ctx_t ctx)
{
    res[0] = n_randtest(state) % NMOD8_CTX(ctx).n;
    return GR_SUCCESS;
}

int
nmod8_write(gr_stream_t out, const nmod8_t x, const gr_ctx_t ctx)
{
    gr_stream_write_ui(out, x[0]);
    return GR_SUCCESS;
}

int
nmod8_zero(nmod8_t x, const gr_ctx_t ctx)
{
    x[0] = 0;
    return GR_SUCCESS;
}

int
nmod8_one(nmod8_t x, const gr_ctx_t ctx)
{
    x[0] = (NMOD8_CTX(ctx).n != 1);
    return GR_SUCCESS;
}

int
nmod8_set_si(nmod8_t res, slong v, const gr_ctx_t ctx)
{
    ulong t;
    nmod_t mod = NMOD8_CTX(ctx);
    t = FLINT_ABS(v);
    NMOD_RED(t, t, mod);
    if (v < 0)
        t = nmod_neg(t, mod);
    res[0] = t;
    return GR_SUCCESS;
}

int
nmod8_set_ui(nmod8_t res, ulong v, const gr_ctx_t ctx)
{
    ulong t;
    nmod_t mod = NMOD8_CTX(ctx);
    NMOD_RED(t, v, mod);
    res[0] = t;
    return GR_SUCCESS;
}

int
nmod8_set_fmpz(nmod8_t res, const fmpz_t v, const gr_ctx_t ctx)
{
    nmod_t mod = NMOD8_CTX(ctx);
    res[0] = fmpz_get_nmod(v, mod);
    return GR_SUCCESS;
}

truth_t
nmod8_is_zero(const nmod8_t x, const gr_ctx_t ctx)
{
    return (x[0] == 0) ? T_TRUE : T_FALSE;
}

truth_t
nmod8_is_one(const nmod8_t x, const gr_ctx_t ctx)
{
    return (x[0] == (NMOD8_CTX(ctx).n != 1)) ? T_TRUE : T_FALSE;
}

truth_t
nmod8_is_neg_one(const nmod8_t x, const gr_ctx_t ctx)
{
    return (x[0] == NMOD8_CTX(ctx).n - 1) ? T_TRUE : T_FALSE;
}

truth_t
nmod8_equal(const nmod8_t x, const nmod8_t y, const gr_ctx_t ctx)
{
    return (x[0] == y[0]) ? T_TRUE : T_FALSE;
}

int
nmod8_set(nmod8_t res, const nmod8_t x, const gr_ctx_t ctx)
{
    res[0] = x[0];
    return GR_SUCCESS;
}

int
nmod8_neg(nmod8_t res, const nmod8_t x, const gr_ctx_t ctx)
{
    res[0] = nmod_neg(x[0], NMOD8_CTX(ctx));
    return GR_SUCCESS;
}

int
nmod8_add(nmod8_t res, const nmod8_t x, const nmod8_t y, const gr_ctx_t ctx)
{
    res[0] = nmod_add(x[0], y[0], NMOD8_CTX(ctx));
    return GR_SUCCESS;
}

int
nmod8_add_si(nmod8_t res, const nmod8_t x, slong y, const gr_ctx_t ctx)
{
    nmod8_t t;
    nmod8_set_si(t, y, ctx);
    return nmod8_add(res, x, t, ctx);
}

int
nmod8_sub(nmod8_t res, const nmod8_t x, const nmod8_t y, const gr_ctx_t ctx)
{
    res[0] = nmod_sub(x[0], y[0], NMOD8_CTX(ctx));
    return GR_SUCCESS;
}

int
nmod8_mul(nmod8_t res, const nmod8_t x, const nmod8_t y, const gr_ctx_t ctx)
{
    res[0] = nmod_mul(x[0], y[0], NMOD8_CTX(ctx));
    return GR_SUCCESS;
}

int
nmod8_mul_si(nmod8_t res, const nmod8_t x, slong y, const gr_ctx_t ctx)
{
    nmod8_t t;
    nmod8_set_si(t, y, ctx);
    return nmod8_mul(res, x, t, ctx);
}

int
nmod8_addmul(nmod8_t res, const nmod8_t x, const nmod8_t y, const gr_ctx_t ctx)
{
    ulong r = res[0];
    NMOD_ADDMUL(r, x[0], y[0], NMOD8_CTX(ctx));
    res[0] = r;
    return GR_SUCCESS;
}

int
nmod8_submul(nmod8_t res, const nmod8_t x, const nmod8_t y, const gr_ctx_t ctx)
{
    ulong r = res[0];
    ulong t = nmod_neg(y[0], NMOD8_CTX(ctx));
    NMOD_ADDMUL(r, x[0], t, NMOD8_CTX(ctx));
    res[0] = r;
    return GR_SUCCESS;
}

int
nmod8_mul_two(nmod8_t res, const nmod8_t x, const gr_ctx_t ctx)
{
    return nmod8_add(res, x, x, ctx);
}

int
nmod8_sqr(nmod8_t res, const nmod8_t x, const gr_ctx_t ctx)
{
    return nmod8_mul(res, x, x, ctx);
}

int
nmod8_inv(nmod8_t res, const nmod8_t x, const gr_ctx_t ctx)
{
    ulong r, g;

    /* todo: also handle -1 fast? */
    if (x[0] == 1)
    {
        res[0] = x[0];
        return GR_SUCCESS;
    }

    g = n_gcdinv(&r, x[0], NMOD8_CTX(ctx).n);
    if (g == 1)
    {
        res[0] = r;
        return GR_SUCCESS;
    }
    else
    {
        res[0] = 0;
        return GR_DOMAIN;
    }
}

int
nmod8_div(nmod8_t res, const nmod8_t x, const nmod8_t y, const gr_ctx_t ctx)
{
    nmod8_t t;
    int status;

    status = nmod8_inv(t, y, ctx);

    if (status == GR_SUCCESS)
        nmod8_mul(res, x, t, ctx);

    return status;
}

int
nmod8_div_si(nmod8_t res, const nmod8_t x, slong y, const gr_ctx_t ctx)
{
    nmod8_t t;
    nmod8_set_si(t, y, ctx);
    return nmod8_div(res, x, t, ctx);
}

int
nmod8_div_ui(nmod8_t res, const nmod8_t x, ulong y, const gr_ctx_t ctx)
{
    nmod8_t t;
    nmod8_set_ui(t, y, ctx);
    return nmod8_div(res, x, t, ctx);
}

int
nmod8_div_fmpz(nmod8_t res, const nmod8_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    nmod8_t t;
    nmod8_set_fmpz(t, y, ctx);
    return nmod8_div(res, x, t, ctx);
}

truth_t
nmod8_is_invertible(const nmod8_t x, const gr_ctx_t ctx)
{
    ulong r, g;
    g = n_gcdinv(&r, x[0], NMOD8_CTX(ctx).n);
    return (g == 1) ? T_TRUE : T_FALSE;
}

truth_t
nmod8_divides(const nmod8_t x, const nmod8_t y, const gr_ctx_t ctx)
{
    ulong t;
    return nmod_divides(&t, y[0], x[0], NMOD8_CTX(ctx)) ? T_TRUE : T_FALSE;
}

int
nmod8_div_nonunique(nmod8_t res, const nmod8_t x, const nmod8_t y, const gr_ctx_t ctx)
{
    nmod8_t t;
    int status;

    status = nmod8_inv(t, y, ctx);

    if (status == GR_SUCCESS)
    {
        nmod8_mul(res, x, t, ctx);
    }
    else
    {
        ulong t2;
        status = nmod_divides(&t2, x[0], y[0], NMOD8_CTX(ctx)) ? GR_SUCCESS : GR_DOMAIN;
        res[0] = t2;
    }

    return status;
}

void
_nmod8_vec_init(nmod8_struct * res, slong len, gr_ctx_t ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        res[i] = 0;
}

void
_nmod8_vec_clear(nmod8_struct * res, slong len, gr_ctx_t ctx)
{
}

int
_nmod8_vec_set(nmod8_struct * res, const nmod8_struct * vec, slong len, gr_ctx_t ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        res[i] = vec[i];

    return GR_SUCCESS;
}

int
_nmod8_vec_neg(nmod8_struct * res, const nmod8_struct * vec, slong len, gr_ctx_t ctx)
{
    slong i;
    nmod_t mod = NMOD8_CTX(ctx);

    for (i = 0; i < len; i++)
        res[i] = nmod_neg(vec[i], mod);

    return GR_SUCCESS;
}

int
_nmod8_vec_add(nmod8_struct * res, const nmod8_struct * vec1, const nmod8_struct * vec2, slong len, gr_ctx_t ctx)
{
    slong i;
    nmod_t mod = NMOD8_CTX(ctx);

    for (i = 0; i < len; i++)
        res[i] = _nmod_add(vec1[i], vec2[i], mod);

    return GR_SUCCESS;
}


int
_nmod8_vec_sub(nmod8_struct * res, const nmod8_struct * vec1, const nmod8_struct * vec2, slong len, gr_ctx_t ctx)
{
    slong i;
    nmod_t mod = NMOD8_CTX(ctx);

    for (i = 0; i < len; i++)
        res[i] = nmod_sub(vec1[i], vec2[i], mod);

    return GR_SUCCESS;
}

int
_nmod8_vec_dot(nmod8_t res, const nmod8_t initial, int subtract, const nmod8_struct * vec1, const nmod8_struct * vec2, slong len, gr_ctx_t ctx)
{
    slong i;
    ulong n, s;

    if (len <= 0)
    {
        if (initial == NULL)
            nmod8_zero(res, ctx);
        else
            nmod8_set(res, initial, ctx);
        return GR_SUCCESS;
    }

    n = NMOD8_CTX(ctx).n;

    if (initial == NULL)
    {
        s = 0;
    }
    else
    {
        if (subtract)
            s = n_negmod(initial[0], n);
        else
            s = initial[0];
    }

    if (len < 65536)
    {
        unsigned int ss = 0;

        for (i = 0; i + 4 < len; i += 4)
        {
            s += (unsigned int) vec1[i + 0] * (unsigned int) vec2[i + 0];
            s += (unsigned int) vec1[i + 1] * (unsigned int) vec2[i + 1];
            s += (unsigned int) vec1[i + 2] * (unsigned int) vec2[i + 2];
            s += (unsigned int) vec1[i + 3] * (unsigned int) vec2[i + 3];
        }

        for ( ; i < len; i++)
            s += (unsigned int) vec1[i] * (unsigned int) vec2[i];

        s += ss;
    }
    else
    {
        ulong ss;

        const dot_params_t params = _nmod_vec_dot_params(len, NMOD8_CTX(ctx));
        NMOD_VEC_DOT(ss, i, len, (ulong) vec1[i], (ulong) vec2[i], NMOD8_CTX(ctx), params);
        s = n_addmod(s, ss, n);
    }

    nmod8_set_ui(res, s, ctx);

    if (subtract && res[0] != 0)
        res[0] = (n - res[0]);

    return GR_SUCCESS;
}

int
_nmod8_vec_dot_rev(nmod8_t res, const nmod8_t initial, int subtract, const nmod8_struct * vec1, const nmod8_struct * vec2, slong len, gr_ctx_t ctx)
{
    slong i;
    ulong n, s;

    if (len <= 0)
    {
        if (initial == NULL)
            nmod8_zero(res, ctx);
        else
            nmod8_set(res, initial, ctx);
        return GR_SUCCESS;
    }

    n = NMOD8_CTX(ctx).n;

    if (initial == NULL)
    {
        s = 0;
    }
    else
    {
        if (subtract)
            s = n_negmod(initial[0], n);
        else
            s = initial[0];
    }

    if (len < 65536)
    {
        unsigned int ss = 0;

        for (i = 0; i + 4 < len; i += 4)
        {
            s += (unsigned int) vec1[i + 0] * (unsigned int) vec2[len - 1 - i - 0];
            s += (unsigned int) vec1[i + 1] * (unsigned int) vec2[len - 1 - i - 1];
            s += (unsigned int) vec1[i + 2] * (unsigned int) vec2[len - 1 - i - 2];
            s += (unsigned int) vec1[i + 3] * (unsigned int) vec2[len - 1 - i - 3];
        }

        for ( ; i < len; i++)
            s += (unsigned int) vec1[i] * (unsigned int) vec2[len - 1 - i];

        s += ss;
    }
    else
    {
        ulong ss;

        const dot_params_t params = _nmod_vec_dot_params(len, NMOD8_CTX(ctx));
        NMOD_VEC_DOT(ss, i, len, (ulong) vec1[i], (ulong) vec2[len - 1 - i], NMOD8_CTX(ctx), params);
        s = n_addmod(s, ss, n);
    }

    nmod8_set_ui(res, s, ctx);

    if (subtract && res[0] != 0)
        res[0] = (n - res[0]);

    return GR_SUCCESS;
}

/* todo: tuning for rectangular matrices */
int
_nmod8_mat_mul(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    if (A->r >= 256 && A->c >= 256 && B->c >= 256)
        return gr_mat_mul_strassen(C, A, B, ctx);
    else
        return gr_mat_mul_classical(C, A, B, ctx);
}


int _nmod8_methods_initialized = 0;

gr_static_method_table _nmod8_methods;

gr_method_tab_input _nmod8_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) nmod8_ctx_write},
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) nmod8_ctx_is_field},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) nmod8_ctx_is_field},
    {GR_METHOD_CTX_IS_FINITE,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_EXACT,    (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_CANONICAL,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_INIT,            (gr_funcptr) nmod8_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) nmod8_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) nmod8_swap},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) nmod8_set_shallow},
    {GR_METHOD_RANDTEST,        (gr_funcptr) nmod8_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) nmod8_write},
    {GR_METHOD_ZERO,            (gr_funcptr) nmod8_zero},
    {GR_METHOD_ONE,             (gr_funcptr) nmod8_one},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) nmod8_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) nmod8_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) nmod8_is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) nmod8_equal},
    {GR_METHOD_SET,             (gr_funcptr) nmod8_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) nmod8_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) nmod8_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) nmod8_set_fmpz},
    {GR_METHOD_NEG,             (gr_funcptr) nmod8_neg},
    {GR_METHOD_ADD,             (gr_funcptr) nmod8_add},
    {GR_METHOD_ADD_SI,          (gr_funcptr) nmod8_add_si},
    {GR_METHOD_SUB,             (gr_funcptr) nmod8_sub},
    {GR_METHOD_MUL,             (gr_funcptr) nmod8_mul},
    {GR_METHOD_MUL_SI,          (gr_funcptr) nmod8_mul_si},
    {GR_METHOD_ADDMUL,          (gr_funcptr) nmod8_addmul},
    {GR_METHOD_SUBMUL,          (gr_funcptr) nmod8_submul},
    {GR_METHOD_MUL_TWO,         (gr_funcptr) nmod8_mul_two},
    {GR_METHOD_SQR,             (gr_funcptr) nmod8_sqr},
    {GR_METHOD_DIV,             (gr_funcptr) nmod8_div},
    {GR_METHOD_DIV_SI,          (gr_funcptr) nmod8_div_si},
    {GR_METHOD_DIV_UI,          (gr_funcptr) nmod8_div_ui},
    {GR_METHOD_DIV_FMPZ,        (gr_funcptr) nmod8_div_fmpz},
    {GR_METHOD_DIV_NONUNIQUE,   (gr_funcptr) nmod8_div_nonunique},
    {GR_METHOD_DIVIDES,         (gr_funcptr) nmod8_divides},
    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) nmod8_is_invertible},
    {GR_METHOD_INV,             (gr_funcptr) nmod8_inv},
    {GR_METHOD_VEC_INIT,        (gr_funcptr) _nmod8_vec_init},
    {GR_METHOD_VEC_CLEAR,       (gr_funcptr) _nmod8_vec_clear},
    {GR_METHOD_VEC_SET,         (gr_funcptr) _nmod8_vec_set},
    {GR_METHOD_VEC_NEG,         (gr_funcptr) _nmod8_vec_neg},
    {GR_METHOD_VEC_ADD,         (gr_funcptr) _nmod8_vec_add},
    {GR_METHOD_VEC_SUB,         (gr_funcptr) _nmod8_vec_sub},
    {GR_METHOD_VEC_DOT,         (gr_funcptr) _nmod8_vec_dot},
    {GR_METHOD_VEC_DOT_REV,     (gr_funcptr) _nmod8_vec_dot_rev},
    {GR_METHOD_MAT_MUL,         (gr_funcptr) _nmod8_mat_mul},
    {0,                         (gr_funcptr) NULL},
};

void
gr_ctx_init_nmod8(gr_ctx_t ctx, unsigned char n)
{
    ctx->which_ring = GR_CTX_NMOD8;
    ctx->sizeof_elem = sizeof(nmod8_struct);
    ctx->size_limit = WORD_MAX;

    nmod_init(NMOD8_CTX_REF(ctx), n);

    ctx->methods = _nmod8_methods;

    if (!_nmod8_methods_initialized)
    {
        gr_method_tab_init(_nmod8_methods, _nmod8_methods_input);
        _nmod8_methods_initialized = 1;
    }
}
