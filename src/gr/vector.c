/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Vectors over generic rings */

#include "ulong_extras.h"
#include "fmpz.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "gr_special.h"

#define ENTRY_CTX(ctx) (VECTOR_CTX(ctx)->base_ring)

int
_gr_vec_check_resize(gr_vec_t res, slong n, gr_ctx_t ctx)
{
    if (VECTOR_CTX(ctx)->all_sizes)
    {
        gr_vec_set_length(res, n, ENTRY_CTX(ctx));
        return GR_SUCCESS;
    }
    else
    {
        if (n != VECTOR_CTX(ctx)->n)
            return GR_DOMAIN;

        gr_vec_set_length(res, n, ENTRY_CTX(ctx));
        return GR_SUCCESS;
    }
}

void
vector_gr_vec_init(gr_vec_t res, gr_ctx_t ctx)
{
    gr_vec_init(res, VECTOR_CTX(ctx)->n, ENTRY_CTX(ctx));
}

int vector_gr_vec_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_ctx_ptr elem_ctx = ENTRY_CTX(ctx);

    if (VECTOR_CTX(ctx)->all_sizes)
    {
        gr_stream_write(out, "Vectors (any length) over ");
    }
    else
    {
        gr_stream_write(out, "Space of length ");
        gr_stream_write_si(out, VECTOR_CTX(ctx)->n);
        gr_stream_write(out, " vectors over ");
    }

    gr_ctx_write(out, elem_ctx);
    return GR_SUCCESS;
}

/* todo: public */
truth_t gr_ctx_vector_gr_vec_is_fixed_size(gr_ctx_t ctx)
{
    return (VECTOR_CTX(ctx)->all_sizes) ? T_FALSE : T_TRUE;
}

truth_t
vector_ctx_is_threadsafe(gr_ctx_t ctx)
{
    return gr_ctx_is_threadsafe(ENTRY_CTX(ctx));
}

void
vector_gr_vec_clear(gr_vec_t res, gr_ctx_t ctx)
{
    gr_vec_clear(res, ENTRY_CTX(ctx));
}

void
vector_gr_vec_swap(gr_vec_t vec1, gr_vec_t vec2, gr_ctx_t ctx)
{
    gr_poly_swap((gr_poly_struct *) vec1, (gr_poly_struct *) vec2, ENTRY_CTX(ctx));
}

int
vector_gr_vec_write(gr_stream_t out, gr_vec_t vec, gr_ctx_t ctx)
{
    return gr_vec_write(out, vec, ENTRY_CTX(ctx));
}

int
vector_gr_vec_randtest(gr_vec_t res, flint_rand_t state, gr_ctx_t ctx)
{
    slong i, n;
    int status;

    if (VECTOR_CTX(ctx)->all_sizes)
    {
        n = n_randint(state, 7);
        gr_vec_set_length(res, n, ENTRY_CTX(ctx));
    }
    else
        n = res->length;


    status = GR_SUCCESS;
    for (i = 0; i < n; i++)
        status |= gr_randtest(gr_vec_entry_ptr(res, i, ENTRY_CTX(ctx)), state, ENTRY_CTX(ctx));

    return status;
}

truth_t
gr_generic_vec_equal(gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx);

truth_t
vector_gr_vec_equal(const gr_vec_t vec1, const gr_vec_t vec2, gr_ctx_t ctx)
{
    slong len1, len2;

    len1 = vec1->length;
    len2 = vec2->length;

    if (len1 != len2)
        return T_FALSE;

    return _gr_vec_equal(vec1->entries, vec2->entries, len1, ENTRY_CTX(ctx));
}

int
vector_gr_vec_set(gr_vec_t res, const gr_vec_t vec, gr_ctx_t ctx)
{
    return gr_vec_set(res, vec, ENTRY_CTX(ctx));
}

/* todo: move out */
static int
_gr_vec_set_ui(gr_ptr res, slong len, ulong x, gr_ctx_t ctx)
{
    gr_method_unary_op_ui f = GR_UNARY_OP_UI(ctx, SET_UI);
    int status;
    slong i, sz;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    for (i = 0; i < len; i++)
        status |= f(GR_ENTRY(res, i, sz), x, ctx);

    return status;
}

static int
_gr_vec_set_si(gr_ptr res, slong len, slong x, gr_ctx_t ctx)
{
    gr_method_unary_op_si f = GR_UNARY_OP_SI(ctx, SET_SI);
    int status;
    slong i, sz;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    for (i = 0; i < len; i++)
        status |= f(GR_ENTRY(res, i, sz), x, ctx);

    return status;
}

static int
_gr_vec_set_fmpz(gr_ptr res, slong len, const fmpz_t x, gr_ctx_t ctx)
{
    gr_method_unary_op_fmpz f = GR_UNARY_OP_FMPZ(ctx, SET_FMPZ);
    int status;
    slong i, sz;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    for (i = 0; i < len; i++)
        status |= f(GR_ENTRY(res, i, sz), x, ctx);

    return status;
}

static int
_gr_vec_set_fmpq(gr_ptr res, slong len, const fmpq_t x, gr_ctx_t ctx)
{
    gr_method_unary_op_fmpq f = GR_UNARY_OP_FMPQ(ctx, SET_FMPQ);
    int status;
    slong i, sz;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    for (i = 0; i < len; i++)
        status |= f(GR_ENTRY(res, i, sz), x, ctx);

    return status;
}


static int
_gr_vec_set_d(gr_ptr res, slong len, double x, gr_ctx_t ctx)
{
    gr_method_unary_op_d f = GR_UNARY_OP_D(ctx, SET_D);
    int status;
    slong i, sz;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    for (i = 0; i < len; i++)
        status |= f(GR_ENTRY(res, i, sz), x, ctx);

    return status;
}

int
vector_gr_vec_set_ui(gr_vec_t res, ulong x, gr_ctx_t ctx)
{
    slong len = VECTOR_CTX(ctx)->n;

    if (VECTOR_CTX(ctx)->all_sizes)
        return GR_DOMAIN;

    /* should not be needed; just a precaution (might want an assert) */
    if (res->length != len)
        gr_vec_set_length(res, len, ENTRY_CTX(ctx));

    return _gr_vec_set_ui(res->entries, len, x, ENTRY_CTX(ctx));
}

int
vector_gr_vec_set_si(gr_vec_t res, slong x, gr_ctx_t ctx)
{
    slong len = VECTOR_CTX(ctx)->n;

    if (VECTOR_CTX(ctx)->all_sizes)
        return GR_DOMAIN;

    /* should not be needed; just a precaution (might want an assert) */
    if (res->length != len)
        gr_vec_set_length(res, len, ENTRY_CTX(ctx));

    return _gr_vec_set_si(res->entries, len, x, ENTRY_CTX(ctx));
}

int
vector_gr_vec_set_fmpz(gr_vec_t res, const fmpz_t x, gr_ctx_t ctx)
{
    slong len = VECTOR_CTX(ctx)->n;

    if (VECTOR_CTX(ctx)->all_sizes)
        return GR_DOMAIN;

    /* should not be needed; just a precaution (might want an assert) */
    if (res->length != len)
        gr_vec_set_length(res, len, ENTRY_CTX(ctx));

    return _gr_vec_set_fmpz(res->entries, len, x, ENTRY_CTX(ctx));
}

int
vector_gr_vec_set_fmpq(gr_vec_t res, const fmpq_t x, gr_ctx_t ctx)
{
    slong len = VECTOR_CTX(ctx)->n;

    if (VECTOR_CTX(ctx)->all_sizes)
        return GR_DOMAIN;

    /* should not be needed; just a precaution (might want an assert) */
    if (res->length != len)
        gr_vec_set_length(res, len, ENTRY_CTX(ctx));

    return _gr_vec_set_fmpq(res->entries, len, x, ENTRY_CTX(ctx));
}

int
vector_gr_vec_set_d(gr_vec_t res, double x, gr_ctx_t ctx)
{
    slong len = VECTOR_CTX(ctx)->n;

    if (VECTOR_CTX(ctx)->all_sizes)
        return GR_DOMAIN;

    /* should not be needed; just a precaution (might want an assert) */
    if (res->length != len)
        gr_vec_set_length(res, len, ENTRY_CTX(ctx));

    return _gr_vec_set_d(res->entries, len, x, ENTRY_CTX(ctx));
}


/* todo: convert from matrices, ...? */
int
vector_gr_vec_set_other(gr_vec_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    if (x_ctx == ctx)
    {
        return vector_gr_vec_set(res, x, ctx);
    }
    else if (x_ctx->which_ring == GR_CTX_GR_VEC)
    {
        const gr_vec_struct * xvec = x;
        slong i;
        int status = GR_SUCCESS;
        slong sz, xsz;

        if (res->length != xvec->length)
        {
            if (VECTOR_CTX(ctx)->all_sizes)
                gr_vec_set_length(res, xvec->length, ENTRY_CTX(ctx));
            else
                return GR_DOMAIN;
        }

        sz = ENTRY_CTX(ctx)->sizeof_elem;
        xsz = ENTRY_CTX(x_ctx)->sizeof_elem;

        for (i = 0; i < res->length; i++)
        {
            status = gr_set_other(GR_ENTRY(res->entries, i, sz),
                        GR_ENTRY(xvec->entries, i, xsz),
                        ENTRY_CTX(x_ctx),
                        ENTRY_CTX(ctx));

            if (status != GR_SUCCESS)
                return status;
        }

        return GR_SUCCESS;
    }

    return GR_UNABLE;
}

int
vector_gr_vec_neg(gr_vec_t res, const gr_vec_t src, gr_ctx_t ctx)
{
    return gr_poly_neg((gr_poly_struct *) res, (gr_poly_struct *) src, ENTRY_CTX(ctx));
}

int
vector_gr_vec_zero(gr_vec_t res, gr_ctx_t ctx) \
{
    slong xlen = VECTOR_CTX(ctx)->n; \

    if (VECTOR_CTX(ctx)->all_sizes) \
        return GR_DOMAIN;

    /* should not be needed; just a precaution (might want an assert) */
    if (res->length != xlen)
        gr_vec_set_length(res, xlen, ENTRY_CTX(ctx));

    return _gr_vec_zero(res->entries, xlen, ENTRY_CTX(ctx));
}

static int
_gr_vec_apply_const(gr_ptr res, gr_method_constant_op f, slong len, gr_ctx_t ctx)
{
    int status;
    slong i, sz;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    for (i = 0; i < len; i++)
        status |= f(GR_ENTRY(res, i, sz), ctx);

    return status;
}

static truth_t
_gr_vec_all_binary_predicate(gr_method_binary_predicate f, gr_srcptr x, gr_srcptr y, slong len, gr_ctx_t ctx)
{
    truth_t res, res1;
    slong i, sz;

    sz = ctx->sizeof_elem;
    res = T_TRUE;

    for (i = 0; i < len; i++)
    {
        res1 = f(GR_ENTRY(x, i, sz), GR_ENTRY(y, i, sz), ctx);

        if (res1 == T_FALSE)
            return T_FALSE;

        if (res1 == T_UNKNOWN)
            res = T_UNKNOWN;
    }

    return res;
}

truth_t
vector_gr_vec_divides(const gr_vec_t x, const gr_vec_t y, gr_ctx_t ctx)
{
    if (x->length != y->length)
        return T_FALSE;

    return _gr_vec_all_binary_predicate(GR_BINARY_PREDICATE(ENTRY_CTX(ctx), DIVIDES), x->entries, y->entries, x->length, ENTRY_CTX(ctx));
}

#define DEF_CONSTANT_OP_FROM_OP(op, OP) \
int \
vector_gr_vec_ ## op(gr_vec_t res, gr_ctx_t ctx) \
{ \
    slong xlen = VECTOR_CTX(ctx)->n; \
 \
    if (VECTOR_CTX(ctx)->all_sizes) \
        return GR_DOMAIN;\
 \
    gr_method_constant_op f = GR_CONSTANT_OP(ENTRY_CTX(ctx), OP); \
 \
    if (res->length != xlen) \
        gr_vec_set_length(res, xlen, ENTRY_CTX(ctx)); \
 \
    return _gr_vec_apply_const(res->entries, f, xlen, ENTRY_CTX(ctx)); \
} \

DEF_CONSTANT_OP_FROM_OP(one, ONE)
DEF_CONSTANT_OP_FROM_OP(neg_one, NEG_ONE)
DEF_CONSTANT_OP_FROM_OP(i, I)
DEF_CONSTANT_OP_FROM_OP(pi, PI)


static int
_gr_vec_apply_unary(gr_ptr res, gr_method_unary_op f, gr_srcptr src, slong len, gr_ctx_t ctx)
{
    int status;
    slong i, sz;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    for (i = 0; i < len; i++)
        status |= f(GR_ENTRY(res, i, sz), GR_ENTRY(src, i, sz), ctx);

    return status;
}

#define DEF_UNARY_OP_FROM_ENTRY_OP(op, OP) \
int \
vector_gr_vec_ ## op(gr_vec_t res, const gr_vec_t x, gr_ctx_t ctx) \
{ \
    slong xlen = x->length; \
    gr_method_unary_op f = GR_UNARY_OP(ENTRY_CTX(ctx), OP); \
 \
    if (res->length != xlen) \
        gr_vec_set_length(res, xlen, ENTRY_CTX(ctx)); \
 \
    return _gr_vec_apply_unary(res->entries, f, x->entries, xlen, ENTRY_CTX(ctx)); \
} \

DEF_UNARY_OP_FROM_ENTRY_OP(inv, INV)
DEF_UNARY_OP_FROM_ENTRY_OP(sqrt, SQRT)
DEF_UNARY_OP_FROM_ENTRY_OP(rsqrt, RSQRT)
DEF_UNARY_OP_FROM_ENTRY_OP(floor, FLOOR)
DEF_UNARY_OP_FROM_ENTRY_OP(ceil, CEIL)
DEF_UNARY_OP_FROM_ENTRY_OP(trunc, TRUNC)
DEF_UNARY_OP_FROM_ENTRY_OP(nint, NINT)
DEF_UNARY_OP_FROM_ENTRY_OP(abs, ABS)
DEF_UNARY_OP_FROM_ENTRY_OP(conj, CONJ)
DEF_UNARY_OP_FROM_ENTRY_OP(re, RE)
DEF_UNARY_OP_FROM_ENTRY_OP(im, IM)
DEF_UNARY_OP_FROM_ENTRY_OP(sgn, SGN)
DEF_UNARY_OP_FROM_ENTRY_OP(csgn, CSGN)

DEF_UNARY_OP_FROM_ENTRY_OP(exp, EXP)
DEF_UNARY_OP_FROM_ENTRY_OP(log, LOG)


#define DEF_BINARY_OP(op) \
int \
vector_gr_vec_ ## op(gr_vec_t res, const gr_vec_t x, const gr_vec_t y, gr_ctx_t ctx) \
{ \
    slong xlen = x->length; \
 \
    if (xlen != y->length) \
        return GR_DOMAIN; \
 \
    if (res->length != xlen) \
        gr_vec_set_length(res, xlen, ENTRY_CTX(ctx)); \
 \
    return _gr_vec_## op(res->entries, x->entries, y->entries, xlen, ENTRY_CTX(ctx)); \
} \
 \
int \
vector_gr_vec_ ## op ## _si(gr_vec_t res, const gr_vec_t x, slong c, gr_ctx_t ctx) \
{ \
    slong xlen = x->length; \
 \
    if (res->length != xlen) \
        gr_vec_set_length(res, xlen, ENTRY_CTX(ctx)); \
 \
    return _gr_vec_ ## op ## _scalar_si(res->entries, x->entries, xlen, c, ENTRY_CTX(ctx)); \
} \
 \
int \
vector_gr_vec_ ## op ## _ui(gr_vec_t res, const gr_vec_t x, ulong c, gr_ctx_t ctx) \
{ \
    slong xlen = x->length; \
 \
    if (res->length != xlen) \
        gr_vec_set_length(res, xlen, ENTRY_CTX(ctx)); \
 \
    return _gr_vec_ ## op ## _scalar_ui(res->entries, x->entries, xlen, c, ENTRY_CTX(ctx)); \
} \
 \
int \
vector_gr_vec_ ## op ## _fmpz(gr_vec_t res, const gr_vec_t x, const fmpz_t c, gr_ctx_t ctx) \
{ \
    slong xlen = x->length; \
 \
    if (res->length != xlen) \
        gr_vec_set_length(res, xlen, ENTRY_CTX(ctx)); \
 \
    return _gr_vec_ ## op ## _scalar_fmpz(res->entries, x->entries, xlen, c, ENTRY_CTX(ctx)); \
} \
 \
int \
vector_gr_vec_ ## op ## _fmpq(gr_vec_t res, const gr_vec_t x, const fmpq_t c, gr_ctx_t ctx) \
{ \
    slong xlen = x->length; \
 \
    if (res->length != xlen) \
        gr_vec_set_length(res, xlen, ENTRY_CTX(ctx)); \
 \
    return _gr_vec_ ## op ## _scalar_fmpq(res->entries, x->entries, xlen, c, ENTRY_CTX(ctx)); \
} \
 \
int \
vector_gr_vec_ ## op ## _other(gr_vec_t res, const gr_vec_t x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx) \
{ \
    slong xlen = x->length; \
 \
    if (y_ctx == ctx) \
    { \
        return vector_gr_vec_ ## op(res, x, y, ctx); \
    } \
    else if (y_ctx == ENTRY_CTX(ctx)) \
    { \
        if (res->length != xlen) \
            gr_vec_set_length(res, xlen, y_ctx); \
 \
        return _gr_vec_ ## op ## _scalar(res->entries, x->entries, xlen, y, y_ctx); \
    } \
    else if (y_ctx->which_ring == GR_CTX_GR_VEC) \
    { \
        const gr_vec_struct * yvec = y; \
        gr_ctx_struct * entry_ctx = ENTRY_CTX(ctx); \
        gr_ctx_struct * y_entry_ctx = ENTRY_CTX(y_ctx); \
 \
        if (xlen != yvec->length) \
            return GR_DOMAIN; \
 \
        if (res->length != xlen) \
            gr_vec_set_length(res, xlen, entry_ctx); \
 \
        return _gr_vec_ ## op ## _other(res->entries, x->entries, yvec->entries, y_entry_ctx, xlen, entry_ctx); \
    } \
    else \
    { \
        gr_ctx_struct * entry_ctx = ENTRY_CTX(ctx); \
 \
        if (res->length != xlen) \
            gr_vec_set_length(res, xlen, entry_ctx); \
 \
        return _gr_vec_ ## op ## _scalar_other(res->entries, x->entries, xlen, y, y_ctx, entry_ctx); \
    } \
} \
int \
vector_gr_vec_other_ ## op(gr_vec_t res, gr_srcptr x, gr_ctx_t x_ctx, const gr_vec_t y, gr_ctx_t ctx) \
{ \
    slong ylen = y->length; \
 \
    if (x_ctx == ctx) \
    { \
        return vector_gr_vec_ ## op(res, x, y, ctx); \
    } \
    else if (x_ctx == ENTRY_CTX(ctx)) \
    { \
        if (res->length != ylen) \
            gr_vec_set_length(res, ylen, x_ctx); \
 \
        return _gr_scalar_ ## op ## _vec(res->entries, x, y->entries, ylen, x_ctx); \
    } \
    else if (x_ctx->which_ring == GR_CTX_GR_VEC) \
    { \
        const gr_vec_struct * xvec = x; \
        gr_ctx_struct * entry_ctx = ENTRY_CTX(ctx); \
        gr_ctx_struct * x_entry_ctx = ENTRY_CTX(x_ctx); \
 \
        if (ylen != xvec->length) \
            return GR_DOMAIN; \
 \
        if (res->length != ylen) \
            gr_vec_set_length(res, ylen, entry_ctx); \
 \
        return _gr_other_ ## op ## _vec(res->entries, xvec->entries, x_entry_ctx, y->entries, ylen, entry_ctx); \
    } \
    else \
    { \
        gr_ctx_struct * entry_ctx = ENTRY_CTX(ctx); \
 \
        if (res->length != ylen) \
            gr_vec_set_length(res, ylen, entry_ctx); \
 \
        return _gr_scalar_other_ ## op ## _vec(res->entries, x, x_ctx, y->entries, ylen, entry_ctx); \
    } \
} \

DEF_BINARY_OP(add)
DEF_BINARY_OP(sub)
DEF_BINARY_OP(mul)
DEF_BINARY_OP(div)
DEF_BINARY_OP(divexact)
DEF_BINARY_OP(pow)

#define DEF_BINARY_OP_NO_TYPE_VARIANTS(op) \
int \
vector_gr_vec_ ## op(gr_vec_t res, const gr_vec_t x, const gr_vec_t y, gr_ctx_t ctx) \
{ \
    slong xlen = x->length; \
 \
    if (xlen != y->length) \
        return GR_DOMAIN; \
 \
    if (res->length != xlen) \
        gr_vec_set_length(res, xlen, ENTRY_CTX(ctx)); \
 \
    return _gr_vec_## op(res->entries, x->entries, y->entries, xlen, ENTRY_CTX(ctx)); \
} \

int
_gr_vec_div_nonunique(gr_ptr res, gr_srcptr x, gr_srcptr y, slong len, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    slong i;

    for (i = 0; i < len; i++)
        status |= gr_div_nonunique(GR_ENTRY(res, i, sz), GR_ENTRY(x, i, sz), GR_ENTRY(y, i, sz), ctx);

    return status;
}

DEF_BINARY_OP_NO_TYPE_VARIANTS(div_nonunique)


/* todo: all versions */

int gr_generic_mul_ui_via_ZZ(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
{
    gr_ctx_t ZZ;
    fmpz_t t;
    int status;

    gr_ctx_init_fmpz(ZZ); /* no need to clear */

    fmpz_init_set_ui(t, y);
    status = gr_mul_other(res, x, t, ZZ, ctx);
    fmpz_clear(t);
    return status;
}


int _gr_vec_methods_initialized = 0;

gr_static_method_table _gr_vec_methods;

gr_method_tab_input _gr_vec_methods_input[] =
{
    {GR_METHOD_CTX_IS_THREADSAFE,    (gr_funcptr) vector_ctx_is_threadsafe},

    {GR_METHOD_CTX_WRITE,   (gr_funcptr) vector_gr_vec_ctx_write},
    {GR_METHOD_INIT,        (gr_funcptr) vector_gr_vec_init},
    {GR_METHOD_CLEAR,       (gr_funcptr) vector_gr_vec_clear},
    {GR_METHOD_SWAP,        (gr_funcptr) vector_gr_vec_swap},
    {GR_METHOD_RANDTEST,    (gr_funcptr) vector_gr_vec_randtest},
    {GR_METHOD_WRITE,       (gr_funcptr) vector_gr_vec_write},
    {GR_METHOD_EQUAL,       (gr_funcptr) vector_gr_vec_equal},
    {GR_METHOD_SET,         (gr_funcptr) vector_gr_vec_set},
    {GR_METHOD_SET_OTHER,   (gr_funcptr) vector_gr_vec_set_other},
    {GR_METHOD_SET_UI,      (gr_funcptr) vector_gr_vec_set_ui},
    {GR_METHOD_SET_SI,      (gr_funcptr) vector_gr_vec_set_si},
    {GR_METHOD_SET_FMPZ,    (gr_funcptr) vector_gr_vec_set_fmpz},
    {GR_METHOD_SET_FMPQ,    (gr_funcptr) vector_gr_vec_set_fmpq},
    {GR_METHOD_SET_D,       (gr_funcptr) vector_gr_vec_set_d},

    {GR_METHOD_ZERO,        (gr_funcptr) vector_gr_vec_zero},
    {GR_METHOD_ONE,         (gr_funcptr) vector_gr_vec_one},
    {GR_METHOD_NEG_ONE,     (gr_funcptr) vector_gr_vec_neg_one},

    {GR_METHOD_NEG,         (gr_funcptr) vector_gr_vec_neg},

    {GR_METHOD_ADD,         (gr_funcptr) vector_gr_vec_add},
    {GR_METHOD_ADD_UI,      (gr_funcptr) vector_gr_vec_add_ui},
    {GR_METHOD_ADD_SI,      (gr_funcptr) vector_gr_vec_add_si},
    {GR_METHOD_ADD_FMPZ,    (gr_funcptr) vector_gr_vec_add_fmpz},
    {GR_METHOD_ADD_FMPQ,    (gr_funcptr) vector_gr_vec_add_fmpq},
    {GR_METHOD_ADD_OTHER,   (gr_funcptr) vector_gr_vec_add_other},
    {GR_METHOD_OTHER_ADD,   (gr_funcptr) vector_gr_vec_other_add},

    {GR_METHOD_SUB,         (gr_funcptr) vector_gr_vec_sub},
    {GR_METHOD_SUB_UI,      (gr_funcptr) vector_gr_vec_sub_ui},
    {GR_METHOD_SUB_SI,      (gr_funcptr) vector_gr_vec_sub_si},
    {GR_METHOD_SUB_FMPZ,    (gr_funcptr) vector_gr_vec_sub_fmpz},
    {GR_METHOD_SUB_FMPQ,    (gr_funcptr) vector_gr_vec_sub_fmpq},
    {GR_METHOD_SUB_OTHER,   (gr_funcptr) vector_gr_vec_sub_other},
    {GR_METHOD_OTHER_SUB,   (gr_funcptr) vector_gr_vec_other_sub},

    {GR_METHOD_MUL,         (gr_funcptr) vector_gr_vec_mul},
    {GR_METHOD_MUL_UI,      (gr_funcptr) vector_gr_vec_mul_ui},
    {GR_METHOD_MUL_SI,      (gr_funcptr) vector_gr_vec_mul_si},
    {GR_METHOD_MUL_FMPZ,    (gr_funcptr) vector_gr_vec_mul_fmpz},
    {GR_METHOD_MUL_FMPQ,    (gr_funcptr) vector_gr_vec_mul_fmpq},
    {GR_METHOD_MUL_OTHER,   (gr_funcptr) vector_gr_vec_mul_other},
    {GR_METHOD_OTHER_MUL,   (gr_funcptr) vector_gr_vec_other_mul},

    {GR_METHOD_INV,         (gr_funcptr) vector_gr_vec_inv},

    {GR_METHOD_DIV,         (gr_funcptr) vector_gr_vec_div},
    {GR_METHOD_DIV_UI,      (gr_funcptr) vector_gr_vec_div_ui},
    {GR_METHOD_DIV_SI,      (gr_funcptr) vector_gr_vec_div_si},
    {GR_METHOD_DIV_FMPZ,    (gr_funcptr) vector_gr_vec_div_fmpz},
    {GR_METHOD_DIV_FMPQ,    (gr_funcptr) vector_gr_vec_div_fmpq},
    {GR_METHOD_DIV_OTHER,   (gr_funcptr) vector_gr_vec_div_other},
    {GR_METHOD_OTHER_DIV,   (gr_funcptr) vector_gr_vec_other_div},

    {GR_METHOD_DIV_NONUNIQUE,   (gr_funcptr) vector_gr_vec_div_nonunique},
    {GR_METHOD_DIVIDES,         (gr_funcptr) vector_gr_vec_divides},

    {GR_METHOD_DIVEXACT,        (gr_funcptr) vector_gr_vec_divexact},
    {GR_METHOD_DIVEXACT_UI,     (gr_funcptr) vector_gr_vec_divexact_ui},
    {GR_METHOD_DIVEXACT_SI,     (gr_funcptr) vector_gr_vec_divexact_si},
    {GR_METHOD_DIVEXACT_FMPZ,   (gr_funcptr) vector_gr_vec_divexact_fmpz},
    {GR_METHOD_DIVEXACT_FMPQ,   (gr_funcptr) vector_gr_vec_divexact_fmpq},
    {GR_METHOD_DIVEXACT_OTHER,  (gr_funcptr) vector_gr_vec_divexact_other},
    {GR_METHOD_OTHER_DIVEXACT,  (gr_funcptr) vector_gr_vec_other_divexact},

    {GR_METHOD_POW,         (gr_funcptr) vector_gr_vec_pow},
    {GR_METHOD_POW_UI,      (gr_funcptr) vector_gr_vec_pow_ui},
    {GR_METHOD_POW_SI,      (gr_funcptr) vector_gr_vec_pow_si},
    {GR_METHOD_POW_FMPZ,    (gr_funcptr) vector_gr_vec_pow_fmpz},
    {GR_METHOD_POW_FMPQ,    (gr_funcptr) vector_gr_vec_pow_fmpq},
    {GR_METHOD_POW_OTHER,   (gr_funcptr) vector_gr_vec_pow_other},
    {GR_METHOD_OTHER_POW,   (gr_funcptr) vector_gr_vec_other_pow},

    {GR_METHOD_SQRT,            (gr_funcptr) vector_gr_vec_sqrt},
    {GR_METHOD_RSQRT,           (gr_funcptr) vector_gr_vec_rsqrt},
    {GR_METHOD_FLOOR,           (gr_funcptr) vector_gr_vec_floor},
    {GR_METHOD_CEIL,            (gr_funcptr) vector_gr_vec_ceil},
    {GR_METHOD_TRUNC,           (gr_funcptr) vector_gr_vec_trunc},
    {GR_METHOD_NINT,            (gr_funcptr) vector_gr_vec_nint},
    {GR_METHOD_ABS,             (gr_funcptr) vector_gr_vec_abs},
    {GR_METHOD_CONJ,            (gr_funcptr) vector_gr_vec_conj},
    {GR_METHOD_RE,              (gr_funcptr) vector_gr_vec_re},
    {GR_METHOD_IM,              (gr_funcptr) vector_gr_vec_im},
    {GR_METHOD_SGN,             (gr_funcptr) vector_gr_vec_sgn},
    {GR_METHOD_CSGN,            (gr_funcptr) vector_gr_vec_csgn},

    {GR_METHOD_I,             (gr_funcptr) vector_gr_vec_i},
    {GR_METHOD_PI,            (gr_funcptr) vector_gr_vec_pi},

    {GR_METHOD_EXP,            (gr_funcptr) vector_gr_vec_exp},
    {GR_METHOD_LOG,            (gr_funcptr) vector_gr_vec_log},

    {0,                     (gr_funcptr) NULL},
};

void
_gr_ctx_init_vector(gr_ctx_t ctx, gr_ctx_t base_ring, int all_sizes, slong n)
{
    ctx->which_ring = GR_CTX_GR_VEC;
    ctx->sizeof_elem = sizeof(gr_vec_struct);
    ctx->size_limit = WORD_MAX;

    if (n < 0)
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);

    VECTOR_CTX(ctx)->base_ring = (gr_ctx_struct *) base_ring;
    VECTOR_CTX(ctx)->all_sizes = all_sizes;
    VECTOR_CTX(ctx)->n = n;

    ctx->methods = _gr_vec_methods;

    if (!_gr_vec_methods_initialized)
    {
        gr_method_tab_init(_gr_vec_methods, _gr_vec_methods_input);
        _gr_vec_methods_initialized = 1;
    }
}

void gr_ctx_init_vector_gr_vec(gr_ctx_t ctx, gr_ctx_t base_ring)
{
    _gr_ctx_init_vector(ctx, base_ring, 1, 0);
}

void gr_ctx_init_vector_space_gr_vec(gr_ctx_t ctx, gr_ctx_t base_ring, slong n)
{
    _gr_ctx_init_vector(ctx, base_ring, 0, n);
}
