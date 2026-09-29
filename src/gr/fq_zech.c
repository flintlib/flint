/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include "ulong_extras.h"
#include "fmpz.h"
#include "fq_nmod.h"
#include "fq_zech.h"
#include "fq_zech_vec.h"
#include "fq_zech_poly.h"
#include "fq_zech_poly_factor.h"
#include "fq_zech_mat.h"
#include "nmod_poly.h"
#include "fmpz_mod_poly.h"
#include "gr.h"
#include "gr_poly.h"
#include "gr/impl.h"
#include "gr_vec.h"
#include "gr_mat.h"
#include "gr_generic.h"

#define FQ_CTX(ring_ctx) ((fq_zech_ctx_struct *)(GR_CTX_DATA_AS_PTR(ring_ctx)))

static const char * default_var = "a";

/* todo: lots of inlining */

static void
_gr_fq_zech_ctx_clear(gr_ctx_t ctx)
{
    fq_zech_ctx_clear(FQ_CTX(ctx));
    flint_free(GR_CTX_DATA_AS_PTR(ctx));
}

static int
_gr_fq_zech_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    status |= gr_stream_write(out, "GF(");
    status |= gr_stream_write_ui(out, fq_zech_ctx_prime(FQ_CTX(ctx)));
    status |= gr_stream_write(out, "^");
    status |= gr_stream_write_si(out, fq_zech_ctx_degree(FQ_CTX(ctx)));
    status |= gr_stream_write(out, ") (fq_zech)");
    return status;
}

static int _gr_fq_zech_ctx_set_gen_name(gr_ctx_t ctx, const char * s)
{
    slong len;
    len = strlen(s);

    FQ_CTX(ctx)->fq_nmod_ctx->var = flint_realloc(FQ_CTX(ctx)->fq_nmod_ctx->var, len + 1);
    memcpy(FQ_CTX(ctx)->fq_nmod_ctx->var, s, len + 1);
    return GR_SUCCESS;
}

static int _gr_fq_zech_ctx_set_gen_names(gr_ctx_t ctx, const char ** s)
{
    return _gr_fq_zech_ctx_set_gen_name(ctx, s[0]);
}

static int
_gr_fq_zech_ctx_gen_name(char ** name, slong i, gr_ctx_t ctx)
{
    if (i != 0)
        return GR_DOMAIN;

    char * var = FQ_CTX(ctx)->fq_nmod_ctx->var;
    size_t len = strlen(var);
    * name = flint_malloc(len + 1);
    if (* name == NULL)
        return GR_UNABLE;
    strncpy(* name, var, len + 1);

    return GR_SUCCESS;
}

static void
_gr_fq_zech_init(fq_zech_t x, const gr_ctx_t ctx)
{
    fq_zech_init(x, FQ_CTX(ctx));
}

static void
_gr_fq_zech_clear(fq_zech_t x, const gr_ctx_t ctx)
{
    fq_zech_clear(x, FQ_CTX(ctx));
}

static void
_gr_fq_zech_swap(fq_zech_t x, fq_zech_t y, const gr_ctx_t ctx)
{
    fq_zech_t t;
    *t = *x;
    *x = *y;
    *y = *t;
}

static void
_gr_fq_zech_set_shallow(fq_zech_t res, const fq_zech_t x, const gr_ctx_t ctx)
{
    *res = *x;
}

static int
_gr_fq_zech_randtest(fq_zech_t res, flint_rand_t state, const gr_ctx_t ctx)
{
    fq_zech_randtest(res, state, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_zech_write(gr_stream_t out, const fq_zech_t x, const gr_ctx_t ctx)
{
    return gr_stream_write_free(out, fq_zech_get_str_pretty(x, FQ_CTX(ctx)));
}

static int
_gr_fq_zech_zero(fq_zech_t x, const gr_ctx_t ctx)
{
    fq_zech_zero(x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_zech_one(fq_zech_t x, const gr_ctx_t ctx)
{
    fq_zech_one(x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_zech_set_si(fq_zech_t res, slong v, const gr_ctx_t ctx)
{
    fq_zech_set_si(res, v, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_zech_set_ui(fq_zech_t res, ulong v, const gr_ctx_t ctx)
{
    fq_zech_set_ui(res, v, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_zech_set_fmpz(fq_zech_t res, const fmpz_t v, const gr_ctx_t ctx)
{
    fq_zech_set_fmpz(res, v, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static truth_t
_gr_fq_zech_is_zero(const fq_zech_t x, const gr_ctx_t ctx)
{
    return fq_zech_is_zero(x, FQ_CTX(ctx)) ? T_TRUE : T_FALSE;
}

static truth_t
_gr_fq_zech_is_one(const fq_zech_t x, const gr_ctx_t ctx)
{
    return fq_zech_is_one(x, FQ_CTX(ctx)) ? T_TRUE : T_FALSE;
}

static truth_t
_gr_fq_zech_equal(const fq_zech_t x, const fq_zech_t y, const gr_ctx_t ctx)
{
    return fq_zech_equal(x, y, FQ_CTX(ctx)) ? T_TRUE : T_FALSE;
}

static int
_gr_fq_zech_set(fq_zech_t res, const fq_zech_t x, const gr_ctx_t ctx)
{
    fq_zech_set(res, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_zech_neg(fq_zech_t res, const fq_zech_t x, const gr_ctx_t ctx)
{
    fq_zech_neg(res, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_zech_add(fq_zech_t res, const fq_zech_t x, const fq_zech_t y, const gr_ctx_t ctx)
{
    fq_zech_add(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_zech_sub(fq_zech_t res, const fq_zech_t x, const fq_zech_t y, const gr_ctx_t ctx)
{
    fq_zech_sub(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_zech_mul(fq_zech_t res, const fq_zech_t x, const fq_zech_t y, const gr_ctx_t ctx)
{
    fq_zech_mul(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_zech_sqr(fq_zech_t res, const fq_zech_t x, const gr_ctx_t ctx)
{
    fq_zech_mul(res, x, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_zech_mul_two(fq_zech_t res, const fq_zech_t x, const gr_ctx_t ctx)
{
    fq_zech_add(res, x, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_zech_addmul(fq_zech_t res, const fq_zech_t x, const fq_zech_t y, const gr_ctx_t ctx)
{
    fq_zech_t t;
    fq_zech_mul(t, x, y, FQ_CTX(ctx));
    fq_zech_add(res, res, t, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_zech_submul(fq_zech_t res, const fq_zech_t x, const fq_zech_t y, const gr_ctx_t ctx)
{
    fq_zech_t t;
    fq_zech_mul(t, x, y, FQ_CTX(ctx));
    fq_zech_sub(res, res, t, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_zech_mul_si(fq_zech_t res, const fq_zech_t x, slong y, const gr_ctx_t ctx)
{
    fq_zech_mul_si(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

/* fq_zech_mul_ui inlined (used for example by _gr_poly_derivative) */
static int
_gr_fq_zech_mul_ui(fq_zech_t res, const fq_zech_t x, ulong y, const gr_ctx_t ctx)
{
    const fq_zech_ctx_struct * fctx = FQ_CTX(ctx);

    if (y >= fctx->p)
        y = n_mod2_precomp(y, fctx->p, fctx->ppre);

    if (y == 0 || x->value == fctx->qm1)
        res->value = fctx->qm1;
    else
        res->value = n_addmod(x->value, fctx->prime_field_table[y], fctx->qm1);

    return GR_SUCCESS;
}

static int
_gr_fq_zech_mul_fmpz(fq_zech_t res, const fq_zech_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    fq_zech_mul_fmpz(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_zech_inv(fq_zech_t res, const fq_zech_t x, const gr_ctx_t ctx)
{
    if (fq_zech_is_zero(x, FQ_CTX(ctx)))
    {
        return GR_DOMAIN;
    }
    else
    {
        fq_zech_inv(res, x, FQ_CTX(ctx));
        return GR_SUCCESS;
    }
}

static int
_gr_fq_zech_div(fq_zech_t res, const fq_zech_t x, const fq_zech_t y, const gr_ctx_t ctx)
{
    if (fq_zech_is_zero(y, FQ_CTX(ctx)))
    {
        return GR_DOMAIN;
    }
    else
    {
        fq_zech_t t;
        fq_zech_init(t, FQ_CTX(ctx));
        fq_zech_inv(t, y, FQ_CTX(ctx));
        fq_zech_mul(res, x, t, FQ_CTX(ctx));
        fq_zech_clear(t, FQ_CTX(ctx));
        return GR_SUCCESS;
    }
}

#if 0
static int
_gr_fq_zech_sqr(fq_zech_t res, const fq_zech_t x, const gr_ctx_t ctx)
{
    fq_zech_sqr(res, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_zech_pow_ui(fq_zech_t res, const fq_zech_t x, ulong y, const gr_ctx_t ctx)
{
    fq_zech_pow_ui(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_zech_pow_fmpz(fq_zech_t res, const fq_zech_t x, const fmpz_t y, gr_ctx_t ctx)
{
    if (fmpz_sgn(y) < 0)
    {
        return gr_generic_pow_fmpz(res, x, y, ctx);
    }
    else
    {
        fq_zech_pow(res, x, y, FQ_CTX(ctx));
        return GR_SUCCESS;
    }
}
#endif

static truth_t
_gr_fq_zech_is_invertible(const fq_zech_t x, const gr_ctx_t ctx)
{
    return (!fq_zech_is_zero(x, FQ_CTX(ctx))) ? T_TRUE : T_FALSE;
}

static truth_t
_gr_fq_zech_is_square(const fq_zech_t x, const gr_ctx_t ctx)
{
    return fq_zech_is_square(x, FQ_CTX(ctx)) ? T_TRUE : T_FALSE;
}

static int
_gr_fq_zech_sqrt(fq_zech_t res, const fq_zech_t x, const gr_ctx_t ctx)
{
    if (fq_zech_sqrt(res, x, FQ_CTX(ctx)))
    {
        return GR_SUCCESS;
    }
    else
    {
        return GR_DOMAIN;
    }
}

static int
_gr_ctx_fq_zech_prime(fmpz_t p, gr_ctx_t ctx)
{
    fmpz_set_ui(p, fq_zech_ctx_prime(FQ_CTX(ctx)));
    return GR_SUCCESS;
}

static int
_gr_ctx_fq_zech_degree(slong * deg, gr_ctx_t ctx)
{
    *deg = fq_zech_ctx_degree(FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_ctx_fq_zech_order(fmpz_t q, gr_ctx_t ctx)
{
    fmpz_set_ui(q, fq_zech_ctx_order_ui(FQ_CTX(ctx)));
    return GR_SUCCESS;
}

static int
_gr_fq_zech_gen(gr_ptr res, gr_ctx_t ctx)
{
    fq_zech_gen(res, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_zech_frobenius(gr_ptr res, gr_srcptr x, slong e, gr_ctx_t ctx)
{
    fq_zech_frobenius(res, x, e, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_zech_multiplicative_order(fmpz_t res, gr_srcptr x, gr_ctx_t ctx)
{
    int ret;
    ret = fq_zech_multiplicative_order(res, x, FQ_CTX(ctx));

    if (ret == 1)
        return GR_SUCCESS;

    /* todo: better solution? */
    return GR_DOMAIN;
}

static int
_gr_fq_zech_norm(fmpz_t res, gr_srcptr x, gr_ctx_t ctx)
{
    fq_zech_norm(res, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_zech_trace(fmpz_t res, gr_srcptr x, gr_ctx_t ctx)
{
    fq_zech_trace(res, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static truth_t
_gr_fq_zech_is_primitive(gr_srcptr x, gr_ctx_t ctx)
{
    return fq_zech_is_primitive(x, FQ_CTX(ctx)) ? T_TRUE : T_FALSE;
}

static int
_gr_fq_zech_pth_root(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    fq_zech_pth_root(res, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static void
_gr_fq_zech_vec_init(fq_zech_struct * vec, slong len, gr_ctx_t ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        fq_zech_init(vec + i, FQ_CTX(ctx));
}

static void
_gr_fq_zech_vec_clear(fq_zech_struct * vec, slong len, gr_ctx_t ctx)
{
}

static void
_gr_fq_zech_vec_swap(fq_zech_struct * vec1, fq_zech_struct * vec2, slong len, gr_ctx_t ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        fq_zech_swap(vec1 + i, vec2 + i, FQ_CTX(ctx));
}

static int
_gr_fq_zech_vec_set(fq_zech_struct * res, const fq_zech_struct * vec, slong len, gr_ctx_t ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        res[i].value = vec[i].value;

    return GR_SUCCESS;
}

static int
_gr_fq_zech_vec_normalise(slong * res, const fq_zech_struct * vec, slong len, gr_ctx_t ctx)
{
    while (len > 0 && fq_zech_is_zero(vec + len - 1, FQ_CTX(ctx)))
        len--;

    res[0] = len;
    return GR_SUCCESS;
}

static slong
_gr_fq_zech_vec_normalise_weak(const fq_zech_struct * vec, slong len, gr_ctx_t ctx)
{
    while (len > 0 && fq_zech_is_zero(vec + len - 1, FQ_CTX(ctx)))
        len--;

    return len;
}

static int
_gr_fq_zech_vec_mul_scalar(fq_zech_struct * res, const fq_zech_struct * vec, slong len, const fq_zech_t x, gr_ctx_t ctx)
{
    _fq_zech_vec_scalar_mul_fq_zech(res, vec, len, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

/* Zech logarithm of 2^c, or qm1 (the representation of zero) if p = 2. */
static ulong
_fq_zech_two_pow_si_log(slong c, const fq_zech_ctx_t ctx)
{
    ulong log2, e;

    if (ctx->p == 2)
        return ctx->qm1;

    log2 = ctx->prime_field_table[2];
    e = (c >= 0) ? (ulong) c : -(ulong) c;
    e = n_mulmod2(e % ctx->qm1, log2, ctx->qm1);

    return (c >= 0) ? e : n_negmod(e, ctx->qm1);
}

static int
_gr_fq_zech_mul_2exp_si(fq_zech_t res, const fq_zech_t x, slong c, const gr_ctx_t ctx)
{
    const fq_zech_ctx_struct * fctx = FQ_CTX(ctx);
    ulong t;

    if (c == 0)
    {
        res->value = x->value;
        return GR_SUCCESS;
    }

    t = _fq_zech_two_pow_si_log(c, fctx);

    if (t == fctx->qm1)   /* characteristic 2 */
    {
        if (c < 0)
            return GR_DOMAIN;

        res->value = fctx->qm1;
        return GR_SUCCESS;
    }

    res->value = (x->value == fctx->qm1) ? fctx->qm1 : n_addmod(x->value, t, fctx->qm1);
    return GR_SUCCESS;
}

static int
_gr_fq_zech_vec_mul_scalar_2exp_si(fq_zech_struct * res, const fq_zech_struct * vec, slong len, slong c, gr_ctx_t ctx)
{
    const fq_zech_ctx_struct * fctx = FQ_CTX(ctx);
    ulong t, qm1 = fctx->qm1;
    slong i;

    if (c == 0)
    {
        for (i = 0; i < len; i++)
            res[i] = vec[i];
        return GR_SUCCESS;
    }

    t = _fq_zech_two_pow_si_log(c, fctx);

    if (t == qm1)   /* characteristic 2 */
    {
        if (c < 0)
            return GR_DOMAIN;

        for (i = 0; i < len; i++)
            res[i].value = qm1;
        return GR_SUCCESS;
    }

    for (i = 0; i < len; i++)
        res[i].value = (vec[i].value == qm1) ? qm1 : n_addmod(vec[i].value, t, qm1);

    return GR_SUCCESS;
}

/* Dot products with the Zech logarithm arithmetic inlined. The zero
   element is represented by qm1. */

#define FQ_ZECH_ADD_LOG(s, v, qm1, table) \
    do { \
        if ((s) == (qm1)) \
            (s) = (v); \
        else if ((v) != (qm1)) \
        { \
            ulong __c = (table)[n_submod((s), (v), (qm1))]; \
            (s) = (__c == (qm1)) ? (qm1) : n_addmod(__c, (v), (qm1)); \
        } \
    } while (0)

#define FQ_ZECH_DOT_TERM(s, i, vec1, vec2, IDX2) \
    do { \
        ulong __a = (vec1)[i].value; \
        ulong __b = (vec2)[IDX2(i)].value; \
        if (__a != __qm1 && __b != __qm1) \
        { \
            __a = n_addmod(__a, __b, __qm1); \
            FQ_ZECH_ADD_LOG(s, __a, __qm1, __table); \
        } \
    } while (0)

/* Each Zech addition depends on a table lookup, so a single accumulator
   makes the loop latency-bound (especially when the table does not
   fit in L1 or L2 cache); use four independent accumulators for
   long dot products. This is done in separate functions to keep the
   code for short dot products lean. */
#ifndef FQ_ZECH_DOT_UNROLL_CUTOFF
#define FQ_ZECH_DOT_UNROLL_CUTOFF 8
#endif

#define FQ_ZECH_DOT_LOG_4(func, IDX2) \
FLINT_STATIC_NOINLINE ulong \
func(const fq_zech_struct * vec1, const fq_zech_struct * vec2, slong len, \
    const fq_zech_ctx_struct * fctx) \
{ \
    ulong __qm1 = fctx->qm1; \
    const ulong * __table = fctx->zech_log_table; \
    ulong __s, __s1, __s2, __s3; \
    slong __i; \
    __s = __s1 = __s2 = __s3 = __qm1; \
    for (__i = 0; __i + 4 <= len; __i += 4) \
    { \
        FQ_ZECH_DOT_TERM(__s, __i, vec1, vec2, IDX2); \
        FQ_ZECH_DOT_TERM(__s1, __i + 1, vec1, vec2, IDX2); \
        FQ_ZECH_DOT_TERM(__s2, __i + 2, vec1, vec2, IDX2); \
        FQ_ZECH_DOT_TERM(__s3, __i + 3, vec1, vec2, IDX2); \
    } \
    for ( ; __i < len; __i++) \
        FQ_ZECH_DOT_TERM(__s, __i, vec1, vec2, IDX2); \
    FQ_ZECH_ADD_LOG(__s, __s1, __qm1, __table); \
    FQ_ZECH_ADD_LOG(__s2, __s3, __qm1, __table); \
    FQ_ZECH_ADD_LOG(__s, __s2, __qm1, __table); \
    return __s; \
}

#define FQ_ZECH_VEC_DOT(res, initial, subtract, vec1, vec2, len, fctx, IDX2, func4) \
    do { \
        ulong __qm1 = (fctx)->qm1; \
        const ulong * __table = (fctx)->zech_log_table; \
        ulong __s; \
        slong __i; \
        if ((len) < FQ_ZECH_DOT_UNROLL_CUTOFF) \
        { \
            __s = __qm1; \
            for (__i = 0; __i < (len); __i++) \
                FQ_ZECH_DOT_TERM(__s, __i, vec1, vec2, IDX2); \
        } \
        else \
        { \
            __s = func4(vec1, vec2, len, fctx); \
        } \
        if (subtract && __s != __qm1) \
        { \
            __s += (fctx)->qm1o2; \
            if (__s >= __qm1) \
                __s -= __qm1; \
        } \
        if ((initial) != NULL) \
        { \
            ulong __t = (initial)->value; \
            FQ_ZECH_ADD_LOG(__s, __t, __qm1, __table); \
        } \
        (res)->value = __s; \
    } while (0)

/* res += vec * x with the logarithm of x given; the terms are
   independent, so this is throughput-bound rather than latency-bound. */
static void
_fq_zech_vec_addmul_log(fq_zech_struct * res, const fq_zech_struct * vec,
    slong len, ulong xl, const fq_zech_ctx_struct * fctx)
{
    ulong qm1 = fctx->qm1;
    const ulong * table = fctx->zech_log_table;
    ulong a, s;
    slong i;

    if (xl == qm1)
        return;

    for (i = 0; i < len; i++)
    {
        a = vec[i].value;
        if (a != qm1)
        {
            a = n_addmod(a, xl, qm1);
            s = res[i].value;
            FQ_ZECH_ADD_LOG(s, a, qm1, table);
            res[i].value = s;
        }
    }
}

static int
_gr_fq_zech_vec_addmul_scalar(fq_zech_struct * res, const fq_zech_struct * vec, slong len, const fq_zech_t x, gr_ctx_t ctx)
{
    _fq_zech_vec_addmul_log(res, vec, len, x->value, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_zech_vec_submul_scalar(fq_zech_struct * res, const fq_zech_struct * vec, slong len, const fq_zech_t x, gr_ctx_t ctx)
{
    const fq_zech_ctx_struct * fctx = FQ_CTX(ctx);
    ulong xl = x->value;

    /* multiply x by -1 */
    if (xl != fctx->qm1)
    {
        xl += fctx->qm1o2;
        if (xl >= fctx->qm1)
            xl -= fctx->qm1;
    }

    _fq_zech_vec_addmul_log(res, vec, len, xl, fctx);
    return GR_SUCCESS;
}

#define FQ_ZECH_IDX_FWD(i) (i)
#define FQ_ZECH_IDX_REV(i) ((len) - 1 - (i))

FQ_ZECH_DOT_LOG_4(_fq_zech_dot_log_4, FQ_ZECH_IDX_FWD)
FQ_ZECH_DOT_LOG_4(_fq_zech_dot_rev_log_4, FQ_ZECH_IDX_REV)

static int
_gr_fq_zech_vec_dot(fq_zech_struct * res, const fq_zech_struct * initial, int subtract, const fq_zech_struct * vec1, const fq_zech_struct * vec2, slong len, gr_ctx_t ctx)
{
    FQ_ZECH_VEC_DOT(res, initial, subtract, vec1, vec2, len, FQ_CTX(ctx), FQ_ZECH_IDX_FWD, _fq_zech_dot_log_4);
    return GR_SUCCESS;
}

static int
_gr_fq_zech_vec_dot_rev(fq_zech_struct * res, const fq_zech_struct * initial, int subtract, const fq_zech_struct * vec1, const fq_zech_struct * vec2, slong len, gr_ctx_t ctx)
{
    FQ_ZECH_VEC_DOT(res, initial, subtract, vec1, vec2, len, FQ_CTX(ctx), FQ_ZECH_IDX_REV, _fq_zech_dot_rev_log_4);
    return GR_SUCCESS;
}

#undef FQ_ZECH_IDX_FWD
#undef FQ_ZECH_IDX_REV
#undef FQ_ZECH_VEC_DOT
#undef FQ_ZECH_DOT_LOG_4
#undef FQ_ZECH_DOT_TERM
#undef FQ_ZECH_ADD_LOG


/* gcd and xgcd of the fq_zech_poly module */
static int
_gr_fq_zech_poly_gcd(fq_zech_struct * G, slong * lenG, const fq_zech_struct * A, slong lenA, const fq_zech_struct * B, slong lenB, gr_ctx_t ctx)
{
    *lenG = _fq_zech_poly_gcd(G, A, lenA, B, lenB, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_zech_poly_xgcd(slong * lenG, fq_zech_struct * G, fq_zech_struct * S, fq_zech_struct * T, const fq_zech_struct * A, slong lenA, const fq_zech_struct * B, slong lenB, gr_ctx_t ctx)
{
    *lenG = _fq_zech_poly_xgcd(G, S, T, A, lenA, B, lenB, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_zech_poly_mullow(fq_zech_struct * res,
    const fq_zech_struct * poly1, slong len1,
    const fq_zech_struct * poly2, slong len2, slong n, gr_ctx_t ctx)
{
    if (len1 + len2 - 1 == n)
    {
        if (poly1 == poly2 && len1 == len2)
            _fq_zech_poly_sqr(res, poly1, len1, FQ_CTX(ctx));
        else if (len1 >= len2)
            _fq_zech_poly_mul(res, poly1, len1, poly2, len2, FQ_CTX(ctx));
        else
            _fq_zech_poly_mul(res, poly2, len2, poly1, len1, FQ_CTX(ctx));
    }
    else
    {
        if (len1 >= len2)
            _fq_zech_poly_mullow(res, poly1, len1, poly2, len2, n, FQ_CTX(ctx));
        else
            _fq_zech_poly_mullow(res, poly2, len2, poly1, len1, n, FQ_CTX(ctx));
    }

    return GR_SUCCESS;
}

static int
_gr_fq_zech_poly_mulmid(fq_zech_struct * res,
    const fq_zech_struct * poly1, slong len1,
    const fq_zech_struct * poly2, slong len2, slong nlo, slong nhi, gr_ctx_t ctx)
{
    if (nlo == 0)
        return _gr_fq_zech_poly_mullow(res, poly1, len1, poly2, len2, nhi, ctx);

    if (FLINT_MIN(len1, len2) >= FQ_ZECH_POLY_MUL_UNIVARIATE_MIN_LEN(FQ_CTX(ctx)) &&
        _fq_zech_poly_mulmid_want_univariate(len1, len2, nlo, nhi, FQ_CTX(ctx)))
        _fq_zech_poly_mulmid_univariate(res, poly1, len1, poly2, len2, nlo, nhi, FQ_CTX(ctx));
    else
        return _gr_poly_mulmid_classical(res, poly1, len1, poly2, len2, nlo, nhi, ctx);

    return GR_SUCCESS;
}

/* With the inlined Zech arithmetic, basecase division is fast; Newton
   division only wins when both the quotient and the divisor are long
   enough for the products to use univariate multiplication with a good
   margin. */
static int
_gr_fq_zech_poly_divrem(fq_zech_struct * Q, fq_zech_struct * R,
    const fq_zech_struct * A, slong lenA,
    const fq_zech_struct * B, slong lenB, gr_ctx_t ctx)
{
    slong lenQ = lenA - lenB + 1;

    if (FLINT_MIN(lenQ, lenB) >= 4 * FQ_ZECH_POLY_MUL_UNIVARIATE_MIN_LEN(FQ_CTX(ctx)) &&
        FLINT_MIN(lenQ, lenB) >= 4 * _fq_zech_poly_mul_univariate_threshold(FQ_CTX(ctx)))
        return _gr_poly_divrem_newton(Q, R, A, lenA, B, lenB, ctx);
    else
        return _gr_poly_divrem_basecase(Q, R, A, lenA, B, lenB, ctx);
}

/* Horner's rule with the Zech logarithm arithmetic inlined. The zero
   element is represented by qm1. */
static int
_gr_fq_zech_poly_evaluate(fq_zech_t res, const fq_zech_struct * f, slong len,
    const fq_zech_t x, gr_ctx_t ctx)
{
    const fq_zech_ctx_struct * fctx = FQ_CTX(ctx);
    ulong qm1 = fctx->qm1;
    const ulong * table = fctx->zech_log_table;
    ulong s, a, c, xv;
    slong i;

    if (len == 0)
    {
        res->value = qm1;
        return GR_SUCCESS;
    }

    xv = x->value;

    if (len == 1 || xv == qm1)
    {
        res->value = f[0].value;
        return GR_SUCCESS;
    }

    s = f[len - 1].value;

    for (i = len - 2; i >= 0; i--)
    {
        /* s = s * x */
        if (s != qm1)
            s = n_addmod(s, xv, qm1);

        /* s = s + f[i] */
        a = f[i].value;
        if (s == qm1)
        {
            s = a;
        }
        else if (a != qm1)
        {
            c = table[n_submod(s, a, qm1)];
            s = (c == qm1) ? qm1 : n_addmod(c, a, qm1);
        }
    }

    res->value = s;
    return GR_SUCCESS;
}

/* todo: also need the _other version ... ? */
/* todo: implement generically */

static int
_gr_fq_zech_roots_gr_poly(gr_vec_t roots, gr_vec_t mult, const fq_zech_poly_t poly, int flags, gr_ctx_t ctx)
{
    if (poly->length == 0)
        return GR_DOMAIN;

    {
        gr_ctx_t ZZ;
        fq_zech_poly_factor_t fac;
        slong i, num;

        gr_ctx_init_fmpz(ZZ);
        fq_zech_poly_factor_init(fac, FQ_CTX(ctx));
        fq_zech_poly_roots(fac, poly, 1, FQ_CTX(ctx));

        num = fac->num;

        gr_vec_set_length(roots, num, ctx);
        gr_vec_set_length(mult, num, ZZ);

        for (i = 0; i < num; i++)
        {
            fq_zech_neg(gr_vec_entry_ptr(roots, i, ctx), fac->poly[i].coeffs, FQ_CTX(ctx));

            /* work around flint bug: factors can be non-monic */
            if (!fq_zech_is_one(fac->poly[i].coeffs + 1, FQ_CTX(ctx)))
                fq_zech_div(gr_vec_entry_ptr(roots, i, ctx), gr_vec_entry_ptr(roots, i, ctx), fac->poly[i].coeffs + 1, FQ_CTX(ctx));

            fmpz_set_ui(((fmpz *) mult->entries) + i, fac->exp[i]);
        }

        fq_zech_poly_factor_clear(fac, FQ_CTX(ctx));
        gr_ctx_clear(ZZ);
    }

    return GR_SUCCESS;
}

static int
_gr_fq_zech_mat_mul(fq_zech_mat_t res, const fq_zech_mat_t x, const fq_zech_mat_t y, gr_ctx_t ctx)
{
    fq_zech_mat_mul(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_zech_mat_nonsingular_solve_tril(fq_zech_mat_t X, const fq_zech_mat_t L, const fq_zech_mat_t B, int unit, gr_ctx_t ctx)
{
    if (B->r < 64 || B->c < 64)
        return gr_mat_nonsingular_solve_tril_classical((gr_mat_struct *) X, (const gr_mat_struct *) L, (const gr_mat_struct *) B, unit, ctx);
    else
        return gr_mat_nonsingular_solve_tril_recursive((gr_mat_struct *) X, (const gr_mat_struct *) L, (const gr_mat_struct *) B, unit, ctx);
}

static int
_gr_fq_zech_mat_nonsingular_solve_triu(fq_zech_mat_t X, const fq_zech_mat_t U, const fq_zech_mat_t B, int unit, gr_ctx_t ctx)
{
    if (B->r < 64 || B->c < 64)
        return gr_mat_nonsingular_solve_triu_classical((gr_mat_struct *) X, (const gr_mat_struct *) U, (const gr_mat_struct *) B, unit, ctx);
    else
        return gr_mat_nonsingular_solve_triu_recursive((gr_mat_struct *) X, (const gr_mat_struct *) U, (const gr_mat_struct *) B, unit, ctx);
}

static int
_gr_fq_zech_mat_charpoly(fq_zech_struct * res, const fq_zech_mat_t mat, gr_ctx_t ctx)
{
    slong n = mat->r;

    if (n <= 4)
        return _gr_mat_charpoly_berkowitz(res, (const gr_mat_struct *) mat, ctx);
    else
        return _gr_mat_charpoly_danilevsky(res, (const gr_mat_struct *) mat, ctx);
}

/* Vector methods avoiding per-element dispatch through the method table */

static int
_gr_fq_zech_vec_zero(fq_zech_struct * res, slong len, gr_ctx_t ctx)
{
    ulong qm1 = FQ_CTX(ctx)->qm1;
    slong i;

    for (i = 0; i < len; i++)
        res[i].value = qm1;

    return GR_SUCCESS;
}

static int
_gr_fq_zech_vec_add(fq_zech_struct * res, const fq_zech_struct * vec1, const fq_zech_struct * vec2, slong len, gr_ctx_t ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        fq_zech_add(res + i, vec1 + i, vec2 + i, FQ_CTX(ctx));

    return GR_SUCCESS;
}

static int
_gr_fq_zech_vec_sub(fq_zech_struct * res, const fq_zech_struct * vec1, const fq_zech_struct * vec2, slong len, gr_ctx_t ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        fq_zech_sub(res + i, vec1 + i, vec2 + i, FQ_CTX(ctx));

    return GR_SUCCESS;
}

static int
_gr_fq_zech_vec_add_scalar(fq_zech_struct * res, const fq_zech_struct * vec, slong len, const fq_zech_struct * c, gr_ctx_t ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        fq_zech_add(res + i, vec + i, c, FQ_CTX(ctx));

    return GR_SUCCESS;
}

static truth_t
_gr_fq_zech_vec_is_zero(const fq_zech_struct * vec, slong len, gr_ctx_t ctx)
{
    ulong qm1 = FQ_CTX(ctx)->qm1;
    slong i;

    for (i = 0; i < len; i++)
        if (vec[i].value != qm1)
            return T_FALSE;

    return T_TRUE;
}

static truth_t
_gr_fq_zech_vec_equal(const fq_zech_struct * vec1, const fq_zech_struct * vec2, slong len, gr_ctx_t ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        if (vec1[i].value != vec2[i].value)
            return T_FALSE;

    return T_TRUE;
}

int _fq_zech_methods_initialized = 0;

gr_static_method_table _fq_zech_methods;

gr_method_tab_input _fq_zech_methods_input[] =
{
    {GR_METHOD_CTX_CLEAR,       (gr_funcptr) _gr_fq_zech_ctx_clear},
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _gr_fq_zech_ctx_write},
    {GR_METHOD_CTX_SET_GEN_NAME,    (gr_funcptr) _gr_fq_zech_ctx_set_gen_name},
    {GR_METHOD_CTX_SET_GEN_NAMES,   (gr_funcptr) _gr_fq_zech_ctx_set_gen_names},
    {GR_METHOD_CTX_NGENS,       (gr_funcptr) gr_generic_ctx_ngens_1},
    {GR_METHOD_CTX_GEN_NAME,    (gr_funcptr) _gr_fq_zech_ctx_gen_name},
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FINITE,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_ALGEBRAICALLY_CLOSED,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ORDERED_RING,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_EXACT,    (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_CANONICAL,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_INIT,            (gr_funcptr) _gr_fq_zech_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) _gr_fq_zech_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) _gr_fq_zech_swap},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) _gr_fq_zech_set_shallow},
    {GR_METHOD_RANDTEST,        (gr_funcptr) _gr_fq_zech_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) _gr_fq_zech_write},
    {GR_METHOD_ZERO,            (gr_funcptr) _gr_fq_zech_zero},
    {GR_METHOD_ONE,             (gr_funcptr) _gr_fq_zech_one},
    {GR_METHOD_GEN,             (gr_funcptr) _gr_fq_zech_gen},
    {GR_METHOD_GENS,            (gr_funcptr) gr_generic_gens_single},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) _gr_fq_zech_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) _gr_fq_zech_is_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) _gr_fq_zech_equal},
    {GR_METHOD_SET,             (gr_funcptr) _gr_fq_zech_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) _gr_fq_zech_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) _gr_fq_zech_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) _gr_fq_zech_set_fmpz},
    {GR_METHOD_NEG,             (gr_funcptr) _gr_fq_zech_neg},
    {GR_METHOD_ADD,             (gr_funcptr) _gr_fq_zech_add},
    {GR_METHOD_SUB,             (gr_funcptr) _gr_fq_zech_sub},
    {GR_METHOD_MUL,             (gr_funcptr) _gr_fq_zech_mul},
    {GR_METHOD_ADDMUL,          (gr_funcptr) _gr_fq_zech_addmul},
    {GR_METHOD_SUBMUL,          (gr_funcptr) _gr_fq_zech_submul},
    {GR_METHOD_SQR,             (gr_funcptr) _gr_fq_zech_sqr},
    {GR_METHOD_MUL_TWO,         (gr_funcptr) _gr_fq_zech_mul_two},
    {GR_METHOD_MUL_UI,          (gr_funcptr) _gr_fq_zech_mul_ui},
    {GR_METHOD_MUL_SI,          (gr_funcptr) _gr_fq_zech_mul_si},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) _gr_fq_zech_mul_fmpz},
    {GR_METHOD_MUL_2EXP_SI,     (gr_funcptr) _gr_fq_zech_mul_2exp_si},
    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _gr_fq_zech_is_invertible},
    {GR_METHOD_INV,             (gr_funcptr) _gr_fq_zech_inv},
    {GR_METHOD_DIV,             (gr_funcptr) _gr_fq_zech_div},
    {GR_METHOD_IS_SQUARE,       (gr_funcptr) _gr_fq_zech_is_square},
    {GR_METHOD_SQRT,            (gr_funcptr) _gr_fq_zech_sqrt},

    {GR_METHOD_IS_SQUARE,       (gr_funcptr) _gr_fq_zech_is_square},
    {GR_METHOD_SQRT,            (gr_funcptr) _gr_fq_zech_sqrt},

    {GR_METHOD_CTX_FQ_PRIME,            (gr_funcptr) _gr_ctx_fq_zech_prime},
    {GR_METHOD_CTX_FQ_DEGREE,           (gr_funcptr) _gr_ctx_fq_zech_degree},
    {GR_METHOD_CTX_FQ_ORDER,            (gr_funcptr) _gr_ctx_fq_zech_order},
    {GR_METHOD_FQ_FROBENIUS,            (gr_funcptr) _gr_fq_zech_frobenius},
    {GR_METHOD_FQ_MULTIPLICATIVE_ORDER, (gr_funcptr) _gr_fq_zech_multiplicative_order},
    {GR_METHOD_FQ_NORM,                 (gr_funcptr) _gr_fq_zech_norm},
    {GR_METHOD_FQ_TRACE,                (gr_funcptr) _gr_fq_zech_trace},
    {GR_METHOD_FQ_IS_PRIMITIVE,         (gr_funcptr) _gr_fq_zech_is_primitive},
    {GR_METHOD_FQ_PTH_ROOT,             (gr_funcptr) _gr_fq_zech_pth_root},

    {GR_METHOD_VEC_INIT,            (gr_funcptr) _gr_fq_zech_vec_init},
    {GR_METHOD_VEC_CLEAR,           (gr_funcptr) _gr_fq_zech_vec_clear},
    {GR_METHOD_VEC_SET,             (gr_funcptr) _gr_fq_zech_vec_set},
    {GR_METHOD_VEC_ZERO,            (gr_funcptr) _gr_fq_zech_vec_zero},
    {GR_METHOD_VEC_ADD,             (gr_funcptr) _gr_fq_zech_vec_add},
    {GR_METHOD_VEC_SUB,             (gr_funcptr) _gr_fq_zech_vec_sub},
    {GR_METHOD_VEC_ADD_SCALAR,      (gr_funcptr) _gr_fq_zech_vec_add_scalar},
    {GR_METHOD_VEC_IS_ZERO,         (gr_funcptr) _gr_fq_zech_vec_is_zero},
    {GR_METHOD_VEC_EQUAL,           (gr_funcptr) _gr_fq_zech_vec_equal},
    {GR_METHOD_VEC_SWAP,            (gr_funcptr) _gr_fq_zech_vec_swap},
    {GR_METHOD_VEC_NORMALISE,       (gr_funcptr) _gr_fq_zech_vec_normalise},
    {GR_METHOD_VEC_NORMALISE_WEAK,  (gr_funcptr) _gr_fq_zech_vec_normalise_weak},
    {GR_METHOD_VEC_MUL_SCALAR,            (gr_funcptr) _gr_fq_zech_vec_mul_scalar},
    {GR_METHOD_VEC_ADDMUL_SCALAR,            (gr_funcptr) _gr_fq_zech_vec_addmul_scalar},
    {GR_METHOD_VEC_SUBMUL_SCALAR,            (gr_funcptr) _gr_fq_zech_vec_submul_scalar},
    {GR_METHOD_VEC_MUL_SCALAR_2EXP_SI,       (gr_funcptr) _gr_fq_zech_vec_mul_scalar_2exp_si},
    {GR_METHOD_VEC_DOT,                      (gr_funcptr) _gr_fq_zech_vec_dot},
    {GR_METHOD_VEC_DOT_REV,                  (gr_funcptr) _gr_fq_zech_vec_dot_rev},

    {GR_METHOD_POLY_MULLOW,     (gr_funcptr) _gr_fq_zech_poly_mullow},
    {GR_METHOD_POLY_MULMID,     (gr_funcptr) _gr_fq_zech_poly_mulmid},
    {GR_METHOD_POLY_DIVREM,     (gr_funcptr) _gr_fq_zech_poly_divrem},
    {GR_METHOD_POLY_EVALUATE,   (gr_funcptr) _gr_fq_zech_poly_evaluate},
    {GR_METHOD_POLY_GCD,        (gr_funcptr) _gr_fq_zech_poly_gcd},
    {GR_METHOD_POLY_XGCD,       (gr_funcptr) _gr_fq_zech_poly_xgcd},

    {GR_METHOD_POLY_FACTOR,     (gr_funcptr) _gr_poly_factor_finite_field_method},
    {GR_METHOD_POLY_ROOTS,      (gr_funcptr) _gr_fq_zech_roots_gr_poly},

    {GR_METHOD_MAT_MUL,         (gr_funcptr) _gr_fq_zech_mat_mul},
    {GR_METHOD_MAT_NONSINGULAR_SOLVE_TRIL,      (gr_funcptr) _gr_fq_zech_mat_nonsingular_solve_tril},
    {GR_METHOD_MAT_NONSINGULAR_SOLVE_TRIU,      (gr_funcptr) _gr_fq_zech_mat_nonsingular_solve_triu},
    {GR_METHOD_MAT_CHARPOLY,    (gr_funcptr) _gr_fq_zech_mat_charpoly},
    {0,                         (gr_funcptr) NULL},
};

void
_gr_ctx_init_fq_zech_from_ref(gr_ctx_t ctx, const void * fq_zech_ctx)
{
    ctx->which_ring = GR_CTX_FQ_ZECH;
    ctx->sizeof_elem = sizeof(fq_zech_struct);
    GR_CTX_DATA_AS_PTR(ctx) = (fq_zech_ctx_struct *) fq_zech_ctx;
    ctx->size_limit = WORD_MAX;
    ctx->methods = _fq_zech_methods;

    if (!_fq_zech_methods_initialized)
    {
        gr_method_tab_init(_fq_zech_methods, _fq_zech_methods_input);
        _fq_zech_methods_initialized = 1;
    }
}

void
gr_ctx_init_fq_zech(gr_ctx_t ctx, ulong p, slong d, const char * var)
{
    fq_zech_ctx_struct * fq_zech_ctx;

    fq_zech_ctx = flint_malloc(sizeof(fq_zech_ctx_struct));
    fq_zech_ctx_init_ui(fq_zech_ctx, p, d, var == NULL ? default_var : var);

    _gr_ctx_init_fq_zech_from_ref(ctx, fq_zech_ctx);
}

int
gr_ctx_init_fq_zech_modulus_nmod_poly(gr_ctx_t ctx, const nmod_poly_t modulus, const char * var)
{
    fq_zech_ctx_struct * fq_zech_ctx;
    fq_nmod_ctx_struct * fq_nmod_ctx;

    fq_nmod_ctx = flint_malloc(sizeof(fq_nmod_ctx_struct));
    fq_zech_ctx = flint_malloc(sizeof(fq_zech_ctx_struct));

    fq_nmod_ctx_init_modulus(fq_nmod_ctx, modulus, var == NULL ? default_var : var);

    if (fq_zech_ctx_init_fq_nmod_ctx_check(fq_zech_ctx, fq_nmod_ctx))
    {
        fq_zech_ctx->owns_fq_nmod_ctx = 1;
        _gr_ctx_init_fq_zech_from_ref(ctx, fq_zech_ctx);
        return GR_SUCCESS;
    }
    else
    {
        fq_nmod_ctx_clear(fq_nmod_ctx);
        flint_free(fq_zech_ctx);
        flint_free(fq_nmod_ctx);
        return GR_DOMAIN;
    }
}

int
gr_ctx_init_fq_zech_modulus_fmpz_mod_poly(gr_ctx_t ctx, const fmpz_mod_poly_t modulus, fmpz_mod_ctx_t mod_ctx, const char * var)
{
    nmod_poly_t nmodulus;
    int status;

    if (!fmpz_abs_fits_ui(mod_ctx->n))
        return GR_UNABLE;

    nmod_poly_init(nmodulus, fmpz_get_ui(mod_ctx->n));
    fmpz_mod_poly_get_nmod_poly(nmodulus, modulus);
    status = gr_ctx_init_fq_zech_modulus_nmod_poly(ctx, nmodulus, var);
    nmod_poly_clear(nmodulus);
    return status;
}
