/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "nmod_mat.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_mat.h"
#include "gr_poly.h"
#include "gr_generic.h"

static void
_gr_nmod_redc_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Integers mod ");
    gr_stream_write_ui(out, GR_NMOD_REDC_N(ctx));
    if (GR_NMOD_REDC_IS_FAST(ctx))
        gr_stream_write(out, " (nmod_redc_fast)");
    else
        gr_stream_write(out, " (nmod_redc)");
}

static truth_t
_gr_nmod_redc_ctx_is_field(gr_ctx_t ctx)
{
    if (GR_NMOD_REDC_IS_PRIME(ctx) != T_UNKNOWN)
        return GR_NMOD_REDC_IS_PRIME(ctx);

    return n_is_prime(GR_NMOD_REDC_N(ctx)) ? T_TRUE : T_FALSE;
}

static int
_gr_nmod_redc_ctx_set_is_field(gr_ctx_t ctx, truth_t is_field)
{
    GR_NMOD_REDC_IS_PRIME(ctx) = is_field;
    return GR_SUCCESS;
}

static void
_gr_nmod_redc_init(ulong * x, gr_ctx_t FLINT_UNUSED(ctx))
{
    x[0] = 0;
}

static void
_gr_nmod_redc_clear(ulong * FLINT_UNUSED(x), gr_ctx_t FLINT_UNUSED(ctx))
{
}

static void
_gr_nmod_redc_swap(ulong * x, ulong * y, gr_ctx_t FLINT_UNUSED(ctx))
{
    ulong t;
    t = *x;
    *x = *y;
    *y = t;
}

static void
_gr_nmod_redc_set_shallow(ulong * res, const ulong * x, gr_ctx_t FLINT_UNUSED(ctx))
{
    *res = *x;
}

static int
_gr_nmod_redc_randtest(ulong * res, flint_rand_t state, gr_ctx_t ctx)
{
    res[0] = n_randtest(state) % (GR_NMOD_REDC_IS_FAST(ctx) ? 2 * GR_NMOD_REDC_N(ctx) : GR_NMOD_REDC_N(ctx));
    return GR_SUCCESS;
}

static int
_gr_nmod_redc_write(gr_stream_t out, const ulong * x, gr_ctx_t ctx)
{
    gr_stream_write_ui(out, nmod_redc_get_nmod(x[0], GR_NMOD_REDC_CTX(ctx)));
    return GR_SUCCESS;
}

static int
_gr_nmod_redc_zero(ulong * x, gr_ctx_t FLINT_UNUSED(ctx))
{
    x[0] = 0;
    return GR_SUCCESS;
}

static int
_gr_nmod_redc_one(ulong * x, gr_ctx_t ctx)
{
    x[0] = GR_NMOD_REDC_ONE(ctx);
    return GR_SUCCESS;
}

static ulong
nmod_redc_set_si(slong v, nmod_redc_ctx_t ctx)
{
    ulong w = nmod_redc_set_ui(FLINT_UABS(v), ctx);
    if (v < 0)
        return nmod_redc_neg(w, ctx);
    else
        return w;
}

static ulong
nmod_redc_set_fmpz(const fmpz_t x, nmod_redc_ctx_t ctx)
{
    /* todo: conversion in redc form */
    ulong v = fmpz_get_nmod(x, GR_NMOD_REDC_MOD(ctx));
    return nmod_redc_set_nmod(v, GR_NMOD_REDC_CTX(ctx));
}

static int
_gr_nmod_redc_set_si(ulong * res, slong v, gr_ctx_t ctx)
{
    res[0] = nmod_redc_set_si(v, GR_NMOD_REDC_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_nmod_redc_set_ui(ulong * res, ulong v, gr_ctx_t ctx)
{
    res[0] = nmod_redc_set_ui(v, GR_NMOD_REDC_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_nmod_redc_set_fmpz(ulong * res, const fmpz_t x, gr_ctx_t ctx)
{
    res[0] = nmod_redc_set_fmpz(x, GR_NMOD_REDC_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_nmod_redc_get_fmpz(fmpz_t res, const ulong * x, const gr_ctx_t ctx)
{
    /* For now handle both redc and redc_fast in this method */
    ulong v = nmod_redc_fast_normalise(x[0], GR_NMOD_REDC_CTX(ctx));
    v = nmod_redc_get_nmod(v, GR_NMOD_REDC_CTX(ctx));
    fmpz_set_ui(res, v);
    return GR_SUCCESS;
}

static int
_gr_nmod_redc_inv(ulong * res, const ulong * x, const gr_ctx_t ctx)
{
    ulong v, w, r, g;

    v = x[0];

    /* For now handle both redc and redc_fast in this method */
    v = nmod_redc_fast_normalise(v, GR_NMOD_REDC_CTX(ctx));

    /* todo: also handle -1 fast? */
    if (v == GR_NMOD_REDC_ONE(ctx))
    {
        res[0] = v;
        return GR_SUCCESS;
    }

    /* Todo: maybe precompute redc(R) and do the inversion in Montgomery space
       to save a conversion. */
    w = nmod_redc_get_nmod(v, GR_NMOD_REDC_CTX(ctx));

    g = n_gcdinv(&r, w, GR_NMOD_REDC_N(ctx));

    if (g == 1)
    {
        res[0] = nmod_redc_set_nmod(r, GR_NMOD_REDC_CTX(ctx));
        return GR_SUCCESS;
    }
    else
    {
        res[0] = 0;
        return GR_DOMAIN;
    }
}

static int
_gr_nmod_redc_set_other(ulong * res, gr_ptr v, gr_ctx_t v_ctx, gr_ctx_t ctx)
{
    if (v_ctx->which_ring == GR_CTX_NMOD)
    {
        if (GR_NMOD_REDC_N(ctx) != NMOD_CTX(v_ctx).n)
            return GR_DOMAIN;

        res[0] = nmod_redc_set_nmod(((ulong *) v)[0], GR_NMOD_REDC_CTX(ctx));
        return GR_SUCCESS;
    }

    if (v_ctx->which_ring == GR_CTX_NMOD_REDC || v_ctx->which_ring == GR_CTX_NMOD_REDC_FAST)
    {
        if (NMOD_CTX(ctx).n != GR_NMOD_REDC_N(v_ctx))
            return GR_DOMAIN;

        ulong c = ((ulong *) v)[0];

        if (!GR_NMOD_REDC_IS_FAST(ctx))
            c = nmod_redc_fast_normalise(c, GR_NMOD_REDC_CTX(ctx));

        res[0] = c;
        return GR_SUCCESS;
    }

    {
        ulong c;
        int status;
        gr_ctx_t nmod_ctx;
        _gr_ctx_init_nmod(nmod_ctx, &GR_NMOD_REDC_MOD(ctx));

        status = gr_set_other(&c, v, v_ctx, nmod_ctx);

        if (status != GR_SUCCESS)
            return status;

        res[0] = nmod_redc_set_nmod(c, GR_NMOD_REDC_CTX(ctx));
        return GR_SUCCESS;
    }

}

static truth_t
_gr_nmod_redc_is_zero(const ulong * x, gr_ctx_t FLINT_UNUSED(ctx))
{
    return (x[0] == 0) ? T_TRUE : T_FALSE;
}

static truth_t
_gr_nmod_redc_fast_is_zero(const ulong * x, gr_ctx_t ctx)
{
    return (x[0] == 0 || x[0] == GR_NMOD_REDC_N(ctx)) ? T_TRUE : T_FALSE;
}

static truth_t
_gr_nmod_redc_is_one(const ulong * x, gr_ctx_t ctx)
{
    return (x[0] == GR_NMOD_REDC_ONE(ctx)) ? T_TRUE : T_FALSE;
}

static truth_t
_gr_nmod_redc_fast_is_one(const ulong * x, gr_ctx_t ctx)
{
    return (x[0] == GR_NMOD_REDC_ONE(ctx) || x[0] == GR_NMOD_REDC_ONE(ctx) + GR_NMOD_REDC_N(ctx)) ? T_TRUE : T_FALSE;
}

static truth_t
_gr_nmod_redc_is_neg_one(const ulong * x, gr_ctx_t ctx)
{
    return (x[0] == GR_NMOD_REDC_N(ctx) - GR_NMOD_REDC_ONE(ctx)) ? T_TRUE : T_FALSE;
}

static truth_t
_gr_nmod_redc_fast_is_neg_one(const ulong * x, gr_ctx_t ctx)
{
    return (x[0] == GR_NMOD_REDC_N(ctx) - GR_NMOD_REDC_ONE(ctx) ||
            x[0] == 2 * GR_NMOD_REDC_N(ctx) - GR_NMOD_REDC_ONE(ctx)) ? T_TRUE : T_FALSE;
}

static truth_t
_gr_nmod_redc_equal(const ulong * x, const ulong * y, gr_ctx_t FLINT_UNUSED(ctx))
{
    return (x[0] == y[0]) ? T_TRUE : T_FALSE;
}

static truth_t
_gr_nmod_redc_fast_equal(const ulong * x, const ulong * y, gr_ctx_t ctx)
{
    return (nmod_redc_fast_normalise(x[0], GR_NMOD_REDC_CTX(ctx)) == nmod_redc_fast_normalise(y[0], GR_NMOD_REDC_CTX(ctx))) ? T_TRUE : T_FALSE;
}

static int
_gr_nmod_redc_set(ulong * res, const ulong * x, gr_ctx_t FLINT_UNUSED(ctx))
{
    res[0] = x[0];
    return GR_SUCCESS;
}

static int
_gr_nmod_redc_neg(ulong * res, const ulong * x, gr_ctx_t ctx)
{
    res[0] = nmod_redc_neg(x[0], GR_NMOD_REDC_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_nmod_redc_fast_neg(ulong * res, const ulong * x, gr_ctx_t ctx)
{
    res[0] = nmod_redc_fast_neg(x[0], GR_NMOD_REDC_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_nmod_redc_add(ulong * res, const ulong * x, const ulong * y, gr_ctx_t ctx)
{
    res[0] = nmod_redc_add(x[0], y[0], GR_NMOD_REDC_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_nmod_redc_fast_add(ulong * res, const ulong * x, const ulong * y, gr_ctx_t ctx)
{
    res[0] = nmod_redc_fast_add(x[0], y[0], GR_NMOD_REDC_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_nmod_redc_sub(ulong * res, const ulong * x, const ulong * y, gr_ctx_t ctx)
{
    res[0] = nmod_redc_sub(x[0], y[0], GR_NMOD_REDC_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_nmod_redc_fast_sub(ulong * res, const ulong * x, const ulong * y, gr_ctx_t ctx)
{
    res[0] = nmod_redc_fast_sub(x[0], y[0], GR_NMOD_REDC_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_nmod_redc_mul(ulong * res, const ulong * x, const ulong * y, gr_ctx_t ctx)
{
    res[0] = nmod_redc_mul(x[0], y[0], GR_NMOD_REDC_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_nmod_redc_fast_mul(ulong * res, const ulong * x, const ulong * y, gr_ctx_t ctx)
{
    res[0] = nmod_redc_fast_mul(x[0], y[0], GR_NMOD_REDC_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_nmod_redc_div(ulong * res, const ulong * x, const ulong * y, gr_ctx_t ctx)
{
    ulong t;
    int status;

    status = _gr_nmod_redc_inv(&t, y, ctx);
    if (status != GR_SUCCESS)
        return status;

    res[0] = nmod_redc_mul(x[0], t, GR_NMOD_REDC_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_nmod_redc_fast_div(ulong * res, const ulong * x, const ulong * y, gr_ctx_t ctx)
{
    ulong t;
    int status;

    status = _gr_nmod_redc_inv(&t, y, ctx);
    if (status != GR_SUCCESS)
        return status;

    res[0] = nmod_redc_fast_mul(x[0], t, GR_NMOD_REDC_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_nmod_redc_sqr(ulong * res, const ulong * x, gr_ctx_t ctx)
{
    res[0] = nmod_redc_mul(x[0], x[0], GR_NMOD_REDC_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_nmod_redc_fast_sqr(ulong * res, const ulong * x, gr_ctx_t ctx)
{
    res[0] = nmod_redc_fast_mul(x[0], x[0], GR_NMOD_REDC_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_nmod_redc_pow_ui(ulong * res, const ulong * x, ulong e, gr_ctx_t ctx)
{
    if (e == 0)
        res[0] = GR_NMOD_REDC_ONE(ctx);
    else
        res[0] = _nmod_redc_pow_ui(*x, e, GR_NMOD_REDC_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_nmod_redc_fast_pow_ui(ulong * res, const ulong * x, ulong e, gr_ctx_t ctx)
{
    if (e == 0)
        res[0] = GR_NMOD_REDC_ONE(ctx);
    else
        res[0] = _nmod_redc_fast_pow_ui(*x, e, GR_NMOD_REDC_CTX(ctx));
    return GR_SUCCESS;
}

static void
_gr_nmod_redc_vec_init(ulong * res, slong len, gr_ctx_t FLINT_UNUSED(ctx))
{
    slong i;

    for (i = 0; i < len; i++)
        res[i] = 0;
}

static void
_gr_nmod_redc_vec_clear(ulong * FLINT_UNUSED(res), slong FLINT_UNUSED(len), gr_ctx_t FLINT_UNUSED(ctx))
{
}

static int
_gr_nmod_redc_vec_set(ulong * res, const ulong * vec, slong len, gr_ctx_t FLINT_UNUSED(ctx))
{
    slong i;

    for (i = 0; i < len; i++)
        res[i] = vec[i];

    return GR_SUCCESS;
}

static int
_gr_nmod_redc_vec_normalise(slong * res, const ulong * vec, slong len, gr_ctx_t FLINT_UNUSED(ctx))
{
    while (len > 0 && vec[len - 1] == 0)
        len--;

    res[0] = len;
    return GR_SUCCESS;
}

static int
_gr_nmod_redc_fast_vec_normalise(slong * res, const ulong * vec, slong len, gr_ctx_t ctx)
{
    ulong n = GR_NMOD_REDC_N(ctx);

    while (len > 0 && (vec[len - 1] == 0 || vec[len - 1] == n))
        len--;

    res[0] = len;
    return GR_SUCCESS;
}

static slong
_gr_nmod_redc_vec_normalise_weak(const ulong * vec, slong len, gr_ctx_t FLINT_UNUSED(ctx))
{
    while (len > 0 && vec[len - 1] == 0)
        len--;

    return len;
}

static slong
_gr_nmod_redc_fast_vec_normalise_weak(const ulong * vec, slong len, gr_ctx_t ctx)
{
    ulong n = GR_NMOD_REDC_N(ctx);

    while (len > 0 && (vec[len - 1] == 0 || vec[len - 1] == n))
        len--;

    return len;
}

static int
_gr_nmod_redc_vec_neg(ulong * res, const ulong * vec, slong len, gr_ctx_t ctx)
{
    slong i;
    ulong n = GR_NMOD_REDC_N(ctx);

    for (i = 0; i < len; i++)
        res[i] = n_negmod(vec[i], n);

    return GR_SUCCESS;
}


static int
_gr_nmod_redc_fast_vec_neg(ulong * res, const ulong * vec, slong len, gr_ctx_t ctx)
{
    slong i;
    ulong n = GR_NMOD_REDC_N(ctx);

    for (i = 0; i < len; i++)
        res[i] = n_negmod(vec[i], 2 * n);

    return GR_SUCCESS;
}

static
ulong _n_mod_add(ulong a, ulong b, ulong n)
{
   const ulong sum = a + b;
   return sum - n + ((((slong)(sum - n))>>(FLINT_BITS - 1)) & n);
}

static
ulong _n_mod_sub(ulong a, ulong b, ulong n)
{
   const ulong diff = a - b;
   return  ((((slong)diff)>>(FLINT_BITS - 1)) & n) + diff;
}


static int
_gr_nmod_redc_vec_add(ulong * res, const ulong * vec1, const ulong * vec2, slong len, gr_ctx_t ctx)
{
    slong i;
    ulong n = GR_NMOD_REDC_N(ctx);

    if (GR_NMOD_REDC_MOD(ctx).norm)
        for (i = 0; i < len; i++)
            res[i] = _n_mod_add(vec1[i], vec2[i], n);
    else
        for (i = 0; i < len; i++)
            res[i] = n_addmod(vec1[i], vec2[i], n);

    return GR_SUCCESS;
}

static int
_gr_nmod_redc_fast_vec_add(ulong * res, const ulong * vec1, const ulong * vec2, slong len, gr_ctx_t ctx)
{
    slong i;
    ulong n = GR_NMOD_REDC_N(ctx);

    for (i = 0; i < len; i++)
        res[i] = _n_mod_add(vec1[i], vec2[i], 2 * n);

    return GR_SUCCESS;
}

static int
_gr_nmod_redc_vec_sub(ulong * res, const ulong * vec1, const ulong * vec2, slong len, gr_ctx_t ctx)
{
    slong i;
    ulong n = GR_NMOD_REDC_N(ctx);

    if (GR_NMOD_REDC_MOD(ctx).norm)
        for (i = 0; i < len; i++)
            res[i] = _n_mod_sub(vec1[i], vec2[i], n);
    else
        for (i = 0; i < len; i++)
            res[i] = n_submod(vec1[i], vec2[i], n);

    return GR_SUCCESS;
}

static int
_gr_nmod_redc_fast_vec_sub(ulong * res, const ulong * vec1, const ulong * vec2, slong len, gr_ctx_t ctx)
{
    slong i;
    ulong n = GR_NMOD_REDC_N(ctx);

    for (i = 0; i < len; i++)
        res[i] = _n_mod_sub(vec1[i], vec2[i], 2 * n);

    return GR_SUCCESS;
}

static int
_gr_nmod_redc_vec_mul(ulong * res, const ulong * vec1, const ulong * vec2, slong len, gr_ctx_t ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        res[i] = nmod_redc_mul(vec1[i], vec2[i], GR_NMOD_REDC_CTX(ctx));

    return GR_SUCCESS;
}

static int
_gr_nmod_redc_fast_vec_mul(ulong * res, const ulong * vec1, const ulong * vec2, slong len, gr_ctx_t ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        res[i] = nmod_redc_fast_mul(vec1[i], vec2[i], GR_NMOD_REDC_CTX(ctx));

    return GR_SUCCESS;
}

static int
_gr_nmod_redc_vec_mul_scalar(ulong * res, const ulong * vec1, slong len, const ulong * c, gr_ctx_t ctx)
{
    slong i;
    ulong d = *c;

    for (i = 0; i < len; i++)
        res[i] = nmod_redc_mul(vec1[i], d, GR_NMOD_REDC_CTX(ctx));

    return GR_SUCCESS;
}

static int
_gr_nmod_redc_fast_vec_mul_scalar(ulong * res, const ulong * vec1, slong len, const ulong * c, gr_ctx_t ctx)
{
    slong i;
    ulong d = *c;

    for (i = 0; i < len; i++)
        res[i] = nmod_redc_fast_mul(vec1[i], d, GR_NMOD_REDC_CTX(ctx));

    return GR_SUCCESS;
}

static int
_gr_nmod_redc_scalar_mul_vec(ulong * res, ulong * c, const ulong * vec1, slong len, gr_ctx_t ctx)
{
    return _gr_nmod_redc_vec_mul_scalar(res, vec1, len, c, ctx);
}

static int
_gr_nmod_redc_fast_scalar_mul_vec(ulong * res, ulong * c, const ulong * vec1, slong len, gr_ctx_t ctx)
{
    return _gr_nmod_redc_fast_vec_mul_scalar(res, vec1, len, c, ctx);
}

static int
_gr_nmod_redc_vec_mul_scalar_si(ulong * res, const ulong * vec1, slong len, slong c, gr_ctx_t ctx)
{
    ulong d = nmod_redc_set_si(c, GR_NMOD_REDC_CTX(ctx));
    return _gr_nmod_redc_vec_mul_scalar(res, vec1, len, &d, ctx);
}

static int
_gr_nmod_redc_fast_vec_mul_scalar_si(ulong * res, const ulong * vec1, slong len, slong c, gr_ctx_t ctx)
{
    ulong d = nmod_redc_set_si(c, GR_NMOD_REDC_CTX(ctx));
    return _gr_nmod_redc_fast_vec_mul_scalar(res, vec1, len, &d, ctx);
}

static int
_gr_nmod_redc_vec_mul_scalar_ui(ulong * res, const ulong * vec1, slong len, ulong c, gr_ctx_t ctx)
{
    ulong d = nmod_redc_set_ui(c, GR_NMOD_REDC_CTX(ctx));
    return _gr_nmod_redc_vec_mul_scalar(res, vec1, len, &d, ctx);
}

static int
_gr_nmod_redc_fast_vec_mul_scalar_ui(ulong * res, const ulong * vec1, slong len, ulong c, gr_ctx_t ctx)
{
    ulong d = nmod_redc_set_ui(c, GR_NMOD_REDC_CTX(ctx));
    return _gr_nmod_redc_fast_vec_mul_scalar(res, vec1, len, &d, ctx);
}

static int
_gr_nmod_redc_vec_mul_scalar_fmpz(ulong * res, const ulong * vec1, slong len, const fmpz_t c, gr_ctx_t ctx)
{
    ulong d = nmod_redc_set_fmpz(c, GR_NMOD_REDC_CTX(ctx));
    return _gr_nmod_redc_vec_mul_scalar(res, vec1, len, &d, ctx);
}

static int
_gr_nmod_redc_fast_vec_mul_scalar_fmpz(ulong * res, const ulong * vec1, slong len, const fmpz_t c, gr_ctx_t ctx)
{
    ulong d = nmod_redc_set_fmpz(c, GR_NMOD_REDC_CTX(ctx));
    return _gr_nmod_redc_fast_vec_mul_scalar(res, vec1, len, &d, ctx);
}

static int
_gr_nmod_redc_vec_addmul_scalar(ulong * res, const ulong * vec1, slong len, const ulong * c, gr_ctx_t ctx)
{
    slong i;
    ulong n = GR_NMOD_REDC_N(ctx);
    ulong d = *c;

    if (GR_NMOD_REDC_MOD(ctx).norm)
        for (i = 0; i < len; i++)
            res[i] = _n_mod_add(res[i], nmod_redc_mul(vec1[i], d, GR_NMOD_REDC_CTX(ctx)), n);
    else
        for (i = 0; i < len; i++)
            res[i] = n_addmod(res[i], nmod_redc_mul(vec1[i], d, GR_NMOD_REDC_CTX(ctx)), n);

    return GR_SUCCESS;
}

static int
_gr_nmod_redc_vec_submul_scalar(ulong * res, const ulong * vec1, slong len, const ulong * c, gr_ctx_t ctx)
{
    ulong d = nmod_redc_neg(*c, GR_NMOD_REDC_CTX(ctx));
    return _gr_nmod_redc_vec_addmul_scalar(res, vec1, len, &d, ctx);
}

static int
_gr_nmod_redc_fast_vec_addmul_scalar(ulong * res, const ulong * vec1, slong len, const ulong * c, gr_ctx_t ctx)
{
    slong i;
    ulong n = GR_NMOD_REDC_N(ctx);
    ulong d = *c;

    for (i = 0; i < len; i++)
        res[i] = _n_mod_add(res[i], nmod_redc_fast_mul(vec1[i], d, GR_NMOD_REDC_CTX(ctx)), 2 * n);

    return GR_SUCCESS;
}

static int
_gr_nmod_redc_fast_vec_submul_scalar(ulong * res, const ulong * vec1, slong len, const ulong * c, gr_ctx_t ctx)
{
    ulong d = nmod_redc_fast_neg(*c, GR_NMOD_REDC_CTX(ctx));
    return _gr_nmod_redc_fast_vec_addmul_scalar(res, vec1, len, &d, ctx);
}

static int
_gr_nmod_redc_vec_addmul_scalar_si(ulong * res, const ulong * vec1, slong len, slong c, gr_ctx_t ctx)
{
    ulong d = nmod_redc_set_si(c, GR_NMOD_REDC_CTX(ctx));
    return _gr_nmod_redc_vec_addmul_scalar(res, vec1, len, &d, ctx);
}

static int
_gr_nmod_redc_fast_vec_addmul_scalar_si(ulong * res, const ulong * vec1, slong len, slong c, gr_ctx_t ctx)
{
    ulong d = nmod_redc_set_si(c, GR_NMOD_REDC_CTX(ctx));
    return _gr_nmod_redc_fast_vec_addmul_scalar(res, vec1, len, &d, ctx);
}

static int
_gr_nmod_redc_vec_submul_scalar_si(ulong * res, const ulong * vec1, slong len, slong c, gr_ctx_t ctx)
{
    ulong d = nmod_redc_set_si(c, GR_NMOD_REDC_CTX(ctx));
    return _gr_nmod_redc_vec_submul_scalar(res, vec1, len, &d, ctx);
}

static int
_gr_nmod_redc_fast_vec_submul_scalar_si(ulong * res, const ulong * vec1, slong len, slong c, gr_ctx_t ctx)
{
    ulong d = nmod_redc_set_si(c, GR_NMOD_REDC_CTX(ctx));
    return _gr_nmod_redc_fast_vec_submul_scalar(res, vec1, len, &d, ctx);
}

static int
_gr_nmod_redc_vec_product(ulong * res, const ulong * vec, slong len, gr_ctx_t ctx)
{
    if (len <= 1)
    {
        if (len == 1)
            res[0] = vec[0];
        else
            res[0] = GR_NMOD_REDC_ONE(ctx);
    }
    else
    {
        ulong p;
        slong i;
        p = nmod_redc_mul(vec[0], vec[1], GR_NMOD_REDC_CTX(ctx));
        for (i = 2; i < len; i++)
            p = nmod_redc_mul(p, vec[i], GR_NMOD_REDC_CTX(ctx));
        res[0] = p;
    }

    return GR_SUCCESS;
}

static int
_gr_nmod_redc_fast_vec_product(ulong * res, const ulong * vec, slong len, gr_ctx_t ctx)
{
    if (len <= 1)
    {
        if (len == 1)
            res[0] = vec[0];
        else
            res[0] = GR_NMOD_REDC_ONE(ctx);
    }
    else
    {
        ulong p;
        slong i;
        p = nmod_redc_fast_mul(vec[0], vec[1], GR_NMOD_REDC_CTX(ctx));
        for (i = 2; i < len; i++)
            p = nmod_redc_fast_mul(p, vec[i], GR_NMOD_REDC_CTX(ctx));
        res[0] = p;
    }

    return GR_SUCCESS;
}


/* If len is short, try to benefit from fast redc modular reduction.
   Otherwise, fall back on _nmod_vec_dot which might do something more
   clever. The cutoff hasn't been carefully tuned.

   To do: combine the _nmod_vec_dot and redc techniques. */
#define NMOD_REDC_DOT_CUTOFF 40

#define NMOD_REDC_DOT_1(s, i, vec1_i, vec2_i, len, ctx) \
    do \
    { \
        ull_t t2, s2; \
        i = 0; \
        s2 = ull_u_mul_u(vec1_i, vec2_i); \
        for (i = 1; i < len; i++) \
        { \
            t2 = ull_u_mul_u(vec1_i, vec2_i); \
            s2 = ull_add(s2, t2); \
        } \
        (s) = n_ll_redc(s2, GR_NMOD_REDC_N(ctx), GR_NMOD_REDC_NRED(ctx)); \
    } while (0)

#define NMOD_REDC_DOT_2(s, i, vec1_i, vec2_i, len, ctx) \
    do \
    { \
        ull_t t2, s2; \
        i = 0; \
        s2 = ull_u_mul_u(vec1_i, vec2_i); \
        for (i = 1; i < len; i++) \
        { \
            t2 = ull_u_mul_u(vec1_i, vec2_i); \
            s2 = ull_add(s2, t2); \
        } \
        ulong hi, lo; \
        lo = ull_lo(s2); \
        hi = ull_hi(s2); \
        NMOD_RED2(hi, UWORD(0), hi, GR_NMOD_REDC_MOD(ctx)); \
        s2 = ull(hi, lo); \
        (s) = n_ll_redc(s2, GR_NMOD_REDC_N(ctx), GR_NMOD_REDC_NRED(ctx)); \
    } while (0)

#define NMOD_REDC_DOT_3(s, i, vec1_i, vec2_i, len, ctx) \
    do \
    { \
        ulong s2, s1, s0, t1, t0; \
        ull_t s3; \
        i = 0; \
        umul_ppmm(s1, s0, vec1_i, vec2_i); \
        s2 = 0; \
        for (i = 1; i < len; i++) \
        { \
            umul_ppmm(t1, t0, vec1_i, vec2_i); \
            add_sssaaaaaa(s2, s1, s0, s2, s1, s0, 0, t1, t0); \
        } \
        FLINT_ASSERT(s2 < GR_NMOD_REDC_N(ctx)); \
        NMOD_RED2(s1, s2, s1, GR_NMOD_REDC_MOD(ctx)); \
        s3 = ull(s1, s0); \
        s = n_ll_redc(s3, GR_NMOD_REDC_N(ctx), GR_NMOD_REDC_NRED(ctx)); \
    } while (0)

#define NMOD_REDC_FAST_DOT_1(s, i, vec1_i, vec2_i, len, ctx) \
    do \
    { \
        ull_t t2, s2; \
        i = 0; \
        s2 = ull_u_mul_u(vec1_i, vec2_i); \
        for (i = 1; i < len; i++) \
        { \
            t2 = ull_u_mul_u(vec1_i, vec2_i); \
            s2 = ull_add(s2, t2); \
        } \
        (s) = n_ll_redc_fast(s2, GR_NMOD_REDC_N(ctx), GR_NMOD_REDC_NRED(ctx)); \
    } while (0)

#define NMOD_REDC_FAST_DOT_2(s, i, vec1_i, vec2_i, len, ctx) \
    do \
    { \
        ull_t t2, s2; \
        i = 0; \
        s2 = ull_u_mul_u(vec1_i, vec2_i); \
        for (i = 1; i < len; i++) \
        { \
            t2 = ull_u_mul_u(vec1_i, vec2_i); \
            s2 = ull_add(s2, t2); \
        } \
        ulong hi, lo; \
        lo = ull_lo(s2); \
        hi = ull_hi(s2); \
        NMOD_RED2_NONFULLWORD(hi, UWORD(0), hi, GR_NMOD_REDC_MOD(ctx)); \
        s2 = ull(hi, lo); \
        (s) = n_ll_redc_fast(s2, GR_NMOD_REDC_N(ctx), GR_NMOD_REDC_NRED(ctx)); \
    } while (0)

#define NMOD_REDC_FAST_DOT_3(s, i, vec1_i, vec2_i, len, ctx) \
    do \
    { \
        ulong s2, s1, s0, t1, t0; \
        ull_t s3; \
        i = 0; \
        umul_ppmm(s1, s0, vec1_i, vec2_i); \
        s2 = 0; \
        for (i = 1; i < len; i++) \
        { \
            umul_ppmm(t1, t0, vec1_i, vec2_i); \
            add_sssaaaaaa(s2, s1, s0, s2, s1, s0, 0, t1, t0); \
        } \
        FLINT_ASSERT(s2 < GR_NMOD_REDC_N(ctx)); \
        NMOD_RED2_NONFULLWORD(s1, s2, s1, GR_NMOD_REDC_MOD(ctx)); \
        s3 = ull(s1, s0); \
        s = n_ll_redc_fast(s3, GR_NMOD_REDC_N(ctx), GR_NMOD_REDC_NRED(ctx)); \
    } while (0)

static int
_gr_nmod_redc_vec_dot(ulong * res, const ulong * initial, int subtract, const ulong * vec1, const ulong * vec2, slong len, gr_ctx_t ctx)
{
    ulong s;
    dot_params_t params;
    nmod_t mod;
    slong i;

    if (len == 0)
    {
        if (initial == NULL)
            *res = 0;
        else
            *res = *initial;
        return GR_SUCCESS;
    }

    if (len < NMOD_REDC_DOT_CUTOFF)
    {
        ulong nbits, sbits;

        nbits = NMOD_BITS(GR_NMOD_REDC_MOD(ctx));
        /* sum <= len * (n-1)^2 */
        sbits = 2 * nbits + FLINT_BIT_COUNT(len);

        if (sbits < FLINT_BITS + nbits)
        {
            NMOD_REDC_DOT_1(s, i, vec1[i], vec2[i], len, ctx);
        }
        else if (sbits <= 2 * FLINT_BITS)
        {
            NMOD_REDC_DOT_2(s, i, vec1[i], vec2[i], len, ctx);
        }
        else
        {
            NMOD_REDC_DOT_3(s, i, vec1[i], vec2[i], len, ctx);
        }
    }
    else
    {
        mod = GR_NMOD_REDC_MOD(ctx);
        params = _nmod_vec_dot_params(len, mod);
        s = _nmod_vec_dot(vec1, vec2, len, mod, params);
        s = nmod_redc_get_nmod(s, GR_NMOD_REDC_CTX(ctx));
    }

    if (initial == NULL)
    {
        mod = GR_NMOD_REDC_MOD(ctx);

        if (subtract)
            s = nmod_neg(s, mod);
    }
    else
    {
        mod = GR_NMOD_REDC_MOD(ctx);

        if (subtract)
            s = nmod_sub(initial[0], s, mod);
        else
            s = nmod_add(initial[0], s, mod);
    }

    *res = s;

    return GR_SUCCESS;
}

static int
_gr_nmod_redc_vec_dot_rev(ulong * res, const ulong * initial, int subtract, const ulong * vec1, const ulong * vec2, slong len, gr_ctx_t ctx)
{
    ulong s;
    dot_params_t params;
    nmod_t mod;
    slong i;

    if (len == 0)
    {
        if (initial == NULL)
            *res = 0;
        else
            *res = *initial;
        return GR_SUCCESS;
    }

    if (len < NMOD_REDC_DOT_CUTOFF)
    {
        ulong nbits, sbits;

        nbits = NMOD_BITS(GR_NMOD_REDC_MOD(ctx));
        /* sum <= len * (n-1)^2 */
        sbits = 2 * nbits + FLINT_BIT_COUNT(len);

        if (sbits < FLINT_BITS + nbits)
        {
            NMOD_REDC_DOT_1(s, i, vec1[i], vec2[len - 1 - i], len, ctx);
        }
        else if (sbits <= 2 * FLINT_BITS)
        {
            NMOD_REDC_DOT_2(s, i, vec1[i], vec2[len - 1 - i], len, ctx);
        }
        else
        {
            NMOD_REDC_DOT_3(s, i, vec1[i], vec2[len - 1 - i], len, ctx);
        }
    }
    else
    {
        mod = GR_NMOD_REDC_MOD(ctx);
        params = _nmod_vec_dot_params(len, mod);
        s = _nmod_vec_dot_rev(vec1, vec2, len, mod, params);
        s = nmod_redc_get_nmod(s, GR_NMOD_REDC_CTX(ctx));
    }

    if (initial == NULL)
    {
        mod = GR_NMOD_REDC_MOD(ctx);

        if (subtract)
            s = nmod_neg(s, mod);
    }
    else
    {
        mod = GR_NMOD_REDC_MOD(ctx);

        if (subtract)
            s = nmod_sub(initial[0], s, mod);
        else
            s = nmod_add(initial[0], s, mod);
    }

    *res = s;

    return GR_SUCCESS;
}

static int
_gr_nmod_redc_fast_vec_dot(ulong * res, const ulong * initial, int subtract, const ulong * vec1, const ulong * vec2, slong len, gr_ctx_t ctx)
{
    ulong s;
    dot_params_t params;
    nmod_t mod;
    slong i;

    if (len == 0)
    {
        if (initial == NULL)
            *res = 0;
        else
            *res = *initial;
        return GR_SUCCESS;
    }

    if (len < NMOD_REDC_DOT_CUTOFF)
    {
        ulong nbits, sbits;

        nbits = NMOD_BITS(GR_NMOD_REDC_MOD(ctx));
        /* sum <= len * (2n-1)^2 */
        sbits = 2 * (nbits + 1) + FLINT_BIT_COUNT(len);

        if (sbits < FLINT_BITS + nbits)
        {
            NMOD_REDC_FAST_DOT_1(s, i, vec1[i], vec2[i], len, ctx);
        }
        else if (sbits <= 2 * FLINT_BITS)
        {
            NMOD_REDC_FAST_DOT_2(s, i, vec1[i], vec2[i], len, ctx);
        }
        else
        {
            NMOD_REDC_FAST_DOT_3(s, i, vec1[i], vec2[i], len, ctx);
        }

        if (initial != NULL || subtract)
        {
            mod = GR_NMOD_REDC_MOD(ctx);
            mod.n = mod.n * 2;
            mod.ninv = mod.ninv;
            mod.norm = mod.norm - 1;

            if (initial == NULL)
            {
                if (subtract)
                    s = nmod_neg(s, mod);
            }
            else
            {
                if (subtract)
                    s = nmod_sub(initial[0], s, mod);
                else
                    s = nmod_add(initial[0], s, mod);
            }
        }
    }
    else
    {
        mod = GR_NMOD_REDC_MOD(ctx);
        mod.n = mod.n * 2;
        mod.ninv = mod.ninv;
        mod.norm = mod.norm - 1;

        params = _nmod_vec_dot_params(len, mod);
        s = _nmod_vec_dot(vec1, vec2, len, mod, params);
        s = nmod_redc_get_nmod(s, GR_NMOD_REDC_CTX(ctx));

        if (initial == NULL)
        {
            if (subtract)
                s = nmod_neg(s, mod);
        }
        else
        {
            if (subtract)
                s = nmod_sub(initial[0], s, mod);
            else
                s = nmod_add(initial[0], s, mod);
        }
    }

    *res = s;

    return GR_SUCCESS;
}

static int
_gr_nmod_redc_fast_vec_dot_rev(ulong * res, const ulong * initial, int subtract, const ulong * vec1, const ulong * vec2, slong len, gr_ctx_t ctx)
{
    ulong s;
    dot_params_t params;
    nmod_t mod;
    slong i;

    if (len == 0)
    {
        if (initial == NULL)
            *res = 0;
        else
            *res = *initial;
        return GR_SUCCESS;
    }

    if (len < NMOD_REDC_DOT_CUTOFF)
    {
        ulong nbits, sbits;

        nbits = NMOD_BITS(GR_NMOD_REDC_MOD(ctx));
        /* sum <= len * (2n-1)^2 */
        sbits = 2 * (nbits + 1) + FLINT_BIT_COUNT(len);

        if (sbits < FLINT_BITS + nbits)
        {
            NMOD_REDC_FAST_DOT_1(s, i, vec1[i], vec2[len - 1 - i], len, ctx);
        }
        else if (sbits <= 2 * FLINT_BITS)
        {
            NMOD_REDC_FAST_DOT_2(s, i, vec1[i], vec2[len - 1 - i], len, ctx);
        }
        else
        {
            NMOD_REDC_FAST_DOT_3(s, i, vec1[i], vec2[len - 1 - i], len, ctx);
        }

        if (initial != NULL || subtract)
        {
            mod = GR_NMOD_REDC_MOD(ctx);
            mod.n = mod.n * 2;
            mod.ninv = mod.ninv;
            mod.norm = mod.norm - 1;

            if (initial == NULL)
            {
                if (subtract)
                    s = nmod_neg(s, mod);
            }
            else
            {
                if (subtract)
                    s = nmod_sub(initial[0], s, mod);
                else
                    s = nmod_add(initial[0], s, mod);
            }
        }
    }
    else
    {
        mod = GR_NMOD_REDC_MOD(ctx);
        mod.n = mod.n * 2;
        mod.ninv = mod.ninv;
        mod.norm = mod.norm - 1;

        params = _nmod_vec_dot_params(len, mod);
        s = _nmod_vec_dot_rev(vec1, vec2, len, mod, params);
        s = nmod_redc_get_nmod(s, GR_NMOD_REDC_CTX(ctx));

        if (initial == NULL)
        {
            if (subtract)
                s = nmod_neg(s, mod);
        }
        else
        {
            if (subtract)
                s = nmod_sub(initial[0], s, mod);
            else
                s = nmod_add(initial[0], s, mod);
        }
    }

    *res = s;

    return GR_SUCCESS;
}

static int
_gr_nmod_redc_poly_mullow(ulong * res,
    const ulong * poly1, slong len1,
    const ulong * poly2, slong len2, slong n, gr_ctx_t ctx)
{
    slong i, j, m, n1, n2;

    len1 = FLINT_MIN(len1, n);
    len2 = FLINT_MIN(len2, n);
    m = FLINT_MIN(len1, len2);

    if (m < NMOD_REDC_DOT_CUTOFF)
    {
        ulong nbits, sbits;

        nbits = NMOD_BITS(GR_NMOD_REDC_MOD(ctx));
        sbits = 2 * (nbits + 1) + FLINT_BIT_COUNT(m);

        if (poly1 == poly2 && len1 == len2)
        {
            ulong cc;

            res[0] = nmod_redc_mul(poly1[0], poly1[0], GR_NMOD_REDC_CTX(ctx));

            /* todo: avoid the extra reductions */

            if (sbits < FLINT_BITS + nbits)
            {
                for (i = 1; i < FLINT_MIN(n, 2 * len1 - 2); i++)
                {
                    n1 = FLINT_MAX(0, i - len1 + 1);
                    n2 = FLINT_MIN(len1 - 1, (i + 1) / 2 - 1);

                    NMOD_REDC_DOT_1(cc, j, poly1[n1 + j], poly1[i - n2 + (n2 - n1 + 1) - 1 - j], n2 - n1 + 1, ctx);
                    cc = nmod_redc_add(cc, cc, GR_NMOD_REDC_CTX(ctx));
                    if (i % 2 == 0 && i / 2 < len1)
                        cc = nmod_redc_add(cc, nmod_redc_mul(poly1[i / 2], poly1[i / 2], GR_NMOD_REDC_CTX(ctx)), GR_NMOD_REDC_CTX(ctx));

                    res[i] = cc;
                }
            }
            else if (sbits <= 2 * FLINT_BITS)
            {
                for (i = 1; i < FLINT_MIN(n, 2 * len1 - 2); i++)
                {
                    n1 = FLINT_MAX(0, i - len1 + 1);
                    n2 = FLINT_MIN(len1 - 1, (i + 1) / 2 - 1);

                    NMOD_REDC_DOT_2(cc, j, poly1[n1 + j], poly1[i - n2 + (n2 - n1 + 1) - 1 - j], n2 - n1 + 1, ctx);
                    cc = nmod_redc_add(cc, cc, GR_NMOD_REDC_CTX(ctx));
                    if (i % 2 == 0 && i / 2 < len1)
                        cc = nmod_redc_add(cc, nmod_redc_mul(poly1[i / 2], poly1[i / 2], GR_NMOD_REDC_CTX(ctx)), GR_NMOD_REDC_CTX(ctx));

                    res[i] = cc;
                }
            }
            else
            {
                for (i = 1; i < FLINT_MIN(n, 2 * len1 - 2); i++)
                {
                    n1 = FLINT_MAX(0, i - len1 + 1);
                    n2 = FLINT_MIN(len1 - 1, (i + 1) / 2 - 1);

                    NMOD_REDC_DOT_3(cc, j, poly1[n1 + j], poly1[i - n2 + (n2 - n1 + 1) - 1 - j], n2 - n1 + 1, ctx);
                    cc = nmod_redc_add(cc, cc, GR_NMOD_REDC_CTX(ctx));
                    if (i % 2 == 0 && i / 2 < len1)
                        cc = nmod_redc_add(cc, nmod_redc_mul(poly1[i / 2], poly1[i / 2], GR_NMOD_REDC_CTX(ctx)), GR_NMOD_REDC_CTX(ctx));

                    res[i] = cc;
                }
            }

            if (n >= 2 * len1 - 1)
                res[2 * len1 - 2] = nmod_redc_mul(poly1[len1 - 1], poly1[len1 - 1], GR_NMOD_REDC_CTX(ctx));
        }
        else
        {
            if (sbits < FLINT_BITS + nbits)
            {
                for (i = 0; i < n; i++)
                {
                    n1 = FLINT_MIN(len1 - 1, i);
                    n2 = FLINT_MIN(len2 - 1, i);
                    NMOD_REDC_DOT_1(res[i], j, poly1[i - n2 + j], poly2[i + n2 - i - j], n1 + n2 - i + 1, ctx);
                }
            }
            else if (sbits <= 2 * FLINT_BITS)
            {
                for (i = 0; i < n; i++)
                {
                    n1 = FLINT_MIN(len1 - 1, i);
                    n2 = FLINT_MIN(len2 - 1, i);
                    NMOD_REDC_DOT_2(res[i], j, poly1[i - n2 + j], poly2[i + n2 - i - j], n1 + n2 - i + 1, ctx);
                }
            }
            else
            {
                for (i = 0; i < n; i++)
                {
                    n1 = FLINT_MIN(len1 - 1, i);
                    n2 = FLINT_MIN(len2 - 1, i);
                    NMOD_REDC_DOT_3(res[i], j, poly1[i - n2 + j], poly2[i + n2 - i - j], n1 + n2 - i + 1, ctx);
                }
            }
        }

        return GR_SUCCESS;
    }

    if (len1 + len2 - 1 == n)
    {
        if (len1 >= len2)
            _nmod_poly_mul(res, poly1, len1, poly2, len2, GR_NMOD_REDC_MOD(ctx));
        else
            _nmod_poly_mul(res, poly2, len2, poly1, len1, GR_NMOD_REDC_MOD(ctx));
    }
    else
    {
        if (len1 >= len2)
            _nmod_poly_mullow(res, poly1, len1, poly2, len2, n, GR_NMOD_REDC_MOD(ctx));
        else
            _nmod_poly_mullow(res, poly2, len2, poly1, len1, n, GR_NMOD_REDC_MOD(ctx));
    }

    for (i = 0; i < n; i++)
        res[i] = nmod_redc_get_nmod(res[i], GR_NMOD_REDC_CTX(ctx));

    return GR_SUCCESS;
}

static int
_gr_nmod_redc_fast_poly_mullow(ulong * res,
    const ulong * poly1, slong len1,
    const ulong * poly2, slong len2, slong n, gr_ctx_t ctx)
{
    nn_ptr t1, t2;
    slong i, j, m, n1, n2, alloc;

    len1 = FLINT_MIN(len1, n);
    len2 = FLINT_MIN(len2, n);
    m = FLINT_MIN(len1, len2);

    if (m < NMOD_REDC_DOT_CUTOFF)
    {
        ulong nbits, sbits;

        nbits = NMOD_BITS(GR_NMOD_REDC_MOD(ctx));
        sbits = 2 * (nbits + 1) + FLINT_BIT_COUNT(m);

        if (poly1 == poly2 && len1 == len2)
        {
            ulong cc;

            res[0] = nmod_redc_fast_mul(poly1[0], poly1[0], GR_NMOD_REDC_CTX(ctx));

            /* todo: avoid the extra reductions */

            if (sbits < FLINT_BITS + nbits)
            {
                for (i = 1; i < FLINT_MIN(n, 2 * len1 - 2); i++)
                {
                    n1 = FLINT_MAX(0, i - len1 + 1);
                    n2 = FLINT_MIN(len1 - 1, (i + 1) / 2 - 1);

                    NMOD_REDC_FAST_DOT_1(cc, j, poly1[n1 + j], poly1[i - n2 + (n2 - n1 + 1) - 1 - j], n2 - n1 + 1, ctx);
                    cc = nmod_redc_fast_add(cc, cc, GR_NMOD_REDC_CTX(ctx));
                    if (i % 2 == 0 && i / 2 < len1)
                        cc = nmod_redc_fast_add(cc, nmod_redc_fast_mul(poly1[i / 2], poly1[i / 2], GR_NMOD_REDC_CTX(ctx)), GR_NMOD_REDC_CTX(ctx));

                    res[i] = cc;
                }
            }
            else if (sbits <= 2 * FLINT_BITS)
            {
                for (i = 1; i < FLINT_MIN(n, 2 * len1 - 2); i++)
                {
                    n1 = FLINT_MAX(0, i - len1 + 1);
                    n2 = FLINT_MIN(len1 - 1, (i + 1) / 2 - 1);

                    NMOD_REDC_FAST_DOT_2(cc, j, poly1[n1 + j], poly1[i - n2 + (n2 - n1 + 1) - 1 - j], n2 - n1 + 1, ctx);
                    cc = nmod_redc_fast_add(cc, cc, GR_NMOD_REDC_CTX(ctx));
                    if (i % 2 == 0 && i / 2 < len1)
                        cc = nmod_redc_fast_add(cc, nmod_redc_mul(poly1[i / 2], poly1[i / 2], GR_NMOD_REDC_CTX(ctx)), GR_NMOD_REDC_CTX(ctx));

                    res[i] = cc;
                }
            }
            else
            {
                for (i = 1; i < FLINT_MIN(n, 2 * len1 - 2); i++)
                {
                    n1 = FLINT_MAX(0, i - len1 + 1);
                    n2 = FLINT_MIN(len1 - 1, (i + 1) / 2 - 1);

                    NMOD_REDC_FAST_DOT_3(cc, j, poly1[n1 + j], poly1[i - n2 + (n2 - n1 + 1) - 1 - j], n2 - n1 + 1, ctx);
                    cc = nmod_redc_fast_add(cc, cc, GR_NMOD_REDC_CTX(ctx));
                    if (i % 2 == 0 && i / 2 < len1)
                        cc = nmod_redc_fast_add(cc, nmod_redc_mul(poly1[i / 2], poly1[i / 2], GR_NMOD_REDC_CTX(ctx)), GR_NMOD_REDC_CTX(ctx));

                    res[i] = cc;
                }
            }

            if (n >= 2 * len1 - 1)
                res[2 * len1 - 2] = nmod_redc_fast_mul(poly1[len1 - 1], poly1[len1 - 1], GR_NMOD_REDC_CTX(ctx));
        }
        else
        {
            if (sbits < FLINT_BITS + nbits)
            {
                for (i = 0; i < n; i++)
                {
                    n1 = FLINT_MIN(len1 - 1, i);
                    n2 = FLINT_MIN(len2 - 1, i);
                    NMOD_REDC_DOT_1(res[i], j, poly1[i - n2 + j], poly2[i + n2 - i - j], n1 + n2 - i + 1, ctx);
                }
            }
            else if (sbits <= 2 * FLINT_BITS)
            {
                for (i = 0; i < n; i++)
                {
                    n1 = FLINT_MIN(len1 - 1, i);
                    n2 = FLINT_MIN(len2 - 1, i);
                    NMOD_REDC_DOT_2(res[i], j, poly1[i - n2 + j], poly2[i + n2 - i - j], n1 + n2 - i + 1, ctx);
                }
            }
            else
            {
                for (i = 0; i < n; i++)
                {
                    n1 = FLINT_MIN(len1 - 1, i);
                    n2 = FLINT_MIN(len2 - 1, i);
                    NMOD_REDC_DOT_3(res[i], j, poly1[i - n2 + j], poly2[i + n2 - i - j], n1 + n2 - i + 1, ctx);
                }
            }
        }

        return GR_SUCCESS;
    }

    if (poly1 == poly2 && len1 == len2)
    {
        len1 = FLINT_MIN(len1, n);
        alloc = len1;

        t1 = GR_TMP_ALLOC(alloc * sizeof(ulong));

        for (i = 0; i < len1; i++)
            t1[i] = nmod_redc_fast_normalise(poly1[i], GR_NMOD_REDC_CTX(ctx));

        if (len1 + len2 - 1 == n)
            _nmod_poly_mul(res, t1, len1, t1, len1, GR_NMOD_REDC_MOD(ctx));
        else
            _nmod_poly_mullow(res, t1, len1, t1, len1, n, GR_NMOD_REDC_MOD(ctx));
    }
    else
    {
        len1 = FLINT_MIN(len1, n);
        len2 = FLINT_MIN(len2, n);
        alloc = len1 + len2;

        t1 = GR_TMP_ALLOC(alloc * sizeof(ulong));
        t2 = t1 + len1;

        for (i = 0; i < len1; i++)
            t1[i] = nmod_redc_fast_normalise(poly1[i], GR_NMOD_REDC_CTX(ctx));
        for (i = 0; i < len2; i++)
            t2[i] = nmod_redc_fast_normalise(poly2[i], GR_NMOD_REDC_CTX(ctx));

        if (len1 + len2 - 1 == n)
        {
            if (len1 >= len2)
                _nmod_poly_mul(res, t1, len1, t2, len2, GR_NMOD_REDC_MOD(ctx));
            else
                _nmod_poly_mul(res, t2, len2, t1, len1, GR_NMOD_REDC_MOD(ctx));
        }
        else
        {
            if (len1 >= len2)
                _nmod_poly_mullow(res, t1, len1, t2, len2, n, GR_NMOD_REDC_MOD(ctx));
            else
                _nmod_poly_mullow(res, t2, len2, t1, len1, n, GR_NMOD_REDC_MOD(ctx));
        }
    }

    for (i = 0; i < n; i++)
        res[i] = nmod_redc_get_nmod(res[i], GR_NMOD_REDC_CTX(ctx));

    GR_TMP_FREE(t1, alloc * sizeof(ulong));

    return GR_SUCCESS;
}

static void _nmod_redc_fast_poly_divrem_q1_preinv1(nn_ptr Q, nn_ptr R,
                          nn_srcptr A, slong lenA, nn_srcptr B, slong lenB,
                          ulong invL, nmod_redc_ctx_t ctx)
{
    slong i;
    ulong t, q0, q1;

    FLINT_ASSERT(lenA == lenB + 1);
    FLINT_ASSERT(lenB >= 2);

    q1 = nmod_redc_fast_mul(A[lenA-1], invL, ctx);
    t  = nmod_redc_fast_mul(q1, B[lenB-2], ctx);
    t  = nmod_redc_fast_sub(t, A[lenA-2], ctx);
    q0 = nmod_redc_fast_mul(t, invL, ctx);
    Q[0] = nmod_redc_fast_neg(q0, ctx);
    Q[1] = q1;
    q1 = nmod_redc_fast_neg(q1, ctx);
    R[0] = nmod_redc_fast_add(A[0], nmod_redc_fast_mul(q0, B[0], ctx), ctx);

    ulong n = GR_NMOD_REDC_CTX(ctx)->mod.n;
    ulong nred = GR_NMOD_REDC_CTX(ctx)->nred;

    for (i = 1; i < lenB - 1; i++)
    {
        ull_t s, v;
        s = ull_u_mul_u(q1, B[i - 1]);
        v = ull_u_mul_u(q0, B[i]);
        s = ull_add(s, v);
        FLINT_ASSERT(ull_hi(s) < n);
        R[i] = _n_mod_add(A[i], n_ll_redc_fast(s, n, nred), 2 * n);
    }
}

static int _gr_nmod_redc_poly_divrem_q1(nn_ptr Q, nn_ptr R,
                          nn_srcptr A, slong lenA, nn_srcptr B, slong lenB,
                          gr_ctx_t ctx)
{
    ulong invB;
    int status;

    status = _gr_nmod_redc_inv(&invB, B + lenB - 1, ctx);

    if (status != GR_SUCCESS)
        return status;

    _nmod_redc_fast_poly_divrem_q1_preinv1(Q, R, A, lenA, B, lenB, invB, GR_NMOD_REDC_CTX(ctx));
    return GR_SUCCESS;
}

static int _gr_nmod_redc_fast_poly_divrem(nn_ptr Q, nn_ptr R, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, gr_ctx_t ctx)
{
    if (lenB == 1)
        return _gr_poly_divrem_basecase(Q, R, A, lenA, B, lenB, ctx);

    if (lenA == lenB + 1)
        return _gr_nmod_redc_poly_divrem_q1(Q, R, A, lenA, B, lenB, ctx);

    if (lenB >= 20 && lenA - lenB >= 20)
    {
        nn_ptr tA, tB;
        slong i, lenQ, lenR;
        int status;

        tA = GR_TMP_ALLOC((lenA + lenB) * sizeof(ulong));
        tB = tA + lenA;

        /* todo: avoid some conversions? */

        for (i = 0; i < lenA; i++)
        {
            tA[i] = nmod_redc_get_nmod(A[i], GR_NMOD_REDC_CTX(ctx));
            FLINT_ASSERT(tA[i] < GR_NMOD_REDC_N(ctx));
        }
        for (i = 0; i < lenB; i++)
        {
            tB[i] = nmod_redc_get_nmod(B[i], GR_NMOD_REDC_CTX(ctx));
            FLINT_ASSERT(tB[i] < GR_NMOD_REDC_N(ctx));
        }

        gr_ctx_t modctx;
        _gr_ctx_init_nmod(modctx, &GR_NMOD_REDC_MOD(ctx));
        status = _gr_poly_divrem(Q, R, tA, lenA, tB, lenB, modctx);

        lenQ = lenA - lenB + 1;
        lenR = lenB - 1;

        /* If unsuccessful, Q[i] or R[i] could contain invalid data for
           nmod_redc_set_nmod. */
        if (status == GR_SUCCESS)
        {
            for (i = 0; i < lenQ; i++)
            {
                FLINT_ASSERT(Q[i] < GR_NMOD_REDC_N(ctx));
                Q[i] = nmod_redc_set_nmod(Q[i], GR_NMOD_REDC_CTX(ctx));
            }
            for (i = 0; i < lenR; i++)
            {
                FLINT_ASSERT(R[i] < GR_NMOD_REDC_N(ctx));
                R[i] = nmod_redc_set_nmod(R[i], GR_NMOD_REDC_CTX(ctx));
            }
        }

        GR_TMP_FREE(tA, (lenA + lenB) * sizeof(ulong));

        return status;
    }
    else
    {
        return _gr_poly_divrem_newton(Q, R, A, lenA, B, lenB, ctx);
    }
}



int __gr_nmod_redc_methods_initialized = 0;
int __gr_nmod_redc_fast_methods_initialized = 0;

gr_static_method_table __gr_nmod_redc_methods;
gr_static_method_table __gr_nmod_redc_fast_methods;

#pragma GCC diagnostic ignored "-Wcast-function-type"
#define GR_FUNCPTR_CAST (gr_funcptr)

gr_method_tab_input __gr_nmod_redc_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,       GR_FUNCPTR_CAST _gr_nmod_redc_ctx_write},
    {GR_METHOD_CTX_IS_RING,     GR_FUNCPTR_CAST gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, GR_FUNCPTR_CAST gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  GR_FUNCPTR_CAST _gr_nmod_redc_ctx_is_field},
    {GR_METHOD_CTX_IS_FIELD,            GR_FUNCPTR_CAST _gr_nmod_redc_ctx_is_field},
    {GR_METHOD_CTX_IS_FINITE,
                                GR_FUNCPTR_CAST gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC,
                                GR_FUNCPTR_CAST gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_EXACT,    GR_FUNCPTR_CAST gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_CANONICAL,
                                GR_FUNCPTR_CAST gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_SET_IS_FIELD,GR_FUNCPTR_CAST _gr_nmod_redc_ctx_set_is_field},
    {GR_METHOD_INIT,            GR_FUNCPTR_CAST _gr_nmod_redc_init},
    {GR_METHOD_CLEAR,           GR_FUNCPTR_CAST _gr_nmod_redc_clear},
    {GR_METHOD_SWAP,            GR_FUNCPTR_CAST _gr_nmod_redc_swap},
    {GR_METHOD_SET_SHALLOW,     GR_FUNCPTR_CAST _gr_nmod_redc_set_shallow},
    {GR_METHOD_RANDTEST,        GR_FUNCPTR_CAST _gr_nmod_redc_randtest},
    {GR_METHOD_WRITE,           GR_FUNCPTR_CAST _gr_nmod_redc_write},
    {GR_METHOD_ZERO,            GR_FUNCPTR_CAST _gr_nmod_redc_zero},
    {GR_METHOD_ONE,             GR_FUNCPTR_CAST _gr_nmod_redc_one},
    {GR_METHOD_IS_ZERO,         GR_FUNCPTR_CAST _gr_nmod_redc_is_zero},
    {GR_METHOD_IS_ONE,          GR_FUNCPTR_CAST _gr_nmod_redc_is_one},
    {GR_METHOD_IS_NEG_ONE,      GR_FUNCPTR_CAST _gr_nmod_redc_is_neg_one},
    {GR_METHOD_EQUAL,           GR_FUNCPTR_CAST _gr_nmod_redc_equal},
    {GR_METHOD_SET,             GR_FUNCPTR_CAST _gr_nmod_redc_set},
    {GR_METHOD_SET_SI,          GR_FUNCPTR_CAST _gr_nmod_redc_set_si},
    {GR_METHOD_SET_UI,          GR_FUNCPTR_CAST _gr_nmod_redc_set_ui},
    {GR_METHOD_SET_FMPZ,        GR_FUNCPTR_CAST _gr_nmod_redc_set_fmpz},
    {GR_METHOD_SET_OTHER,       GR_FUNCPTR_CAST _gr_nmod_redc_set_other},
    {GR_METHOD_GET_FMPZ,        GR_FUNCPTR_CAST _gr_nmod_redc_get_fmpz},
    {GR_METHOD_NEG,             GR_FUNCPTR_CAST _gr_nmod_redc_neg},
    {GR_METHOD_ADD,             GR_FUNCPTR_CAST _gr_nmod_redc_add},
/*
    {GR_METHOD_ADD_SI,          GR_FUNCPTR_CAST _gr_nmod_redc_add_si},
    {GR_METHOD_ADD_UI,          GR_FUNCPTR_CAST _gr_nmod_redc_add_ui},
*/
    {GR_METHOD_SUB,             GR_FUNCPTR_CAST _gr_nmod_redc_sub},
/*
    {GR_METHOD_SUB_SI,          GR_FUNCPTR_CAST _gr_nmod_redc_sub_si},
    {GR_METHOD_SUB_UI,          GR_FUNCPTR_CAST _gr_nmod_redc_sub_ui},
*/
    {GR_METHOD_MUL,             GR_FUNCPTR_CAST _gr_nmod_redc_mul},
/*
    {GR_METHOD_MUL_SI,          GR_FUNCPTR_CAST _gr_nmod_redc_mul_si},
    {GR_METHOD_MUL_UI,          GR_FUNCPTR_CAST _gr_nmod_redc_mul_ui},
    {GR_METHOD_MUL_FMPZ,        GR_FUNCPTR_CAST _gr_nmod_redc_mul_fmpz},
    {GR_METHOD_ADDMUL,          GR_FUNCPTR_CAST _gr_nmod_redc_addmul},
    {GR_METHOD_SUBMUL,          GR_FUNCPTR_CAST _gr_nmod_redc_submul},
    {GR_METHOD_MUL_TWO,         GR_FUNCPTR_CAST _gr_nmod_redc_mul_two},
    {GR_METHOD_MUL_2EXP_SI,     GR_FUNCPTR_CAST _gr_nmod_redc_mul_2exp_si},
*/
    {GR_METHOD_SQR,             GR_FUNCPTR_CAST _gr_nmod_redc_sqr},
    {GR_METHOD_DIV,             GR_FUNCPTR_CAST _gr_nmod_redc_div},
/*
    {GR_METHOD_DIV_SI,          GR_FUNCPTR_CAST _gr_nmod_redc_div_si},
    {GR_METHOD_DIV_UI,          GR_FUNCPTR_CAST _gr_nmod_redc_div_ui},
    {GR_METHOD_DIV_FMPZ,        GR_FUNCPTR_CAST _gr_nmod_redc_div_fmpz},
    {GR_METHOD_DIV_NONUNIQUE,   GR_FUNCPTR_CAST _gr_nmod_redc_div_nonunique},
    {GR_METHOD_DIVIDES,         GR_FUNCPTR_CAST _gr_nmod_redc_divides},
    {GR_METHOD_IS_INVERTIBLE,   GR_FUNCPTR_CAST _gr_nmod_redc_is_invertible},
*/
    {GR_METHOD_INV,             GR_FUNCPTR_CAST _gr_nmod_redc_inv},
    {GR_METHOD_POW_UI,          GR_FUNCPTR_CAST _gr_nmod_redc_pow_ui},
/*
    {GR_METHOD_POW_SI,          GR_FUNCPTR_CAST _gr_nmod_redc_pow_si},
    {GR_METHOD_POW_FMPZ,        GR_FUNCPTR_CAST _gr_nmod_redc_pow_fmpz},
    {GR_METHOD_IS_SQUARE,       GR_FUNCPTR_CAST _gr_nmod_redc_is_square},
    {GR_METHOD_SQRT,            GR_FUNCPTR_CAST _gr_nmod_redc_sqrt},
*/
    {GR_METHOD_VEC_INIT,        GR_FUNCPTR_CAST _gr_nmod_redc_vec_init},
    {GR_METHOD_VEC_CLEAR,       GR_FUNCPTR_CAST _gr_nmod_redc_vec_clear},
    {GR_METHOD_VEC_SET,         GR_FUNCPTR_CAST _gr_nmod_redc_vec_set},
    {GR_METHOD_VEC_NORMALISE,   GR_FUNCPTR_CAST _gr_nmod_redc_vec_normalise},
    {GR_METHOD_VEC_NORMALISE_WEAK,   GR_FUNCPTR_CAST _gr_nmod_redc_vec_normalise_weak},
    {GR_METHOD_VEC_NEG,         GR_FUNCPTR_CAST _gr_nmod_redc_vec_neg},
    {GR_METHOD_VEC_ADD,         GR_FUNCPTR_CAST _gr_nmod_redc_vec_add},
    {GR_METHOD_VEC_SUB,         GR_FUNCPTR_CAST _gr_nmod_redc_vec_sub},
    {GR_METHOD_VEC_MUL,         GR_FUNCPTR_CAST _gr_nmod_redc_vec_mul},
    {GR_METHOD_VEC_MUL_SCALAR,      GR_FUNCPTR_CAST _gr_nmod_redc_vec_mul_scalar},
    {GR_METHOD_VEC_MUL_SCALAR_SI,   GR_FUNCPTR_CAST _gr_nmod_redc_vec_mul_scalar_si},
    {GR_METHOD_VEC_MUL_SCALAR_UI,   GR_FUNCPTR_CAST _gr_nmod_redc_vec_mul_scalar_ui},
    {GR_METHOD_VEC_MUL_SCALAR_FMPZ, GR_FUNCPTR_CAST _gr_nmod_redc_vec_mul_scalar_fmpz},
/*
    {GR_METHOD_VEC_MUL_SCALAR_2EXP_SI,   GR_FUNCPTR_CAST _gr_nmod_redc_vec_mul_scalar_2exp_si},
*/

    {GR_METHOD_SCALAR_MUL_VEC,      GR_FUNCPTR_CAST _gr_nmod_redc_scalar_mul_vec},
    {GR_METHOD_VEC_ADDMUL_SCALAR,        GR_FUNCPTR_CAST _gr_nmod_redc_vec_addmul_scalar},
    {GR_METHOD_VEC_ADDMUL_SCALAR_SI,     GR_FUNCPTR_CAST _gr_nmod_redc_vec_addmul_scalar_si},
    {GR_METHOD_VEC_SUBMUL_SCALAR,        GR_FUNCPTR_CAST _gr_nmod_redc_vec_submul_scalar},
    {GR_METHOD_VEC_SUBMUL_SCALAR_SI,     GR_FUNCPTR_CAST _gr_nmod_redc_vec_submul_scalar_si},
/*
    {GR_METHOD_VEC_SUM,         GR_FUNCPTR_CAST _gr_nmod_redc_vec_sum},
*/
    {GR_METHOD_VEC_PRODUCT,     GR_FUNCPTR_CAST _gr_nmod_redc_vec_product},
    {GR_METHOD_VEC_DOT,         GR_FUNCPTR_CAST _gr_nmod_redc_vec_dot},
    {GR_METHOD_VEC_DOT_REV,     GR_FUNCPTR_CAST _gr_nmod_redc_vec_dot_rev},
/*
    {GR_METHOD_VEC_RECIPROCALS, GR_FUNCPTR_CAST _gr_nmod_redc_vec_reciprocals},
*/
    {GR_METHOD_POLY_MULLOW,     GR_FUNCPTR_CAST _gr_nmod_redc_poly_mullow},
/*
    {GR_METHOD_POLY_DIVREM,     GR_FUNCPTR_CAST _gr_nmod_redc_poly_divrem},
*/
/*
    {GR_METHOD_POLY_DIVEXACT,   GR_FUNCPTR_CAST _gr_nmod_redc_poly_divexact},
    {GR_METHOD_POLY_INV_SERIES, GR_FUNCPTR_CAST _gr_nmod_redc_poly_inv_series},
    {GR_METHOD_POLY_INV_SERIES_BASECASE, GR_FUNCPTR_CAST _gr_redc_nmod_poly_inv_series_basecase},
    {GR_METHOD_POLY_DIV_SERIES, GR_FUNCPTR_CAST _gr_nmod_redc_poly_div_series},
    {GR_METHOD_POLY_DIV_SERIES_BASECASE, GR_FUNCPTR_CAST _gr_nmod_redc_poly_div_series_basecase},
    {GR_METHOD_POLY_RSQRT_SERIES, GR_FUNCPTR_CAST _gr_nmod_redc_poly_rsqrt_series},
    {GR_METHOD_POLY_SQRT_SERIES,  GR_FUNCPTR_CAST _gr_nmod_redc_poly_sqrt_series},
    {GR_METHOD_POLY_EXP_SERIES,  GR_FUNCPTR_CAST _gr_nmod_redc_poly_exp_series},
    {GR_METHOD_POLY_ROOTS,      GR_FUNCPTR_CAST _gr_nmod_redc_roots_gr_poly},
    {GR_METHOD_MAT_MUL,         GR_FUNCPTR_CAST _gr_nmod_redc_mat_mul},
*/

    {0,                         GR_FUNCPTR_CAST NULL},
};

gr_method_tab_input __gr_nmod_redc_fast_methods_input[] =
{
    {GR_METHOD_IS_ZERO,         GR_FUNCPTR_CAST _gr_nmod_redc_fast_is_zero},
    {GR_METHOD_IS_ONE,          GR_FUNCPTR_CAST _gr_nmod_redc_fast_is_one},
    {GR_METHOD_IS_NEG_ONE,      GR_FUNCPTR_CAST _gr_nmod_redc_fast_is_neg_one},
    {GR_METHOD_EQUAL,           GR_FUNCPTR_CAST _gr_nmod_redc_fast_equal},
    {GR_METHOD_NEG,             GR_FUNCPTR_CAST _gr_nmod_redc_fast_neg},
    {GR_METHOD_ADD,             GR_FUNCPTR_CAST _gr_nmod_redc_fast_add},
    {GR_METHOD_SUB,             GR_FUNCPTR_CAST _gr_nmod_redc_fast_sub},
    {GR_METHOD_MUL,             GR_FUNCPTR_CAST _gr_nmod_redc_fast_mul},
    {GR_METHOD_DIV,             GR_FUNCPTR_CAST _gr_nmod_redc_fast_div},
    {GR_METHOD_SQR,             GR_FUNCPTR_CAST _gr_nmod_redc_fast_sqr},
    {GR_METHOD_POW_UI,          GR_FUNCPTR_CAST _gr_nmod_redc_fast_pow_ui},
    {GR_METHOD_VEC_NORMALISE,   GR_FUNCPTR_CAST _gr_nmod_redc_fast_vec_normalise},
    {GR_METHOD_VEC_NORMALISE_WEAK,   GR_FUNCPTR_CAST _gr_nmod_redc_fast_vec_normalise_weak},
    {GR_METHOD_VEC_NEG,         GR_FUNCPTR_CAST _gr_nmod_redc_fast_vec_neg},
    {GR_METHOD_VEC_ADD,         GR_FUNCPTR_CAST _gr_nmod_redc_fast_vec_add},
    {GR_METHOD_VEC_SUB,         GR_FUNCPTR_CAST _gr_nmod_redc_fast_vec_sub},
    {GR_METHOD_VEC_MUL,         GR_FUNCPTR_CAST _gr_nmod_redc_fast_vec_mul},
    {GR_METHOD_VEC_MUL_SCALAR,      GR_FUNCPTR_CAST _gr_nmod_redc_fast_vec_mul_scalar},
    {GR_METHOD_SCALAR_MUL_VEC,      GR_FUNCPTR_CAST _gr_nmod_redc_fast_scalar_mul_vec},
    {GR_METHOD_VEC_MUL_SCALAR_SI,   GR_FUNCPTR_CAST _gr_nmod_redc_fast_vec_mul_scalar_si},
    {GR_METHOD_VEC_MUL_SCALAR_UI,   GR_FUNCPTR_CAST _gr_nmod_redc_fast_vec_mul_scalar_ui},
    {GR_METHOD_VEC_MUL_SCALAR_FMPZ, GR_FUNCPTR_CAST _gr_nmod_redc_fast_vec_mul_scalar_fmpz},
    {GR_METHOD_VEC_ADDMUL_SCALAR,        GR_FUNCPTR_CAST _gr_nmod_redc_fast_vec_addmul_scalar},
    {GR_METHOD_VEC_ADDMUL_SCALAR_SI,     GR_FUNCPTR_CAST _gr_nmod_redc_fast_vec_addmul_scalar_si},
    {GR_METHOD_VEC_SUBMUL_SCALAR,        GR_FUNCPTR_CAST _gr_nmod_redc_fast_vec_submul_scalar},
    {GR_METHOD_VEC_SUBMUL_SCALAR_SI,     GR_FUNCPTR_CAST _gr_nmod_redc_fast_vec_submul_scalar_si},
    {GR_METHOD_VEC_PRODUCT,     GR_FUNCPTR_CAST _gr_nmod_redc_fast_vec_product},
    {GR_METHOD_VEC_DOT,         GR_FUNCPTR_CAST _gr_nmod_redc_fast_vec_dot},
    {GR_METHOD_VEC_DOT_REV,     GR_FUNCPTR_CAST _gr_nmod_redc_fast_vec_dot_rev},
    {GR_METHOD_POLY_MULLOW,     GR_FUNCPTR_CAST _gr_nmod_redc_fast_poly_mullow},
    {GR_METHOD_POLY_DIVREM,     GR_FUNCPTR_CAST _gr_nmod_redc_fast_poly_divrem},
    {0,                         GR_FUNCPTR_CAST NULL},
};

int
gr_ctx_init_nmod_redc(gr_ctx_t ctx, ulong n)
{
    if (n == 0)
        return GR_DOMAIN;
    /* Exclude mod 1 to avoid worrying about some degenerate cases. */
    if (n % 2 == 0 || n == 1)
        return GR_UNABLE;

    ctx->which_ring = GR_CTX_NMOD_REDC;
    ctx->sizeof_elem = sizeof(ulong);
    ctx->size_limit = WORD_MAX;

    nmod_redc_ctx_init_ui(GR_NMOD_REDC_CTX(ctx), n);
    GR_NMOD_REDC_ONE(ctx) = nmod_redc_set_ui(1, GR_NMOD_REDC_CTX(ctx));
    GR_NMOD_REDC_IS_PRIME(ctx) = T_UNKNOWN;

    ctx->methods = __gr_nmod_redc_methods;

    if (!__gr_nmod_redc_methods_initialized)
    {
        gr_method_tab_init(__gr_nmod_redc_methods, __gr_nmod_redc_methods_input);
        __gr_nmod_redc_methods_initialized = 1;
    }

    return GR_SUCCESS;
}

int
gr_ctx_init_nmod_redc_fast(gr_ctx_t ctx, ulong n)
{
    if (n == 0)
        return GR_DOMAIN;
    /* Exclude mod 1 to avoid worrying about some degenerate cases. */
    if (n % 2 == 0 || n == 1)
        return GR_UNABLE;
    if (n >= UWORD(1) << (FLINT_BITS - 2))
        return GR_UNABLE;

    ctx->which_ring = GR_CTX_NMOD_REDC_FAST;
    ctx->sizeof_elem = sizeof(ulong);
    ctx->size_limit = WORD_MAX;

    nmod_redc_ctx_init_ui(GR_NMOD_REDC_CTX(ctx), n);
    GR_NMOD_REDC_ONE(ctx) = nmod_redc_set_ui(1, GR_NMOD_REDC_CTX(ctx));
    GR_NMOD_REDC_IS_PRIME(ctx) = T_UNKNOWN;

    ctx->methods = __gr_nmod_redc_fast_methods;

    if (!__gr_nmod_redc_fast_methods_initialized)
    {
        gr_method_tab_init(__gr_nmod_redc_fast_methods, __gr_nmod_redc_methods_input);
        gr_method_tab_extend(__gr_nmod_redc_fast_methods, __gr_nmod_redc_fast_methods_input);
        __gr_nmod_redc_fast_methods_initialized = 1;
    }

    return GR_SUCCESS;
}

