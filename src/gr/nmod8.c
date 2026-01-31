/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdint.h>
#include "fmpz.h"
//#include "fmpq.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "gr.h"
#include "gr_poly.h"
#include "gr_mat.h"

#define NMOD8_CTX_REF(ring_ctx) (((nmod_t *)((ring_ctx))))
#define NMOD8_CTX(ring_ctx) (*NMOD8_CTX_REF(ring_ctx))

typedef uint8_t nmod8_t[1];

static void
nmod8_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Integers mod ");
    gr_stream_write_ui(out, NMOD8_CTX(ctx).n);
    gr_stream_write(out, " (nmod8)");
}

/* todo: n_is_prime is fast, but this should still be cached
   or use a fixed table lookup */
static truth_t
nmod8_ctx_is_field(const gr_ctx_t ctx)
{
    return n_is_prime(NMOD8_CTX(ctx).n) ? T_TRUE : T_FALSE;
}

static void
nmod8_init(nmod8_t x, const gr_ctx_t ctx)
{
    x[0] = 0;
}

static void
nmod8_clear(nmod8_t x, const gr_ctx_t ctx)
{
}

static void
nmod8_swap(nmod8_t x, nmod8_t y, const gr_ctx_t ctx)
{
    nmod8_t t;
    *t = *x;
    *x = *y;
    *y = *t;
}

static void
nmod8_set_shallow(nmod8_t res, const nmod8_t x, const gr_ctx_t ctx)
{
    *res = *x;
}

static int
nmod8_randtest(nmod8_t res, flint_rand_t state, const gr_ctx_t ctx)
{
    res[0] = n_randtest(state) % NMOD8_CTX(ctx).n;
    return GR_SUCCESS;
}

static int
nmod8_write(gr_stream_t out, const nmod8_t x, const gr_ctx_t ctx)
{
    gr_stream_write_ui(out, x[0]);
    return GR_SUCCESS;
}

static int
nmod8_zero(nmod8_t x, const gr_ctx_t ctx)
{
    x[0] = 0;
    return GR_SUCCESS;
}

static int
nmod8_one(nmod8_t x, const gr_ctx_t ctx)
{
    x[0] = (NMOD8_CTX(ctx).n != 1);
    return GR_SUCCESS;
}

static int
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

static int
nmod8_set_ui(nmod8_t res, ulong v, const gr_ctx_t ctx)
{
    ulong t;
    nmod_t mod = NMOD8_CTX(ctx);
    NMOD_RED(t, v, mod);
    res[0] = t;
    return GR_SUCCESS;
}

static int
nmod8_set_fmpz(nmod8_t res, const fmpz_t v, const gr_ctx_t ctx)
{
    nmod_t mod = NMOD8_CTX(ctx);
    res[0] = fmpz_get_nmod(v, mod);
    return GR_SUCCESS;
}

static truth_t
nmod8_is_zero(const nmod8_t x, const gr_ctx_t ctx)
{
    return (x[0] == 0) ? T_TRUE : T_FALSE;
}

static truth_t
nmod8_is_one(const nmod8_t x, const gr_ctx_t ctx)
{
    return (x[0] == (NMOD8_CTX(ctx).n != 1)) ? T_TRUE : T_FALSE;
}

static truth_t
nmod8_is_neg_one(const nmod8_t x, const gr_ctx_t ctx)
{
    return (x[0] == NMOD8_CTX(ctx).n - 1) ? T_TRUE : T_FALSE;
}

static truth_t
nmod8_equal(const nmod8_t x, const nmod8_t y, const gr_ctx_t ctx)
{
    return (x[0] == y[0]) ? T_TRUE : T_FALSE;
}

static int
nmod8_set(nmod8_t res, const nmod8_t x, const gr_ctx_t ctx)
{
    res[0] = x[0];
    return GR_SUCCESS;
}

static int
nmod8_neg(nmod8_t res, const nmod8_t x, const gr_ctx_t ctx)
{
    res[0] = nmod_neg(x[0], NMOD8_CTX(ctx));
    return GR_SUCCESS;
}

static int
nmod8_add(nmod8_t res, const nmod8_t x, const nmod8_t y, const gr_ctx_t ctx)
{
    res[0] = nmod_add(x[0], y[0], NMOD8_CTX(ctx));
    return GR_SUCCESS;
}

static int
nmod8_add_si(nmod8_t res, const nmod8_t x, slong y, const gr_ctx_t ctx)
{
    nmod8_t t;
    nmod8_set_si(t, y, ctx);
    return nmod8_add(res, x, t, ctx);
}

static int
nmod8_sub(nmod8_t res, const nmod8_t x, const nmod8_t y, const gr_ctx_t ctx)
{
    res[0] = nmod_sub(x[0], y[0], NMOD8_CTX(ctx));
    return GR_SUCCESS;
}

static int
nmod8_mul(nmod8_t res, const nmod8_t x, const nmod8_t y, const gr_ctx_t ctx)
{
    res[0] = nmod_mul(x[0], y[0], NMOD8_CTX(ctx));
    return GR_SUCCESS;
}

static int
nmod8_mul_si(nmod8_t res, const nmod8_t x, slong y, const gr_ctx_t ctx)
{
    nmod8_t t;
    nmod8_set_si(t, y, ctx);
    return nmod8_mul(res, x, t, ctx);
}

static int
nmod8_addmul(nmod8_t res, const nmod8_t x, const nmod8_t y, const gr_ctx_t ctx)
{
    ulong r = res[0];
    NMOD_ADDMUL(r, x[0], y[0], NMOD8_CTX(ctx));
    res[0] = r;
    return GR_SUCCESS;
}

static int
nmod8_submul(nmod8_t res, const nmod8_t x, const nmod8_t y, const gr_ctx_t ctx)
{
    ulong r = res[0];
    ulong t = nmod_neg(y[0], NMOD8_CTX(ctx));
    NMOD_ADDMUL(r, x[0], t, NMOD8_CTX(ctx));
    res[0] = r;
    return GR_SUCCESS;
}

static int
nmod8_mul_two(nmod8_t res, const nmod8_t x, const gr_ctx_t ctx)
{
    return nmod8_add(res, x, x, ctx);
}

static int
nmod8_sqr(nmod8_t res, const nmod8_t x, const gr_ctx_t ctx)
{
    return nmod8_mul(res, x, x, ctx);
}

static int
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

static int
nmod8_div(nmod8_t res, const nmod8_t x, const nmod8_t y, const gr_ctx_t ctx)
{
    nmod8_t t;
    int status;

    status = nmod8_inv(t, y, ctx);

    if (status == GR_SUCCESS)
        nmod8_mul(res, x, t, ctx);

    return status;
}

static int
nmod8_div_si(nmod8_t res, const nmod8_t x, slong y, const gr_ctx_t ctx)
{
    nmod8_t t;
    nmod8_set_si(t, y, ctx);
    return nmod8_div(res, x, t, ctx);
}

static int
nmod8_div_ui(nmod8_t res, const nmod8_t x, ulong y, const gr_ctx_t ctx)
{
    nmod8_t t;
    nmod8_set_ui(t, y, ctx);
    return nmod8_div(res, x, t, ctx);
}

static int
nmod8_div_fmpz(nmod8_t res, const nmod8_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    nmod8_t t;
    nmod8_set_fmpz(t, y, ctx);
    return nmod8_div(res, x, t, ctx);
}

static truth_t
nmod8_is_invertible(const nmod8_t x, const gr_ctx_t ctx)
{
    ulong r, g;
    g = n_gcdinv(&r, x[0], NMOD8_CTX(ctx).n);
    return (g == 1) ? T_TRUE : T_FALSE;
}

static truth_t
nmod8_divides(const nmod8_t x, const nmod8_t y, const gr_ctx_t ctx)
{
    ulong t;
    return nmod_divides(&t, y[0], x[0], NMOD8_CTX(ctx)) ? T_TRUE : T_FALSE;
}

static int
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

static void
_nmod8_vec_init(uint8_t * res, slong len, gr_ctx_t ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        res[i] = 0;
}

static void
_nmod8_vec_clear(uint8_t * res, slong len, gr_ctx_t ctx)
{
}

static int
_nmod8_vec_zero(uint8_t * res, slong len, gr_ctx_t ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        res[i] = 0;

    return GR_SUCCESS;
}

static int
_nmod8_vec_set(uint8_t * res, const uint8_t * vec, slong len, gr_ctx_t ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        res[i] = vec[i];

    return GR_SUCCESS;
}

static uint8_t _nmod8_add(uint8_t a, uint8_t b, uint8_t n)
{
   const uint8_t sum = a + b;
   return sum - n + ((((int8_t)(sum - n))>>7) & n);
}

static uint8_t _nmod8_sub(uint8_t a, uint8_t b, uint8_t n)
{
   const uint8_t diff = a - b;
   return  ((((int8_t)diff)>>7) & n) + diff;
}

static uint16_t _nmod16_add(uint16_t a, uint16_t b, uint16_t n)
{
   const uint16_t sum = a + b;
   return sum - n + ((((int16_t)(sum - n))>>15) & n);
}

static uint16_t _nmod16_sub(uint16_t a, uint16_t b, uint16_t n)
{
   const uint16_t diff = a - b;
   return  ((((int16_t)diff)>>15) & n) + diff;
}

static int
_nmod8_vec_neg(uint8_t * res, const uint8_t * vec, slong len, gr_ctx_t ctx)
{
    slong i;
    nmod_t mod = NMOD8_CTX(ctx);

    if (NMOD_BITS(mod) == 8)
    {
        for (i = 0; i < len; i++)
            res[i] = _nmod16_sub(0, vec[i], mod.n);
    }
    else
    {
        for (i = 0; i < len; i++)
            res[i] = _nmod8_sub(0, vec[i], mod.n);
    }

    return GR_SUCCESS;
}

static int
_nmod8_vec_add(uint8_t * res, const uint8_t * vec1, const uint8_t * vec2, slong len, gr_ctx_t ctx)
{
    slong i;
    nmod_t mod = NMOD8_CTX(ctx);

    if (NMOD_BITS(mod) == 8)
    {
        for (i = 0; i < len; i++)
            res[i] = _nmod16_add(vec1[i], vec2[i], mod.n);
    }
    else
    {
        for (i = 0; i < len; i++)
            res[i] = _nmod8_add(vec1[i], vec2[i], mod.n);
    }

    return GR_SUCCESS;
}

static int
_nmod8_vec_sub(uint8_t * res, const uint8_t * vec1, const uint8_t * vec2, slong len, gr_ctx_t ctx)
{
    slong i;
    nmod_t mod = NMOD8_CTX(ctx);

    if (NMOD_BITS(mod) == 8)
    {
        for (i = 0; i < len; i++)
            res[i] = _nmod16_sub(vec1[i], vec2[i], mod.n);
    }
    else
    {
        for (i = 0; i < len; i++)
            res[i] = _nmod8_sub(vec1[i], vec2[i], mod.n);
    }

    return GR_SUCCESS;
}

static int
_nmod8_vec_mul(uint8_t * res, const uint8_t * vec1, const uint8_t * vec2, slong len, gr_ctx_t ctx)
{
    slong i;
    nmod_t mod = NMOD8_CTX(ctx);

    for (i = 0; i < len; i++)
        res[i] = nmod_mul(vec1[i], vec2[i], mod);

    return GR_SUCCESS;
}

static int
_nmod8_vec_dot(nmod8_t res, const nmod8_t initial, int subtract, const uint8_t * vec1, const uint8_t * vec2, slong len, gr_ctx_t ctx)
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

static int
_nmod8_vec_dot_rev(nmod8_t res, const nmod8_t initial, int subtract, const uint8_t * vec1, const uint8_t * vec2, slong len, gr_ctx_t ctx)
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

static uint32_t _mod1_preinvert(uint32_t n)
{
    return (~((uint32_t) 0)) / n + 1;
}

static uint32_t _mod1(uint32_t x, uint32_t n, uint32_t ninv)
{
    uint32_t l = ninv * x;
    return (((uint64_t) l) * n) >> 32;
}

static uint32_t _mulmod1(uint32_t a, uint32_t b, uint32_t n, uint32_t ninv)
{
    return _mod1(a * b, n, ninv);
}

static uint32_t _addmulmod1(uint32_t s, uint32_t a, uint32_t b, uint32_t n, uint32_t ninv)
{
    return _mod1(s + a * b, n, ninv);
}

static int
_nmod8_vec_mul_scalar(uint8_t * res, const uint8_t * vec, slong len, const uint8_t * c, gr_ctx_t ctx)
{
    ulong d = *c;
    slong i;

    if (len > 5)
    {
        uint32_t n = NMOD8_CTX(ctx).n, ninv = _mod1_preinvert(NMOD8_CTX(ctx).n);

        for (i = 0; i < len; i++)
            res[i] = _mulmod1(vec[i], d, n, ninv);
    }
    else
    {
        for (i = 0; i < len; i++)
            res[i] = nmod_mul(vec[i], d, NMOD8_CTX(ctx));
    }

    return GR_SUCCESS;
}

static int
_nmod8_vec_addmul_scalar(uint8_t * res, const uint8_t * vec, slong len, const uint8_t * c, gr_ctx_t ctx)
{
    ulong d = *c;
    slong i;

    if (len > 5)
    {
        uint32_t n = NMOD8_CTX(ctx).n, ninv = _mod1_preinvert(NMOD8_CTX(ctx).n);

        for (i = 0; i < len; i++)
            res[i] = _addmulmod1(res[i], vec[i], d, n, ninv);
    }
    else
    {
        for (i = 0; i < len; i++)
            res[i] = nmod_addmul(res[i], vec[i], d, NMOD8_CTX(ctx));
    }

    return GR_SUCCESS;
}

static int
_nmod8_vec_submul_scalar(uint8_t * res, const uint8_t * vec, slong len, const uint8_t * c, gr_ctx_t ctx)
{
    uint8_t d;
    d = nmod_neg(*c, NMOD8_CTX(ctx));
    return _nmod8_vec_addmul_scalar(res, vec, len, &d, ctx);
}

static int
_nmod8_poly_mullow(uint8_t * res, const uint8_t * A, slong Alen, const uint8_t * B, slong Blen, slong trunc, gr_ctx_t ctx)
{
    nn_ptr TR, TA, TB;
    slong i, alloc;
    int squaring;

    Alen = FLINT_MIN(Alen, trunc);
    Blen = FLINT_MIN(Blen, trunc);

    /* todo: tune this */
    if (Alen < 10 || Blen < 10)
        return _gr_poly_mullow_classical(res, A, Alen, B, Blen, trunc, ctx);

    squaring = (A == B) && (Alen == Blen);

    alloc = squaring ? (Alen + trunc) : (Alen + Blen + trunc);
    alloc *= sizeof(ulong);
    TR = GR_TMP_ALLOC(alloc);
    TA = TR + trunc;
    TB = TA + Alen;

    for (i = 0; i < Alen; i++)
        TA[i] = A[i];

    if (squaring)
    {
        _nmod_poly_mullow(TR, TA, Alen, TA, Alen, trunc, NMOD8_CTX(ctx));
    }
    else
    {
        for (i = 0; i < Blen; i++)
            TB[i] = B[i];

        if (Alen >= Blen)
            _nmod_poly_mullow(TR, TA, Alen, TB, Blen, trunc, NMOD8_CTX(ctx));
        else
            _nmod_poly_mullow(TR, TB, Blen, TA, Alen, trunc, NMOD8_CTX(ctx));
    }

    for (i = 0; i < trunc; i++)
        res[i] = TR[i];

    GR_TMP_FREE(TR, alloc);
    return GR_SUCCESS;
}

/* todo: tuning for rectangular matrices */
static int
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
    {GR_METHOD_VEC_ZERO,        (gr_funcptr) _nmod8_vec_zero},
    {GR_METHOD_VEC_NEG,         (gr_funcptr) _nmod8_vec_neg},
    {GR_METHOD_VEC_ADD,         (gr_funcptr) _nmod8_vec_add},
    {GR_METHOD_VEC_SUB,         (gr_funcptr) _nmod8_vec_sub},
    {GR_METHOD_VEC_MUL,         (gr_funcptr) _nmod8_vec_mul},
    {GR_METHOD_VEC_MUL_SCALAR,     (gr_funcptr) _nmod8_vec_mul_scalar},
    {GR_METHOD_VEC_ADDMUL_SCALAR,  (gr_funcptr) _nmod8_vec_addmul_scalar},
    {GR_METHOD_VEC_SUBMUL_SCALAR,  (gr_funcptr) _nmod8_vec_submul_scalar},
    {GR_METHOD_VEC_DOT,         (gr_funcptr) _nmod8_vec_dot},
    {GR_METHOD_VEC_DOT_REV,     (gr_funcptr) _nmod8_vec_dot_rev},
    {GR_METHOD_POLY_MULLOW,     (gr_funcptr) _nmod8_poly_mullow},
    {GR_METHOD_MAT_MUL,         (gr_funcptr) _nmod8_mat_mul},
    {0,                         (gr_funcptr) NULL},
};

int
gr_ctx_init_nmod8(gr_ctx_t ctx, ulong n)
{
    if (n == 0)
        return GR_DOMAIN;
    if (n > 255)
        return GR_UNABLE;

    ctx->which_ring = GR_CTX_NMOD8;
    ctx->sizeof_elem = sizeof(uint8_t);
    ctx->size_limit = WORD_MAX;

    nmod_init(NMOD8_CTX_REF(ctx), n);

    ctx->methods = _nmod8_methods;

    if (!_nmod8_methods_initialized)
    {
        gr_method_tab_init(_nmod8_methods, _nmod8_methods_input);
        _nmod8_methods_initialized = 1;
    }

    return GR_SUCCESS;
}
