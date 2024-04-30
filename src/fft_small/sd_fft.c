/*
    Copyright (C) 2022 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fft_small.h"

/*
    N is supposed to be a good fit for the number of points to process per loop
    in the radix 4 butterflies.

        16x 4-wide AVX registers  => N = 8
        32x 2-wide NEON registers => N = 8
        32x 8-wide AVX512 registers  => N = ?
*/

#define N 8
#define VECND vec8d
#define VECNOP(op) vec8d_##op

/********************* forward butterfly **************************************
    b0 = a0 + w*a1
    b1 = a0 - w*a1
*/

#define RADIX_2_FORWARD_PARAM_J_IS_Z(V, Q) \
    V n    = V##_set_d(Q->p); \
    V ninv = V##_set_d(Q->pinv);

#define RADIX_2_FORWARD_MOTH_J_IS_Z(V, X0, X1) \
{ \
    V x0, x1; \
    x0 = V##_load(X0); \
    x0 = V##_reduce_to_pm1n(x0, n, ninv); \
    x1 = V##_load(X1); \
    x1 = V##_reduce_to_pm1n(x1, n, ninv); \
    V##_store(X0, V##_add(x0, x1)); \
    V##_store(X1, V##_sub(x0, x1)); \
}

#define RADIX_2_FORWARD_PARAM_J_IS_NZ(V, Q, j_r, j_bits) \
    V w = V##_set_d(Q->w2tab[j_bits][j_r]); \
    V n    = V##_set_d(Q->p); \
    V ninv = V##_set_d(Q->pinv);

#define RADIX_2_FORWARD_MOTH_J_IS_NZ(V, X0, X1) \
{ \
    V x0, x1; \
    x0 = V##_load(X0); \
    x0 = V##_reduce_to_pm1n(x0, n, ninv); \
    x1 = V##_load(X1); \
    x1 = V##_mulmod(x1, w, n, ninv); \
    V##_store(X0, V##_add(x0, x1)); \
    V##_store(X1, V##_sub(x0, x1)); \
}

/* for when the V arguments above needs "evaluation" */
#define _RADIX_2_FORWARD_PARAM_J_IS_Z(...)  RADIX_2_FORWARD_PARAM_J_IS_Z(__VA_ARGS__)
#define _RADIX_2_FORWARD_MOTH_J_IS_Z(...)   RADIX_2_FORWARD_MOTH_J_IS_Z(__VA_ARGS__)
#define _RADIX_2_FORWARD_PARAM_J_IS_NZ(...) RADIX_2_FORWARD_PARAM_J_IS_NZ(__VA_ARGS__)
#define _RADIX_2_FORWARD_MOTH_J_IS_NZ(...)  RADIX_2_FORWARD_MOTH_J_IS_NZ(__VA_ARGS__)

/**************** forward butterfly with truncation **************************/

#define RADIX_2_FORWARD_MOTH_TRUNC_2_1_J_IS_Z(V, X0, X1) \
{ \
    V x0, x1; \
    x0 = V##_load(X0); \
    x0 = V##_reduce_to_pm1n(x0, n, ninv); \
    x1 = V##_load(X1); \
    x1 = V##_reduce_to_pm1n(x1, n, ninv); \
    V##_store(X0, V##_add(x0, x1)); \
}

#define RADIX_2_FORWARD_MOTH_TRUNC_2_1_J_IS_NZ(V, X0, X1) \
{ \
    V x0, x1; \
    x0 = V##_load(X0); \
    x0 = V##_reduce_to_pm1n(x0, n, ninv); \
    x1 = V##_load(X1); \
    x1 = V##_mulmod(x1, w, n, ninv); \
    V##_store(X0, V##_add(x0, x1)); \
}

#define _RADIX_2_FORWARD_MOTH_TRUNC_2_1_J_IS_Z(...)  RADIX_2_FORWARD_MOTH_TRUNC_2_1_J_IS_Z(__VA_ARGS__)
#define _RADIX_2_FORWARD_MOTH_TRUNC_2_1_J_IS_NZ(...) RADIX_2_FORWARD_MOTH_TRUNC_2_1_J_IS_NZ(__VA_ARGS__)

/********************* forward butterfly **************************************
    b0 = a0 + w^2*a2 +   w*(a1 + w^2*a3)
    b1 = a0 + w^2*a2 -   w*(a1 + w^2*a3)
    b2 = a0 - w^2*a2 + i*w*(a1 - w^2*a3)
    b3 = a0 - w^2*a2 - i*w*(a1 - w^2*a3)
*/

#define RADIX_4_FORWARD_PARAM_J_IS_Z(V, Q) \
    V iw = V##_set_d(Q->w2tab[1][0]); \
    V n    = V##_set_d(Q->p); \
    V ninv = V##_set_d(Q->pinv);

#define RADIX_4_FORWARD_MOTH_J_IS_Z(V, X0, X1, X2, X3) \
{ \
    V x0, x1, x2, x3, y0, y1, y2, y3; \
    x0 = V##_load(X0); \
    x0 = V##_reduce_to_pm1n(x0, n, ninv); \
    x1 = V##_load(X1); \
    x2 = V##_load(X2); \
    x3 = V##_load(X3); \
    x2 = V##_reduce_to_pm1n(x2, n, ninv); \
    x3 = V##_reduce_to_pm1n(x3, n, ninv); \
    y0 = V##_add(x0, x2); \
    y1 = V##_add(x1, x3); \
    y2 = V##_sub(x0, x2); \
    y3 = V##_sub(x1, x3); \
    y1 = V##_reduce_to_pm1n(y1, n, ninv); \
    y3 = V##_mulmod(y3, iw, n, ninv); \
    x0 = V##_add(y0, y1); \
    x1 = V##_sub(y0, y1); \
    x2 = V##_add(y2, y3); \
    x3 = V##_sub(y2, y3); \
    V##_store(X0, x0); \
    V##_store(X1, x1); \
    V##_store(X2, x2); \
    V##_store(X3, x3); \
}

#define RADIX_4_FORWARD_PARAM_J_IS_NZ(V, Q, j_r, j_bits) \
    FLINT_ASSERT(j_bits > 0); \
    V w  = V##_set_d(Q->w2tab[1+j_bits][2*j_r]); \
    V w2 = V##_set_d(Q->w2tab[0+j_bits][j_r]); \
    V iw = V##_set_d(Q->w2tab[1+j_bits][2*j_r+1]); \
    V n    = V##_set_d(Q->p); \
    V ninv = V##_set_d(Q->pinv);

#define RADIX_4_FORWARD_MOTH_J_IS_NZ(V, X0, X1, X2, X3) \
{ \
    V x0, x1, x2, x3, y0, y1, y2, y3; \
    x0 = V##_load(X0); \
    x0 = V##_reduce_to_pm1n(x0, n, ninv); \
    x1 = V##_load(X1); \
    x2 = V##_load(X2); \
    x3 = V##_load(X3); \
    x2 = V##_mulmod(x2, w2, n, ninv); \
    x3 = V##_mulmod(x3, w2, n, ninv); \
    y0 = V##_add(x0, x2); \
    y1 = V##_add(x1, x3); \
    y2 = V##_sub(x0, x2); \
    y3 = V##_sub(x1, x3); \
    y1 = V##_mulmod(y1, w, n, ninv); \
    y3 = V##_mulmod(y3, iw, n, ninv); \
    x0 = V##_add(y0, y1); \
    x1 = V##_sub(y0, y1); \
    x2 = V##_add(y2, y3); \
    x3 = V##_sub(y2, y3); \
    V##_store(X0, x0); \
    V##_store(X1, x1); \
    V##_store(X2, x2); \
    V##_store(X3, x3); \
}

#define _RADIX_4_FORWARD_PARAM_J_IS_Z(...)  RADIX_4_FORWARD_PARAM_J_IS_Z(__VA_ARGS__)
#define _RADIX_4_FORWARD_MOTH_J_IS_Z(...)   RADIX_4_FORWARD_MOTH_J_IS_Z(__VA_ARGS__)
#define _RADIX_4_FORWARD_PARAM_J_IS_NZ(...) RADIX_4_FORWARD_PARAM_J_IS_NZ(__VA_ARGS__)
#define _RADIX_4_FORWARD_MOTH_J_IS_NZ(...)  RADIX_4_FORWARD_MOTH_J_IS_NZ(__VA_ARGS__)

/**************** basecase transform of size BLK_SZ **************************/
/*
    The basecases below 4 are disabled because the fft is expected to be
    produced in the slightly-worse-than-bit-reversed order of basecase_4.
*/
#define DEFINE_IT(j_is_0) \
FLINT_FORCE_INLINE void CAT(sd_fft_basecase_4, j_is_0)( \
    const sd_fft_lctx_t Q, \
    double* X, \
    ulong j_r, \
    ulong j_bits) \
{ \
    vec4d n    = vec4d_set_d(Q->p); \
    vec4d ninv = vec4d_set_d(Q->pinv); \
    vec4d w, w2, iw; \
    vec4d x0, x1, x2, x3, y0, y1, y2, y3, u, v; \
 \
    /* will abuse the fact that Q->w2tab[0] points to consecutive entries */ \
    FLINT_ASSERT(SD_FFT_CTX_INIT_DEPTH >= 4); \
 \
    x0 = vec4d_load(X+0); \
    x0 = vec4d_reduce_to_pm1n(x0, n, ninv); \
    x1 = vec4d_load(X+4); \
    x2 = vec4d_load(X+8); \
    x3 = vec4d_load(X+12); \
 \
    if (j_is_0) \
    { \
        iw = vec4d_set_d(Q->w2tab[0][1]); \
 \
        x2 = vec4d_reduce_to_pm1n(x2, n, ninv); \
        x3 = vec4d_reduce_to_pm1n(x3, n, ninv); \
        y0 = vec4d_add(x0, x2); \
        y1 = vec4d_add(x1, x3); \
        y2 = vec4d_sub(x0, x2); \
        y3 = vec4d_sub(x1, x3); \
        y1 = vec4d_reduce_to_pm1n(y1, n, ninv); \
        y3 = vec4d_mulmod(y3, iw, n, ninv); \
        x0 = vec4d_add(y0, y1); \
        x1 = vec4d_sub(y0, y1); \
        x2 = vec4d_add(y2, y3); \
        x3 = vec4d_sub(y2, y3); \
    } \
    else \
    { \
        w  = vec4d_set_d(Q->w2tab[1+j_bits][2*j_r]); \
        w2 = vec4d_set_d(Q->w2tab[0+j_bits][j_r]); \
        iw = vec4d_set_d(Q->w2tab[1+j_bits][2*j_r+1]); \
 \
        x2 = vec4d_mulmod(x2, w2, n, ninv); \
        x3 = vec4d_mulmod(x3, w2, n, ninv); \
        y0 = vec4d_add(x0, x2); \
        y1 = vec4d_add(x1, x3); \
        y2 = vec4d_sub(x0, x2); \
        y3 = vec4d_sub(x1, x3); \
        y1 = vec4d_mulmod(y1, w, n, ninv); \
        y3 = vec4d_mulmod(y3, iw, n, ninv); \
        x0 = vec4d_add(y0, y1); \
        x1 = vec4d_sub(y0, y1); \
        x2 = vec4d_add(y2, y3); \
        x3 = vec4d_sub(y2, y3); \
    } \
 \
    if (j_is_0) \
    { \
        u = vec4d_load_aligned(Q->w2tab[0] + 0); \
        v = vec4d_load_aligned(Q->w2tab[0] + 4); \
        w2 = u; \
    } \
    else \
    { \
        u  = vec4d_load_aligned(Q->w2tab[3+j_bits] + 8*j_r + 0); \
        v  = vec4d_load_aligned(Q->w2tab[3+j_bits] + 8*j_r + 4); \
        w2 = vec4d_load_aligned(Q->w2tab[2+j_bits] + 4*j_r + 0); \
    } \
    w  = vec4d_unpack_lo_permute_0_2_1_3(u, v); \
    iw = vec4d_unpack_hi_permute_0_2_1_3(u, v); \
 \
    VEC4D_TRANSPOSE(x0, x1, x2, x3, x0, x1, x2, x3); \
 \
    x0 = vec4d_reduce_to_pm1n(x0, n, ninv); \
    x2 = vec4d_mulmod(x2, w2, n, ninv); \
    x3 = vec4d_mulmod(x3, w2, n, ninv); \
    y0 = vec4d_add(x0, x2); \
    y1 = vec4d_add(x1, x3); \
    y2 = vec4d_sub(x0, x2); \
    y3 = vec4d_sub(x1, x3); \
    y1 = vec4d_mulmod(y1, w, n, ninv); \
    y3 = vec4d_mulmod(y3, iw, n, ninv); \
    x0 = vec4d_add(y0, y1); \
    x1 = vec4d_sub(y0, y1); \
    x2 = vec4d_add(y2, y3); \
    x3 = vec4d_sub(y2, y3); \
 \
    /* another VEC4D_TRANSPOSE here would put the output in bit-reversed */ \
    /* but this slow down is not necessary */ \
 \
    vec4d_store(X+0, x0); \
    vec4d_store(X+4, x1); \
    vec4d_store(X+8, x2); \
    vec4d_store(X+12, x3); \
}

DEFINE_IT(0)
DEFINE_IT(1)
#undef DEFINE_IT

/* use with n = m-2 and m >= 6 */
#define EXTEND_BASECASE(n, m) \
static void CAT3(sd_fft_basecase, m, 1)(const sd_fft_lctx_t Q, double* X, ulong FLINT_UNUSED(j_r), ulong FLINT_UNUSED(j_bits)) \
{ \
    ulong l = n_pow2(m - 2); \
    _RADIX_4_FORWARD_PARAM_J_IS_Z(VECND, Q) \
    ulong i = 0; do { \
        _RADIX_4_FORWARD_MOTH_J_IS_Z(VECND, X+0*l+i, X+1*l+i, X+2*l+i, X+3*l+i) \
    } while (i += N, i < l); \
    CAT3(sd_fft_basecase, n, 1)(Q, X+0*l, 0, 0); \
    CAT3(sd_fft_basecase, n, 0)(Q, X+1*l, 0, 1); \
    CAT3(sd_fft_basecase, n, 0)(Q, X+2*l, 0, 2); \
    CAT3(sd_fft_basecase, n, 0)(Q, X+3*l, 1, 2); \
} \
static void CAT3(sd_fft_basecase, m, 0)(const sd_fft_lctx_t Q, double* X, ulong j_r, ulong j_bits) \
{ \
    ulong l = n_pow2(m - 2); \
    _RADIX_4_FORWARD_PARAM_J_IS_NZ(VECND, Q, j_r, j_bits) \
    ulong i = 0; do { \
        _RADIX_4_FORWARD_MOTH_J_IS_NZ(VECND, X+0*l+i, X+1*l+i, X+2*l+i, X+3*l+i) \
    } while (i += N, i < l); \
    CAT3(sd_fft_basecase, n, 0)(Q, X+0*l, 4*j_r+0, j_bits+2); \
    CAT3(sd_fft_basecase, n, 0)(Q, X+1*l, 4*j_r+1, j_bits+2); \
    CAT3(sd_fft_basecase, n, 0)(Q, X+2*l, 4*j_r+2, j_bits+2); \
    CAT3(sd_fft_basecase, n, 0)(Q, X+3*l, 4*j_r+3, j_bits+2); \
}
EXTEND_BASECASE(4, 6)
EXTEND_BASECASE(6, 8)
#undef EXTEND_BASECASE

/* parameter 1: j can be zero */
static void sd_fft_base_1(const sd_fft_lctx_t Q, ulong I, ulong j)
{
    ulong j_bits, j_r;
    double* x = sd_fft_lctx_blk_index(Q, I);

    FLINT_ASSERT(8 == LG_BLK_SZ);
    FLINT_ASSERT(256 == BLK_SZ);

    SET_J_BITS_AND_J_R(j_bits, j_r, j);

    if (FLINT_UNLIKELY(j == 0))
        sd_fft_basecase_8_1(Q, x, j_r, j_bits);
    else
        sd_fft_basecase_8_0(Q, x, j_r, j_bits);
}

/* parameter 0: j cannot be zero */
static void sd_fft_base_0(const sd_fft_lctx_t Q, ulong I, ulong j)
{
    ulong j_bits, j_r;
    double* x = sd_fft_lctx_blk_index(Q, I);

    FLINT_ASSERT(j != 0);
    FLINT_ASSERT(8 == LG_BLK_SZ);
    FLINT_ASSERT(256 == BLK_SZ);

    SET_J_BITS_AND_J_R(j_bits, j_r, j);

    sd_fft_basecase_8_0(Q, x, j_r, j_bits);
}


/**************** forward butterfly with truncation **************************/

/* third parameter is j == 0 */
#define DEFINE_IT(itrunc, otrunc) \
static void CAT4(sd_fft_moth_trunc_block, itrunc, otrunc, 1)( \
    const sd_fft_lctx_t Q, \
    ulong FLINT_UNUSED(j_r), ulong FLINT_UNUSED(j_bits), \
    double* X0, double* X1, double* X2, double* X3) \
{ \
    _RADIX_4_FORWARD_PARAM_J_IS_Z(VECND, Q); \
    ulong i = 0; do { \
        VECND x0, x1, x2, x3, y0, y1, y2, y3; \
        x0 = x1 = x2 = x3 = VECNOP(zero)(); \
        if (0 < itrunc) x0 = VECNOP(load)(X0+i); \
        if (0 < itrunc) x0 = VECNOP(reduce_to_pm1n)(x0, n, ninv); \
        if (1 < itrunc) x1 = VECNOP(load)(X1+i); \
        if (2 < itrunc) x2 = VECNOP(load)(X2+i); \
        if (2 < itrunc) x2 = VECNOP(reduce_to_pm1n)(x2, n, ninv); \
        if (3 < itrunc) x3 = VECNOP(load)(X3+i); \
        if (3 < itrunc) x3 = VECNOP(reduce_to_pm1n)(x3, n, ninv); \
        y0 = (2 < itrunc) ? VECNOP(add)(x0, x2) : x0; \
        y1 = (3 < itrunc) ? VECNOP(add)(x1, x3) : x1; \
        y2 = (2 < itrunc) ? VECNOP(sub)(x0, x2) : x0; \
        y3 = (3 < itrunc) ? VECNOP(sub)(x1, x3) : x1; \
        y1 = VECNOP(reduce_to_pm1n)(y1, n, ninv); \
        y3 = VECNOP(mulmod)(y3, iw, n, ninv); \
        x0 = VECNOP(add)(y0, y1); \
        x1 = VECNOP(sub)(y0, y1); \
        x2 = VECNOP(add)(y2, y3); \
        x3 = VECNOP(sub)(y2, y3); \
        if (0 < otrunc) VECNOP(store)(X0+i, x0); \
        if (1 < otrunc) VECNOP(store)(X1+i, x1); \
        if (2 < otrunc) VECNOP(store)(X2+i, x2); \
        if (3 < otrunc) VECNOP(store)(X3+i, x3); \
    } while (i += N, i < BLK_SZ); \
    FLINT_ASSERT(i == BLK_SZ); \
} \
static void CAT4(sd_fft_moth_trunc_block, itrunc, otrunc, 0)( \
    const sd_fft_lctx_t Q, \
    ulong j_r, ulong j_bits, \
    double* X0, double* X1, double* X2, double* X3) \
{ \
    _RADIX_4_FORWARD_PARAM_J_IS_NZ(VECND, Q, j_r, j_bits); \
    ulong i = 0; do { \
        VECND x0, x1, x2, x3, y0, y1, y2, y3; \
        x0 = x1 = x2 = x3 = VECNOP(zero)(); \
        if (0 < itrunc) x0 = VECNOP(load)(X0+i); \
        if (0 < itrunc) x0 = VECNOP(reduce_to_pm1n)(x0, n, ninv); \
        if (1 < itrunc) x1 = VECNOP(load)(X1+i); \
        if (2 < itrunc) x2 = VECNOP(load)(X2+i); \
        if (2 < itrunc) x2 = VECNOP(mulmod)(x2, w2, n, ninv); \
        if (3 < itrunc) x3 = VECNOP(load)(X3+i); \
        if (3 < itrunc) x3 = VECNOP(mulmod)(x3, w2, n, ninv); \
        y0 = (2 < itrunc) ? VECNOP(add)(x0, x2) : x0; \
        y1 = (3 < itrunc) ? VECNOP(add)(x1, x3) : x1; \
        y2 = (2 < itrunc) ? VECNOP(sub)(x0, x2) : x0; \
        y3 = (3 < itrunc) ? VECNOP(sub)(x1, x3) : x1; \
        y1 = VECNOP(mulmod)(y1, w, n, ninv); \
        y3 = VECNOP(mulmod)(y3, iw, n, ninv); \
        x0 = VECNOP(add)(y0, y1); \
        x1 = VECNOP(sub)(y0, y1); \
        x2 = VECNOP(add)(y2, y3); \
        x3 = VECNOP(sub)(y2, y3); \
        if (0 < otrunc) VECNOP(store)(X0+i, x0); \
        if (1 < otrunc) VECNOP(store)(X1+i, x1); \
        if (2 < otrunc) VECNOP(store)(X2+i, x2); \
        if (3 < otrunc) VECNOP(store)(X3+i, x3); \
    } while (i += N, i < BLK_SZ); \
    FLINT_ASSERT(i == BLK_SZ); \
}

DEFINE_IT(2, 1)
DEFINE_IT(2, 2)
DEFINE_IT(2, 3)
DEFINE_IT(2, 4)
DEFINE_IT(3, 1)
DEFINE_IT(3, 2)
DEFINE_IT(3, 3)
DEFINE_IT(3, 4)
DEFINE_IT(4, 1)
DEFINE_IT(4, 2)
DEFINE_IT(4, 3)
DEFINE_IT(4, 4)
#undef DEFINE_IT

/************************ the recursive stuff ********************************/

void sd_fft_main_block(
    const sd_fft_lctx_t Q,
    ulong I, /* starting index */
    ulong S, /* stride */
    ulong k, /* BLK_SZ transforms each of length 2^k */
    ulong j)
{
    ulong j_bits, j_r;

    if (k > 4)
    {
        ulong k1 = k/2;
        ulong k2 = k - k1;

        ulong l2 = n_pow2(k2);
        ulong a = 0; do {
            sd_fft_main_block(Q, I + a*S, S<<k2, k1, j);
        } while (a++, a < l2);

        // row ffts
        ulong l1 = n_pow2(k1);
        ulong b = 0; do {
            sd_fft_main_block(Q, I + (b<<k2)*S, S, k2, (j<<k1) + b);
        } while (b++, b < l1);

        return;
    }

    SET_J_BITS_AND_J_R(j_bits, j_r, j);

    if (k >= 2)
    {
        ulong k1 = 2;
        ulong k2 = k - k1;
        ulong l2 = n_pow2(k2);

        /* column ffts */
        if (FLINT_UNLIKELY(j_bits == 0))
        {
            _RADIX_4_FORWARD_PARAM_J_IS_Z(VECND, Q)
            ulong a = 0; do {
                double* X0 = sd_fft_lctx_blk_index(Q, I+a*S + (S<<k2)*0);
                double* X1 = sd_fft_lctx_blk_index(Q, I+a*S + (S<<k2)*1);
                double* X2 = sd_fft_lctx_blk_index(Q, I+a*S + (S<<k2)*2);
                double* X3 = sd_fft_lctx_blk_index(Q, I+a*S + (S<<k2)*3);
                ulong i = 0; do {
                    _RADIX_4_FORWARD_MOTH_J_IS_Z(VECND, X0+i, X1+i, X2+i, X3+i);
                } while (i += N, i < BLK_SZ);
            } while (a++, a < l2);
        }
        else
        {
            _RADIX_4_FORWARD_PARAM_J_IS_NZ(VECND, Q, j_r, j_bits)
            ulong a = 0; do {
                double* X0 = sd_fft_lctx_blk_index(Q, I+a*S + (S<<k2)*0);
                double* X1 = sd_fft_lctx_blk_index(Q, I+a*S + (S<<k2)*1);
                double* X2 = sd_fft_lctx_blk_index(Q, I+a*S + (S<<k2)*2);
                double* X3 = sd_fft_lctx_blk_index(Q, I+a*S + (S<<k2)*3);
                ulong i = 0; do {
                    _RADIX_4_FORWARD_MOTH_J_IS_NZ(VECND, X0+i, X1+i, X2+i, X3+i);
                } while (i += N, i < BLK_SZ);
            } while (a++, a < l2);
        }

        if (l2 == 1)
            return;

        /* row ffts */
        ulong l1 = n_pow2(k1);
        ulong b = 0; do {
            sd_fft_main_block(Q, I + (b<<k2)*S, S, k2, (j<<k1) + b);
        } while (b++, b < l1);
    }
    else if (k == 1)
    {
        double* X0 = sd_fft_lctx_blk_index(Q, I + S*0);
        double* X1 = sd_fft_lctx_blk_index(Q, I + S*1);
        if (FLINT_UNLIKELY(j_bits == 0))
        {
            _RADIX_2_FORWARD_PARAM_J_IS_Z(VECND, Q)
            ulong i = 0; do {
                _RADIX_2_FORWARD_MOTH_J_IS_Z(VECND, X0+i, X1+i);
            } while (i += N, i < BLK_SZ);
        }
        else
        {
            _RADIX_2_FORWARD_PARAM_J_IS_NZ(VECND, Q, j_r, j_bits)
            ulong i = 0; do {
                _RADIX_2_FORWARD_MOTH_J_IS_NZ(VECND, X0+i, X1+i);
            } while (i += N, i < BLK_SZ);
        }
    }
}



void sd_fft_main(
    const sd_fft_lctx_t Q,
    ulong I,    /* starting index */
    ulong S,    /* stride */
    ulong k,    /* 1 transform of length BLK_SZ*2^k */
    ulong j)    /* twist param */
{
    if (k > 2)
    {
        ulong k1 = k/2;
        ulong k2 = k - k1;

        /* column ffts */
        ulong l2 = n_pow2(k2);
        ulong a = 0; do {
            sd_fft_main_block(Q, I + a*S, S<<k2, k1, j);
        } while (a++, a < l2);

        /* row ffts */
        ulong l1 = n_pow2(k1);
        ulong b = 0; do {
            sd_fft_main(Q, I + b*(S<<k2), S, k2, (j<<k1) + b);
        } while (b++, b < l1);

        return;
    }

    if (k == 2)
    {
        // k1 = 2; k2 = 0
        sd_fft_main_block(Q, I, S, 2, j);
        sd_fft_base_1(Q, I + S*0, 4*j + 0);
        sd_fft_base_0(Q, I + S*1, 4*j + 1);
        sd_fft_base_0(Q, I + S*2, 4*j + 2);
        sd_fft_base_0(Q, I + S*3, 4*j + 3);
    }
    else if (k == 1)
    {
        // k1 = 1; k2 = 0
        sd_fft_main_block(Q, I, S, 1, j);
        sd_fft_base_1(Q, I + S*0, 2*j + 0);
        sd_fft_base_0(Q, I + S*1, 2*j + 1);
    }
    else
    {
        sd_fft_base_1(Q, I, j);
    }
}


void sd_fft_trunc_block(
    const sd_fft_lctx_t Q,
    ulong I, // starting index
    ulong S, // stride
    ulong k, // transform length 2^k
    ulong j,
    ulong itrunc,
    ulong otrunc)
{
    ulong j_bits, j_r;

    FLINT_ASSERT(itrunc <= n_pow2(k));
    FLINT_ASSERT(otrunc <= n_pow2(k));

    if (otrunc < 1)
        return;

    if (itrunc <= 1)
    {
        if (itrunc < 1)
        {
            for (ulong a = 0; a < otrunc; a++)
            {
                double* X0 = sd_fft_lctx_blk_index(Q, I + S*a);
                VECND z = VECNOP(zero)();
                ulong i = 0; do {
                    VECNOP(store)(X0+i, z);
                } while (i += N, i < BLK_SZ);
            }
        }
        else
        {
            double* X0 = sd_fft_lctx_blk_index(Q, I + S*0);
            for (ulong a = 1; a < otrunc; a++)
            {
                double* X1 = sd_fft_lctx_blk_index(Q, I + S*a);
                ulong i = 0; do {
                    VECND u = VECNOP(load)(X0+i);
                    VECNOP(store)(X1+i, u);
                } while (i += N, i < BLK_SZ);
            }
        }

        return;
    }

    if (itrunc == otrunc && otrunc == n_pow2(k))
    {
        sd_fft_main_block(Q, I, S, k, j);
        return;
    }

    if (k > 2)
    {
        ulong k1 = k/2;
        ulong k2 = k - k1;

        ulong l2 = n_pow2(k2);
        ulong n1 = otrunc >> k2;
        ulong n2 = otrunc & (l2 - 1);
        ulong z1 = itrunc >> k2;
        ulong z2 = itrunc & (l2 - 1);
        ulong n1p = n1 + (n2 != 0);
        ulong z2p = n_min(l2, itrunc);

        /* columns */
        for (ulong a = 0; a < z2p; a++)
            sd_fft_trunc_block(Q, I + a*S, S << k2, k1, j, z1 + (a < z2), n1p);

        /* full rows */
        for (ulong b = 0; b < n1; b++)
            sd_fft_trunc_block(Q, I + b*(S << k2), S, k2, (j << k1) + b, z2p, l2);

        /* last partial row */
        if (n2 > 0)
            sd_fft_trunc_block(Q, I + n1*(S << k2), S, k2, (j << k1) + n1, z2p, n2);

        return;
    }

    SET_J_BITS_AND_J_R(j_bits, j_r, j);

    if (k == 2)
    {
#define IT(ii, oo) CAT4(sd_fft_moth_trunc_block, ii, oo, 0), \
                   CAT4(sd_fft_moth_trunc_block, ii, oo, 1)
#define LOOKUP_IT(ii, oo, j_is_zero) tab[(j_is_zero) + 2*((oo)-1 + 4*((ii)-2))]
        static void (*tab[3*4*2])(const sd_fft_lctx_t, ulong, ulong, double*, double*, double*, double*) =
                        {IT(2,1), IT(2,2), IT(2,3), IT(2,4),
                         IT(3,1), IT(3,2), IT(3,3), IT(3,4),
                         IT(4,1), IT(4,2), IT(4,3), IT(4,4)};

        LOOKUP_IT(itrunc, otrunc, j == 0)(Q, j_r, j_bits,
                                            sd_fft_lctx_blk_index(Q, I + S*0),
                                            sd_fft_lctx_blk_index(Q, I + S*1),
                                            sd_fft_lctx_blk_index(Q, I + S*2),
                                            sd_fft_lctx_blk_index(Q, I + S*3));
#undef LOOKUP_IT
#undef IT
    }
    else if (k == 1)
    {
        double* X0 = sd_fft_lctx_blk_index(Q, I + S*0);
        double* X1 = sd_fft_lctx_blk_index(Q, I + S*1);
        FLINT_ASSERT(itrunc == 2);
        FLINT_ASSERT(otrunc == 1);
        if (FLINT_UNLIKELY(j_bits == 0))
        {
            _RADIX_2_FORWARD_PARAM_J_IS_Z(VECND, Q)
            ulong i = 0; do {
                _RADIX_2_FORWARD_MOTH_TRUNC_2_1_J_IS_Z(VECND, X0 + i, X1 + i);
            } while (i += N, i < BLK_SZ);
        }
        else
        {
            _RADIX_2_FORWARD_PARAM_J_IS_NZ(VECND, Q, j_r, j_bits)
            ulong i = 0; do {
                _RADIX_2_FORWARD_MOTH_TRUNC_2_1_J_IS_NZ(VECND, X0 + i, X1 + i);
            } while (i += N, i < BLK_SZ);
        }
    }
}


void sd_fft_trunc(
    const sd_fft_lctx_t Q,
    ulong I, /* starting index */
    ulong S, /* stride */
    ulong k, /* transform length BLK_SZ*2^k */
    ulong j,
    ulong itrunc,   /* actual trunc is BLK_SZ*itrunc */
    ulong otrunc)   /* actual trunc is BLK_SZ*otrunc */
{
    if (otrunc < 1)
        return;

    if (itrunc < 1)
    {
        for (ulong a = 0; a < otrunc; a++)
        {
            double* X0 = sd_fft_lctx_blk_index(Q, I + S*a);
            VECND z = VECNOP(zero)();
            ulong i = 0; do {
                VECNOP(store)(X0 + i, z);
            } while (i += N, i < BLK_SZ);
        }

        return;
    }

    if (k > 2)
    {
        ulong k1 = k/2;
        ulong k2 = k - k1;

        ulong l2 = n_pow2(k2);
        ulong n1 = otrunc >> k2;
        ulong n2 = otrunc & (l2 - 1);
        ulong z1 = itrunc >> k2;
        ulong z2 = itrunc & (l2 - 1);
        ulong n1p = n1 + (n2 != 0);
        ulong z2p = n_min(l2, itrunc);

        /* columns */
        for (ulong a = 0; a < z2p; a++)
            sd_fft_trunc_block(Q, I + a*S, S << k2, k1, j, z1 + (a < z2), n1p);

        /* full rows */
        for (ulong b = 0; b < n1; b++)
            sd_fft_trunc(Q, I + b*(S << k2), S, k2, (j << k1) + b, z2p, l2);

        /* last partial row */
        if (n2 > 0)
            sd_fft_trunc(Q, I + n1*(S << k2), S, k2, (j << k1) + n1, z2p, n2);

        return;
    }

    if (k == 2)
    {
        sd_fft_trunc_block(Q, I, S, 2, j, itrunc, otrunc);
                        sd_fft_base_1(Q, I + S*0, 4*j + 0);
        if (otrunc > 1) sd_fft_base_0(Q, I + S*1, 4*j + 1);
        if (otrunc > 2) sd_fft_base_0(Q, I + S*2, 4*j + 2);
        if (otrunc > 3) sd_fft_base_0(Q, I + S*3, 4*j + 3);
    }
    else if (k == 1)
    {
        sd_fft_trunc_block(Q, I, S, 1, j, itrunc, otrunc);
                        sd_fft_base_1(Q, I + S*0, 2*j + 0);
        if (otrunc > 1) sd_fft_base_0(Q, I + S*1, 2*j + 1);
    }
    else
    {
        sd_fft_base_1(Q, I, j);
    }
}




#undef RADIX_2_FORWARD_PARAM_J_IS_Z
#undef RADIX_2_FORWARD_MOTH_J_IS_Z
#undef RADIX_2_FORWARD_PARAM_J_IS_NZ
#undef RADIX_2_FORWARD_MOTH_J_IS_NZ
#undef RADIX_2_FORWARD_MOTH_TRUNC_2_1_J_IS_Z
#undef RADIX_2_FORWARD_MOTH_TRUNC_2_1_J_IS_NZ
#undef RADIX_4_FORWARD_PARAM_J_IS_Z
#undef RADIX_4_FORWARD_MOTH_J_IS_Z
#undef RADIX_4_FORWARD_PARAM_J_IS_NZ
#undef RADIX_4_FORWARD_MOTH_J_IS_NZ


#undef _RADIX_2_FORWARD_PARAM_J_IS_Z
#undef _RADIX_2_FORWARD_MOTH_J_IS_Z
#undef _RADIX_2_FORWARD_PARAM_J_IS_NZ
#undef _RADIX_2_FORWARD_MOTH_J_IS_NZ
#undef _RADIX_2_FORWARD_MOTH_TRUNC_2_1_J_IS_Z
#undef _RADIX_2_FORWARD_MOTH_TRUNC_2_1_J_IS_NZ
#undef _RADIX_4_FORWARD_PARAM_J_IS_Z
#undef _RADIX_4_FORWARD_MOTH_J_IS_Z
#undef _RADIX_4_FORWARD_PARAM_J_IS_NZ
#undef _RADIX_4_FORWARD_MOTH_J_IS_NZ

#undef N
#undef VECND
#undef VECNOP
