/*
    Copyright (C) 2022 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fft_small.h"

/*
    N is supposed to be a good fit for the number of points to process per loop
    in the radix 4 inverse butterflies.

    M is supposed to be <= N and is used for truncation formulas that have
    a larger register pressure.
*/

#define N 8
#define VECND vec8d
#define VECNOP(op) vec8d_##op

#define M 4
#define VECMD vec4d
#define VECMOP(op) vec4d_##op


/************************** inverse butterfly *********************************
    2*a0 =      (b0 + b1)
    2*a1 = w^-1*(b0 - b1)
    W := -w^-1
*/
#define RADIX_2_REVERSE_PARAM_J_IS_Z(V, Q) \
    V n    = V##_set_d(Q->p); \
    V ninv = V##_set_d(Q->pinv);

#define RADIX_2_REVERSE_MOTH_J_IS_Z(V, X0, X1) \
{ \
    V x0, x1, y0, y1; \
    x0 = V##_load(X0); \
    x1 = V##_load(X1); \
    y0 = V##_add(x0, x1); \
    y1 = V##_sub(x0, x1); \
    y0 = V##_reduce_to_pm1n(y0, n, ninv); \
    y1 = V##_reduce_to_pm1n(y1, n, ninv); \
    V##_store(X0, y0); \
    V##_store(X1, y1); \
}

#define RADIX_2_REVERSE_PARAM_J_IS_NZ(V, Q, j_mr, j_bits) \
    V W = V##_set_d(Q->w2tab[j_bits][j_mr]); \
    V n    = V##_set_d(Q->p); \
    V ninv = V##_set_d(Q->pinv);

#define RADIX_2_REVERSE_MOTH_J_IS_NZ(V, X0, X1) \
{ \
    V x0, x1, y0, y1; \
    x0 = V##_load(X0); \
    x1 = V##_load(X1); \
    y0 = V##_add(x0, x1); \
    y1 = V##_sub(x1, x0); \
    y0 = V##_reduce_to_pm1n(y0, n, ninv); \
    y1 = V##_mulmod(y1, W, n, ninv); \
    V##_store(X0, y0); \
    V##_store(X1, y1); \
}

#define _RADIX_2_REVERSE_PARAM_J_IS_Z(...)  RADIX_2_REVERSE_PARAM_J_IS_Z(__VA_ARGS__)
#define _RADIX_2_REVERSE_MOTH_J_IS_Z(...)   RADIX_2_REVERSE_MOTH_J_IS_Z(__VA_ARGS__)
#define _RADIX_2_REVERSE_PARAM_J_IS_NZ(...) RADIX_2_REVERSE_PARAM_J_IS_NZ(__VA_ARGS__)
#define _RADIX_2_REVERSE_MOTH_J_IS_NZ(...)  RADIX_2_REVERSE_MOTH_J_IS_NZ(__VA_ARGS__)

/****************** inverse butterfly with truncation ************************/

/* the legal functions must be implemented */
#define DEFINE_IT(nn, zz, ff) \
static void CAT4(radix_2_moth_inv_trunc_block, nn, zz, ff)( \
    const sd_fft_lctx_t Q, \
    ulong j, \
    double* X0, double* X1) \
{ \
    int l = 2; \
    flint_printf("function l = %d, n = %d, z = %d, f = %d", l, nn, zz, ff); \
    if (1 <= zz && zz <= 2 && nn <= zz && 1 <= nn+ff && nn+ff <= l) \
        flint_throw(FLINT_ERROR, " is not implemented\n"); \
    else \
        flint_throw(FLINT_ERROR, " does not exist and should not be called\n"); \
}

DEFINE_IT(0,1,0)
DEFINE_IT(0,2,0)
DEFINE_IT(2,1,0)
DEFINE_IT(2,1,1)
DEFINE_IT(2,2,0)
DEFINE_IT(2,2,1)
#undef DEFINE_IT


/* {x0, x1} = {2*x0 - w*x1, x0 - w*x1} */
static void radix_2_moth_inv_trunc_block_1_2_1(
    const sd_fft_lctx_t Q,
    ulong j,
    double* X0, double* X1)
{
    VECND n    = VECNOP(set_d)(Q->p);
    VECND ninv = VECNOP(set_d)(Q->pinv);
    VECND w = VECNOP(set_d)(sd_fft_lctx_w2(Q, j));
    VECND c = VECNOP(set_d)(2);
    ulong i = 0; do {
        VECND a, b, u, v;
        a = VECNOP(load)(X0 + i);
        b = VECNOP(load)(X1 + i);
        b = VECNOP(nmulmod)(b, w, n, ninv);
        u = VECNOP(fmadd)(c, a, b);
        v = VECNOP(add)(a, b);
        u = VECNOP(reduce_to_pm1n)(u, n, ninv);
        v = VECNOP(reduce_to_pm1n)(v, n, ninv);
        VECNOP(store)(X0 + i, u);
        VECNOP(store)(X1 + i, v);
    } while (i += N, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/* {x0} = {2*x0 - w*x1} */
static void radix_2_moth_inv_trunc_block_1_2_0(
    const sd_fft_lctx_t Q,
    ulong j,
    double* X0, double* X1)
{
    VECND n    = VECNOP(set_d)(Q->p);
    VECND ninv = VECNOP(set_d)(Q->pinv);
    VECND w = VECNOP(set_d)(sd_fft_lctx_w2(Q, j));
    VECND c = VECNOP(set_d)(2);
    ulong i = 0; do {
        VECND a, b, u;
        a = VECNOP(load)(X0 + i);
        b = VECNOP(load)(X1 + i);
        b = VECNOP(mulmod)(b, w, n, ninv);
        u = VECNOP(fmsub)(c, a, b);
        u = VECNOP(reduce_to_pm1n)(u, n, ninv);
        VECNOP(store)(X0 + i, u);
    } while (i += N, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/* {x0, x1} = {2*x0, x0} */
static void radix_2_moth_inv_trunc_block_1_1_1(
    const sd_fft_lctx_t Q,
    ulong j,
    double* X0, double* X1)
{
    VECND n    = VECNOP(set_d)(Q->p);
    VECND ninv = VECNOP(set_d)(Q->pinv);
    ulong i = 0; do {
        VECND a, u;
        a = VECNOP(load)(X0 + i);
        u = VECNOP(add)(a, a);
        u = VECNOP(reduce_to_pm1n)(u, n, ninv);
        VECNOP(store)(X0 + i, u);
        VECNOP(store)(X1 + i, a);
    } while (i += N, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/* {x0} = {2*x0} */
static void radix_2_moth_inv_trunc_block_1_1_0(
    const sd_fft_lctx_t Q,
    ulong j,
    double* X0, double* X1)
{
    VECND n    = VECNOP(set_d)(Q->p);
    VECND ninv = VECNOP(set_d)(Q->pinv);
    ulong i = 0; do {
        VECND a;
        a = VECNOP(load)(X0 + i);
        a = VECNOP(add)(a, a);
        a = VECNOP(reduce_to_pm1n)(a, n, ninv);
        VECNOP(store)(X0 + i, a);
    } while (i += N, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/* {x0} = {(x0 + w*x1)/2} */
static void radix_2_moth_inv_trunc_block_0_2_1(
    const sd_fft_lctx_t Q,
    ulong j,
    double* X0, double* X1)
{
    VECND n    = VECNOP(set_d)(Q->p);
    VECND ninv = VECNOP(set_d)(Q->pinv);
    VECND w = VECNOP(set_d)(sd_fft_lctx_w2(Q, j));
    VECND c = VECNOP(set_d)(vec1d_fnmadd(0.5, Q->p, 0.5));
    ulong i = 0; do {
        VECND a, b;
        a = VECNOP(load)(X0 + i);
        b = VECNOP(load)(X1 + i);
        b = VECNOP(mulmod)(b, w, n, ninv);
        a = VECNOP(add)(a, b);
        a = VECNOP(mulmod)(a, c, n, ninv);
        VECNOP(store)(X0 + i, a);
    } while (i += N, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/* {x0} = {(x0)/2} */
static void radix_2_moth_inv_trunc_block_0_1_1(
    const sd_fft_lctx_t Q,
    ulong j,
    double* X0, double* X1)
{
    VECND n    = VECNOP(set_d)(Q->p);
    VECND ninv = VECNOP(set_d)(Q->pinv);
    VECND c = VECNOP(set_d)(vec1d_fnmadd(0.5, Q->p, 0.5));
    ulong i = 0; do {
        VECND a;
        a = VECNOP(load)(X0 + i);
        a = VECNOP(mulmod)(a, c, n, ninv);
        VECNOP(store)(X0 + i, a);
    } while (i += N, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/************************* inverse butterfly **********************************
    4*a0 =            (b0 + b1) +        (b2 + b3)
    4*a1 =       w^-1*(b0 - b1) - i*w^-1*(b2 - b3)
    4*a2 = w^-2*(     (b0 + b1) -        (b2 + b3))
    4*a3 = w^-2*(w^-1*(b0 - b1) + i*w^-1*(b2 - b3))
    W  := -w^-1
    W2 := -w^-2
    IW := i*w^-1
*/
#define RADIX_4_REVERSE_PARAM_J_IS_Z(V, Q) \
    V IW = V##_set_d(Q->w2tab[0][1]); \
    V n    = V##_set_d(Q->p); \
    V ninv = V##_set_d(Q->pinv); \

#define RADIX_4_REVERSE_MOTH_J_IS_Z(V, X0, X1, X2, X3) \
{ \
    V x0, x1, x2, x3, y0, y1, y2, y3; \
    x0 = V##_load(X0); \
    x1 = V##_load(X1); \
    x2 = V##_load(X2); \
    x3 = V##_load(X3); \
    y0 = V##_add(x0, x1); \
    y1 = V##_add(x2, x3); \
    y2 = V##_sub(x0, x1); \
    y3 = V##_sub(x2, x3); \
    y0 = V##_reduce_to_pm1n(y0, n, ninv); \
    y1 = V##_reduce_to_pm1n(y1, n, ninv); \
    y2 = V##_reduce_to_pm1n(y2, n, ninv); \
    y3 = V##_mulmod(y3, IW, n, ninv); \
    V##_store(X0, V##_add(y0, y1)); \
    V##_store(X2, V##_sub(y0, y1)); \
    V##_store(X1, V##_sub(y2, y3)); \
    V##_store(X3, V##_add(y2, y3)); \
}

#define RADIX_4_REVERSE_PARAM_J_IS_NZ(V, Q, j_mr, j_bits) \
    V W  = V##_set_d(Q->w2tab[1+j_bits][2*j_mr+1]); \
    V W2 = V##_set_d(Q->w2tab[0+j_bits][j_mr]); \
    V IW = V##_set_d(Q->w2tab[1+j_bits][2*j_mr+0]); \
    V n    = V##_set_d(Q->p); \
    V ninv = V##_set_d(Q->pinv); \

#define RADIX_4_REVERSE_MOTH_J_IS_NZ(V, X0, X1, X2, X3) \
{ \
    V x0, x1, x2, x3, y0, y1, y2, y3; \
    x0 = V##_load(X0); \
    x1 = V##_load(X1); \
    x2 = V##_load(X2); \
    x3 = V##_load(X3); \
    y0 = V##_add(x0, x1); \
    y1 = V##_add(x2, x3); \
    y2 = V##_sub(x0, x1); \
    y3 = V##_sub(x3, x2); \
    y2 = V##_mulmod(y2, W, n, ninv); \
    y3 = V##_mulmod(y3, IW, n, ninv); \
    x0 = V##_add(y0, y1); \
    x1 = V##_sub(y3, y2); \
    V##_store(X1, x1); \
    x2 = V##_sub(y1, y0); \
    x3 = V##_add(y3, y2); \
    x0 = V##_reduce_to_pm1n(x0, n, ninv); \
    x2 = V##_mulmod(x2, W2, n, ninv); \
    x3 = V##_mulmod(x3, W2, n, ninv); \
    V##_store(X0, x0); \
    V##_store(X2, x2); \
    V##_store(X3, x3); \
}

#define _RADIX_4_REVERSE_PARAM_J_IS_Z(...)  RADIX_4_REVERSE_PARAM_J_IS_Z(__VA_ARGS__)
#define _RADIX_4_REVERSE_MOTH_J_IS_Z(...)   RADIX_4_REVERSE_MOTH_J_IS_Z(__VA_ARGS__)
#define _RADIX_4_REVERSE_PARAM_J_IS_NZ(...) RADIX_4_REVERSE_PARAM_J_IS_NZ(__VA_ARGS__)
#define _RADIX_4_REVERSE_MOTH_J_IS_NZ(...)  RADIX_4_REVERSE_MOTH_J_IS_NZ(__VA_ARGS__)

/************ basecase inverse transform of size BLK_SZ **********************/
/*
    The basecases below 4 are disabled because the ifft expects input in
    slightly-worse-than-bit-reversed order as in basecase_4.
*/
#define DEFINE_IT(j_is_0) \
FLINT_FORCE_INLINE void sd_ifft_basecase_4_##j_is_0(\
    const sd_fft_lctx_t Q, double* X, ulong j_mr, ulong j_bits) \
{ \
    vec4d n    = vec4d_set_d(Q->p); \
    vec4d ninv = vec4d_set_d(Q->pinv); \
    vec4d W, W2, IW, u, v; \
    vec4d x0, x1, x2, x3, y0, y1, y2, y3; \
 \
    x0 = vec4d_load(X + 0); \
    x1 = vec4d_load(X + 4); \
    x2 = vec4d_load(X + 8); \
    x3 = vec4d_load(X + 12); \
 \
    if (j_is_0) \
    { \
        W  = vec4d_set_d4(-Q->w2tab[0][0], Q->w2tab[0][3], Q->w2tab[0][7], Q->w2tab[0][5]); \
        IW = vec4d_set_d4( Q->w2tab[0][1], Q->w2tab[0][2], Q->w2tab[0][6], Q->w2tab[0][4]); \
        W2 = vec4d_set_d4(-Q->w2tab[0][0], Q->w2tab[0][1], Q->w2tab[0][3], Q->w2tab[0][2]); \
    } \
    else \
    { \
        W2 = vec4d_load_aligned(Q->w2tab[2+j_bits] + 4*j_mr);   /* a b c d */ \
        W2 = vec4d_permute_3_2_1_0(W2);                         /* d c b a */ \
        u = vec4d_load_aligned(Q->w2tab[3+j_bits] + 8*j_mr+0);  /* 0 1 2 3 */ \
        v = vec4d_load_aligned(Q->w2tab[3+j_bits] + 8*j_mr+4);  /* 4 5 6 7 */ \
        W  = vec4d_unpackhi_permute_3_1_2_0(u, v);              /* 7 5 3 1 */ \
        IW = vec4d_unpacklo_permute_3_1_2_0(u, v);              /* 6 4 2 0 */ \
    } \
 \
    y0 = vec4d_add(x0, x1); \
    y1 = vec4d_add(x2, x3); \
    y2 = vec4d_sub(x0, x1); \
    y3 = vec4d_sub(x3, x2); \
    y2 = vec4d_mulmod(y2, W, n, ninv); \
    y3 = vec4d_mulmod(y3, IW, n, ninv); \
    x0 = vec4d_add(y0, y1); \
    x1 = vec4d_sub(y3, y2); \
    x2 = vec4d_sub(y1, y0); \
    x3 = vec4d_add(y3, y2); \
    x0 = vec4d_reduce_to_pm1n(x0, n, ninv); \
    x2 = vec4d_mulmod(x2, W2, n, ninv); \
    x3 = vec4d_mulmod(x3, W2, n, ninv); \
 \
    VEC4D_TRANSPOSE(x0, x1, x2, x3, x0, x1, x2, x3); \
 \
    if (j_is_0) \
    { \
        IW = vec4d_set_d(Q->w2tab[0][1]); \
        y0 = vec4d_add(x0, x1); \
        y1 = vec4d_add(x2, x3); \
        y2 = vec4d_sub(x0, x1); \
        y3 = vec4d_sub(x2, x3); \
        y0 = vec4d_reduce_to_pm1n(y0, n, ninv); \
        y1 = vec4d_reduce_to_pm1n(y1, n, ninv); \
        y2 = vec4d_reduce_to_pm1n(y2, n, ninv); \
        y3 = vec4d_mulmod(y3, IW, n, ninv); \
        vec4d_store(X+0, vec4d_add(y0, y1)); \
        vec4d_store(X+8, vec4d_sub(y0, y1)); \
        vec4d_store(X+4, vec4d_sub(y2, y3)); \
        vec4d_store(X+12, vec4d_add(y2, y3)); \
    } \
    else \
    { \
        W  = vec4d_set_d(Q->w2tab[1+j_bits][2*j_mr+1]); \
        IW = vec4d_set_d(Q->w2tab[1+j_bits][2*j_mr+0]); \
        W2 = vec4d_set_d(Q->w2tab[0+j_bits][j_mr]); \
        y0 = vec4d_add(x0, x1); \
        y1 = vec4d_add(x2, x3); \
        y2 = vec4d_sub(x0, x1); \
        y3 = vec4d_sub(x3, x2); \
        y2 = vec4d_mulmod(y2, W, n, ninv); \
        y3 = vec4d_mulmod(y3, IW, n, ninv); \
        x0 = vec4d_add(y0, y1); \
        x1 = vec4d_sub(y3, y2); \
        x2 = vec4d_sub(y1, y0); \
        x3 = vec4d_add(y3, y2); \
        x0 = vec4d_reduce_to_pm1n(x0, n, ninv); \
        x2 = vec4d_mulmod(x2, W2, n, ninv); \
        x3 = vec4d_mulmod(x3, W2, n, ninv); \
        vec4d_store(X+0, x0); \
        vec4d_store(X+4, x1); \
        vec4d_store(X+8, x2); \
        vec4d_store(X+12, x3); \
    } \
}

DEFINE_IT(0)
DEFINE_IT(1)
#undef DEFINE_IT


/* use with n = m-2 and m >= 6 */
#define EXTEND_BASECASE(n, m) \
void CAT3(sd_ifft_basecase, m, 1)(const sd_fft_lctx_t Q, double* X, ulong j_mr, ulong j_bits) \
{ \
    ulong l = n_pow2(m - 2); \
    FLINT_ASSERT(j_bits == 0); \
    CAT3(sd_ifft_basecase, n, 1)(Q, X+0*l, 0, 0); \
    CAT3(sd_ifft_basecase, n, 0)(Q, X+1*l, 0, 1); \
    CAT3(sd_ifft_basecase, n, 0)(Q, X+2*l, 1, 2); \
    CAT3(sd_ifft_basecase, n, 0)(Q, X+3*l, 0, 2); \
    { \
        _RADIX_4_REVERSE_PARAM_J_IS_Z(VECND, Q) \
        ulong i = 0; do { \
            _RADIX_4_REVERSE_MOTH_J_IS_Z(VECND, X+0*l+i, X+1*l+i, X+2*l+i, X+3*l+i) \
        } while (i += N, i < l); \
        FLINT_ASSERT(i == l); \
    } \
} \
void CAT3(sd_ifft_basecase, m, 0)(const sd_fft_lctx_t Q, double* X, ulong j_mr, ulong j_bits) \
{ \
    ulong l = n_pow2(m - 2); \
    FLINT_ASSERT(j_bits != 0); \
    CAT3(sd_ifft_basecase, n, 0)(Q, X+0*l, 4*j_mr+3, 2+j_bits); \
    CAT3(sd_ifft_basecase, n, 0)(Q, X+1*l, 4*j_mr+2, 2+j_bits); \
    CAT3(sd_ifft_basecase, n, 0)(Q, X+2*l, 4*j_mr+1, 2+j_bits); \
    CAT3(sd_ifft_basecase, n, 0)(Q, X+3*l, 4*j_mr+0, 2+j_bits); \
    { \
        _RADIX_4_REVERSE_PARAM_J_IS_NZ(VECND, Q, j_mr, j_bits) \
        ulong i = 0; do { \
            _RADIX_4_REVERSE_MOTH_J_IS_NZ(VECND, X+0*l+i, X+1*l+i, X+2*l+i, X+3*l+i) \
        } while (i += N, i < l); \
        FLINT_ASSERT(i == l); \
    } \
}

EXTEND_BASECASE(4, 6)
EXTEND_BASECASE(6, 8)
#undef EXTEND_BASECASE

/* parameter 1: j can be zero */
void sd_ifft_base_1(const sd_fft_lctx_t Q, ulong I, ulong j)
{
    ulong j_bits, j_mr;
    double* x = sd_fft_lctx_blk_index(Q, I);

    SET_J_BITS_AND_J_MR(j_bits, j_mr, j);

    if (j == 0)
        sd_ifft_basecase_8_1(Q, x, j_mr, j_bits);
    else
        sd_ifft_basecase_8_0(Q, x, j_mr, j_bits);
}

/* parameter 0: j cannot be zero */
void sd_ifft_base_0(const sd_fft_lctx_t Q, ulong I, ulong j)
{
    ulong j_bits, j_mr;
    double* x = sd_fft_lctx_blk_index(Q, I);

    FLINT_ASSERT(j != 0);

    SET_J_BITS_AND_J_MR(j_bits, j_mr, j);

    sd_ifft_basecase_8_0(Q, x, j_mr, j_bits);
}

/***************** inverse butterfy with truncation **************************/

/* the legal function are opt-in, and illegal/unimplemented are NULL */

#define radix_4_moth_inv_trunc_block_0_1_0 NULL
#define radix_4_moth_inv_trunc_block_0_1_1 NULL
#define radix_4_moth_inv_trunc_block_0_2_0 NULL
#define radix_4_moth_inv_trunc_block_0_2_1 NULL
#define radix_4_moth_inv_trunc_block_0_3_0 NULL
#define radix_4_moth_inv_trunc_block_0_3_1 NULL
#define radix_4_moth_inv_trunc_block_0_4_0 NULL
#define radix_4_moth_inv_trunc_block_1_2_0 NULL
#define radix_4_moth_inv_trunc_block_1_2_1 NULL
#define radix_4_moth_inv_trunc_block_1_3_0 NULL
#define radix_4_moth_inv_trunc_block_1_3_1 NULL
#define radix_4_moth_inv_trunc_block_2_1_0 NULL
#define radix_4_moth_inv_trunc_block_2_1_1 NULL
#define radix_4_moth_inv_trunc_block_2_3_0 NULL
#define radix_4_moth_inv_trunc_block_2_3_1 NULL
#define radix_4_moth_inv_trunc_block_3_1_0 NULL
#define radix_4_moth_inv_trunc_block_3_1_1 NULL
#define radix_4_moth_inv_trunc_block_3_2_0 NULL
#define radix_4_moth_inv_trunc_block_3_2_1 NULL
#define radix_4_moth_inv_trunc_block_4_1_0 NULL
#define radix_4_moth_inv_trunc_block_4_1_1 NULL
#define radix_4_moth_inv_trunc_block_4_2_0 NULL
#define radix_4_moth_inv_trunc_block_4_2_1 NULL
#define radix_4_moth_inv_trunc_block_4_3_0 NULL
#define radix_4_moth_inv_trunc_block_4_3_1 NULL
#define radix_4_moth_inv_trunc_block_4_4_0 NULL
#define radix_4_moth_inv_trunc_block_4_4_1 NULL

/*
k = 2, n = 3, z = 4, f = true
[      -r + 1           r + 1         2   r*w^3]
[        2//w           -2//w         0    -w^2]
[(r + 1)//w^2   (-r + 1)//w^2   -2//w^2    -r*w]
[          -r               r         1   r*w^3]
*/
static void radix_4_moth_inv_trunc_block_3_4_1(
    const sd_fft_lctx_t Q,
    ulong j, ulong j_bits,
    double* X0, double* X1, double* X2, double* X3)
{
    ulong j_mr = n_pow2(j_bits) - 1 - j;
    ulong j_r = j & (n_pow2(j_bits)/2 - 1);
    VECMD n    = VECMOP(set_d)(Q->p);
    VECMD ninv = VECMOP(set_d)(Q->pinv);
    double W  = FLINT_UNLIKELY(j == 0) ? -1.0 : Q->w2tab[1+j_bits][2*j_mr+1];
    double W2 = FLINT_UNLIKELY(j == 0) ? -1.0 : Q->w2tab[0+j_bits][j_mr];
    double twoW = vec1d_reduce_pm1n_to_pmhn(-2*W, Q->p);
    VECMD f0 = VECMOP(set_d)(Q->w2tab[0][1]);         /* r */
    VECMD f1 = VECMOP(set_d)(twoW);                   /* 2*w^-1 */
    VECMD f2 = VECMOP(set_d)(2);
    VECMD f3 = VECMOP(set_d)(W2);                     /* -w^-2 */
    double rw = FLINT_UNLIKELY(j == 0) ? Q->w2tab[0][1] : Q->w2tab[1+j_bits][2*j_r+1];
    VECMD fr = VECMOP(set_d)(rw);                         /* r*w */
    VECMD fq = VECMOP(set_d)(Q->w2tab[0+j_bits][j_r]);    /* w^2 */
    VECMD fp_ = VECMOP(mulmod)(fr, fq, n, ninv);
    VECMD fp = VECMOP(reduce_pm1n_to_pmhn)(fp_, n);       /* r*w^3 */
    ulong i = 0; do {
        VECMD a, b, c, d, u, v, p, q, r;
        a = VECMOP(load)(X0+i);
        b = VECMOP(load)(X1+i);
        c = VECMOP(load)(X2+i);
        d = VECMOP(load)(X3+i);
        u = VECMOP(add)(a, b);
        v = VECMOP(sub)(a, b);
        p = VECMOP(mulmod)(d, fp, n, ninv);
        q = VECMOP(mulmod)(d, fq, n, ninv);
        r = VECMOP(mulmod)(d, fr, n, ninv);
        c = VECMOP(reduce_to_pm1n)(c, n, ninv);
        u = VECMOP(reduce_to_pm1n)(u, n, ninv);
        b = VECMOP(mulmod)(v, f1, n, ninv);
        v = VECMOP(mulmod)(v, f0, n, ninv);
        d = VECMOP(sub)(c, v);
        c = VECMOP(fmsub)(f2, c, v);
        a = VECMOP(add)(c, u);
        c = VECMOP(sub)(c, u);
        c = VECMOP(mulmod)(c, f3, n, ninv);
        VECMOP(store)(X0+i, VECMOP(add)(a, p));
        VECMOP(store)(X1+i, VECMOP(sub)(b, q));
        VECMOP(store)(X2+i, VECMOP(sub)(c, r));
        VECMOP(store)(X3+i, VECMOP(add)(d, p));
    } while (i += M, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/*
k = 2, n = 3, z = 4, f = false
[      -r + 1           r + 1         2   r*w^3]
[        2//w           -2//w         0    -w^2]
[(r + 1)//w^2   (-r + 1)//w^2   -2//w^2    -r*w]
*/
static void radix_4_moth_inv_trunc_block_3_4_0(
    const sd_fft_lctx_t Q,
    ulong j, ulong j_bits,
    double* X0, double* X1, double* X2, double* X3)
{
    ulong j_mr = n_pow2(j_bits) - 1 - j;
    ulong j_r = j & (n_pow2(j_bits)/2 - 1);
    VECMD n    = VECMOP(set_d)(Q->p);
    VECMD ninv = VECMOP(set_d)(Q->pinv);
    double W  = FLINT_UNLIKELY(j == 0) ? -1.0 : Q->w2tab[1+j_bits][2*j_mr+1];
    double W2 = FLINT_UNLIKELY(j == 0) ? -1.0 : Q->w2tab[0+j_bits][j_mr];
    double twoW = vec1d_reduce_pm1n_to_pmhn(-2*W, Q->p);
    VECMD f0 = VECMOP(set_d)(Q->w2tab[0][1]);         /* r */
    VECMD f1 = VECMOP(set_d)(twoW);                   /* 2*w^-1 */
    VECMD f2 = VECMOP(set_d)(2);
    VECMD f3 = VECMOP(set_d)(W2);                     /* -w^-2 */
    double rw = FLINT_UNLIKELY(j == 0) ? Q->w2tab[0][1] : Q->w2tab[1+j_bits][2*j_r+1];
    VECMD fr = VECMOP(set_d)(rw);                         /* r*w */
    VECMD fq = VECMOP(set_d)(Q->w2tab[0+j_bits][j_r]);    /* w^2 */
    VECMD fp_ = VECMOP(mulmod)(fr, fq, n, ninv);
    VECMD fp = VECMOP(reduce_pm1n_to_pmhn)(fp_, n);       /* r*w^3 */
    ulong i = 0; do {
        VECMD a, b, c, d, u, v, p, q, r;
        a = VECMOP(load)(X0+i);
        b = VECMOP(load)(X1+i);
        c = VECMOP(load)(X2+i);
        d = VECMOP(load)(X3+i);
        u = VECMOP(add)(a, b);
        v = VECMOP(sub)(a, b);
        p = VECMOP(mulmod)(d, fp, n, ninv);
        q = VECMOP(mulmod)(d, fq, n, ninv);
        r = VECMOP(mulmod)(d, fr, n, ninv);
        c = VECMOP(reduce_to_pm1n)(c, n, ninv);
        u = VECMOP(reduce_to_pm1n)(u, n, ninv);
        b = VECMOP(mulmod)(v, f1, n, ninv);
        v = VECMOP(mulmod)(v, f0, n, ninv);
        c = VECMOP(fmsub)(f2, c, v);
        a = VECMOP(add)(c, u);
        c = VECMOP(sub)(c, u);
        c = VECMOP(mulmod)(c, f3, n, ninv);
        a = VECMOP(add)(a, p);
        b = VECMOP(sub)(b, q);
        c = VECMOP(sub)(c, r);
        VECMOP(store)(X0+i, a);
        VECMOP(store)(X1+i, b);
        VECMOP(store)(X2+i, c);
    } while (i += M, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/*
k = 2, n = 3, z = 3, f = true
[      -r + 1           r + 1         2]
[        2//w           -2//w         0]
[(r + 1)//w^2   (-r + 1)//w^2   -2//w^2]
[          -r               r         1]

    {x0, x1, x3, x4} = {        -r*(x0 - x1) + (x0 + x1) + 2*x2,
                        2*w^-1*    (x0 - x1),
                         -w^-2*(-r*(x0 - x1) - (x0 + x1) + 2*x2),
                                -r*(x0 - x1)             +   x2  }
*/
static void radix_4_moth_inv_trunc_block_3_3_1(
    const sd_fft_lctx_t Q,
    ulong j, ulong j_bits,
    double* X0, double* X1, double* X2, double* X3)
{
    ulong j_mr = n_pow2(j_bits) - 1 - j;
    VECND n    = VECNOP(set_d)(Q->p);
    VECND ninv = VECNOP(set_d)(Q->pinv);
    double W  = (j == 0) ? -1.0 : Q->w2tab[1+j_bits][2*j_mr+1];
    double W2 = (j == 0) ? -1.0 : Q->w2tab[0+j_bits][j_mr];
    double twoW = vec1d_reduce_pm1n_to_pmhn(-2*W, Q->p);
    VECND f0 = VECNOP(set_d)(Q->w2tab[0][1]);  /* r */
    VECND f1 = VECNOP(set_d)(twoW);            /* 2*w^-1 */
    VECND f2 = VECNOP(set_d)(2);
    VECND f3 = VECNOP(set_d)(W2);              /* -w^-2 */
    ulong i = 0; do {
        VECND a, b, c, u, v;
        a = VECNOP(load)(X0+i);
        b = VECNOP(load)(X1+i);
        c = VECNOP(load)(X2+i);
        v = VECNOP(sub)(a, b);
        VECNOP(store)(X1+i, VECNOP(mulmod)(v, f1, n, ninv));
        c = VECNOP(reduce_to_pm1n)(c, n, ninv);
        v = VECNOP(mulmod)(v, f0, n, ninv);
        VECNOP(store)(X3+i, VECNOP(sub)(c, v));
        u = VECNOP(reduce_to_pm1n)(VECNOP(add)(a, b), n, ninv);
        c = VECNOP(fnmadd)(f2, c, v);
        a = VECNOP(sub)(u, c);
        c = VECNOP(add)(u, c);
        c = VECNOP(nmulmod)(c, f3, n, ninv);
        VECNOP(store)(X0+i, a);
        VECNOP(store)(X2+i, c);
    } while (i += N, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/*
k = 2, n = 3, z = 3, f = false
[      -r + 1           r + 1         2]
[        2//w           -2//w         0]
[(r + 1)//w^2   (-r + 1)//w^2   -2//w^2]

    {x0, x1, x3} = {        -r*(x0 - x1) + (x0 + x1) + 2*x2,
                    2*w^-1*(x0 - x1),
                     -w^-2*(-r*(x0 - x1) - (x0 + x1) + 2*x2)}

                 = {        2*x2 - r*v + u,
                    2*w^-1*v,
                     -w^-2*(2*x2 - r*v - u)}
*/
static void radix_4_moth_inv_trunc_block_3_3_0(
    const sd_fft_lctx_t Q,
    ulong j, ulong j_bits,
    double* X0, double* X1, double* X2, double* X3)
{
    ulong j_mr = n_pow2(j_bits) - 1 - j;
    VECND n    = VECNOP(set_d)(Q->p);
    VECND ninv = VECNOP(set_d)(Q->pinv);
    double W  = (j == 0) ? -1.0 : Q->w2tab[1+j_bits][2*j_mr+1];
    double W2 = (j == 0) ? -1.0 : Q->w2tab[0+j_bits][j_mr];
    double twoW = vec1d_reduce_pm1n_to_pmhn(-2*W, Q->p);
    VECND f0 = VECNOP(set_d)(Q->w2tab[0][1]);  /* r */
    VECND f1 = VECNOP(set_d)(twoW);            /* 2*w^-1 */
    VECND f2 = VECNOP(set_d)(2);
    VECND f3 = VECNOP(set_d)(W2);              /* -w^-2 */
    ulong i = 0; do {
        VECND a, b, c, u, v;
        a = VECNOP(load)(X0+i);
        b = VECNOP(load)(X1+i);
        c = VECNOP(load)(X2+i);
        v = VECNOP(sub)(a, b);
        VECNOP(store)(X1+i, VECNOP(mulmod)(v, f1, n, ninv));
        c = VECNOP(reduce_to_pm1n)(c, n, ninv);
        v = VECNOP(mulmod)(v, f0, n, ninv);
        u = VECNOP(add)(a, b);
        u = VECNOP(reduce_to_pm1n)(u, n, ninv);
        c = VECNOP(fnmadd)(f2, c, v);
        a = VECNOP(sub)(u, c);
        c = VECNOP(add)(u, c);
        c = VECNOP(nmulmod)(c, f3, n, ninv);
        VECNOP(store)(X0+i, a);
        VECNOP(store)(X2+i, c);
    } while (i += N, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/*
k = 2, n = 2, z = 4, f = true
[            2                2        -w^2             0]
[         2//w            -2//w           0          -w^2]
[1//2*r + 1//2   -1//2*r + 1//2   -1//2*w^2   -1//2*r*w^3]
*/
static void radix_4_moth_inv_trunc_block_2_4_1(
    const sd_fft_lctx_t Q,
    ulong j, ulong j_bits,
    double* X0, double* X1, double* X2, double* X3)
{
    ulong j_mr = n_pow2(j_bits) - 1 - j;
    ulong j_r = j & (n_pow2(j_bits)/2 - 1);
    VECMD n    = VECMOP(set_d)(Q->p);
    VECMD ninv = VECMOP(set_d)(Q->pinv);
    double W = FLINT_UNLIKELY(j == 0) ? -1.0 : Q->w2tab[1+j_bits][2*j_mr+1];
    double rw = FLINT_UNLIKELY(j == 0) ? Q->w2tab[0][1] : Q->w2tab[1+j_bits][2*j_r+1];
    double w = Q->w2tab[j_bits][j_r];
    double twoW = vec1d_reduce_pm1n_to_pmhn(-2*W, Q->p);
    double rw3 = vec1d_mulmod(w, rw, Q->p, Q->pinv);
    VECMD f0 = VECMOP(set_d)(2);
    VECMD f1 = VECMOP(set_d)(twoW);                                   /* 2*w^-1 */
    VECMD f2 = VECMOP(set_d)(vec1d_fnmadd(0.5, Q->p, 0.5));           /* 1/2 */
    VECMD f3 = VECMOP(set_d)(Q->w2tab[0][1]);                         /* r */
    VECMD f4 = VECMOP(set_d)(w);                                      /* w^2 */
    VECMD f5 = VECMOP(set_d)(vec1d_reduce_pm1n_to_pmhn(rw3, Q->p));   /* r*w^3 */
    ulong i = 0; do {
        VECMD a, b, u, v, s, t, g, h, p, q, r;
        u = VECMOP(load)(X0+i);
        v = VECMOP(load)(X1+i);
        a = VECMOP(load)(X2+i);
        b = VECMOP(load)(X3+i);
        p = VECMOP(mulmod)(a, f4, n, ninv);
        q = VECMOP(mulmod)(b, f4, n, ninv);
        r = VECMOP(mulmod)(b, f5, n, ninv);
        s = VECMOP(add)(u, v);
        s = VECMOP(reduce_to_pm1n)(s, n, ninv);
        t = VECMOP(sub)(u, v);
        g = VECMOP(mulmod)(s, f0, n, ninv);
        h = VECMOP(mulmod)(t, f1, n, ninv);
        t = VECMOP(mulmod)(t, f3, n, ninv);
        VECMOP(store)(X0+i, VECMOP(sub)(g, p));
        VECMOP(store)(X1+i, VECMOP(sub)(h, q));
        u = VECMOP(add)(s, t);
        v = VECMOP(add)(p, r);
        u = VECMOP(sub)(u, v);
        u = VECMOP(mulmod)(u, f2, n, ninv);
        VECMOP(store)(X2+i, u);
    } while (i += M, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/*
k = 2, n = 2, z = 4, f = false
[   2       2   -w^2      0]
[2//w   -2//w      0   -w^2]
*/
static void radix_4_moth_inv_trunc_block_2_4_0(
    const sd_fft_lctx_t Q,
    ulong j, ulong j_bits,
    double* X0, double* X1, double* X2, double* X3)
{
    ulong j_mr = n_pow2(j_bits) - 1 - j;
    ulong j_r = j & (n_pow2(j_bits)/2 - 1);
    VECND n    = VECNOP(set_d)(Q->p);
    VECND ninv = VECNOP(set_d)(Q->pinv);
    double wi = FLINT_UNLIKELY(j == 0) ? -1.0 : Q->w2tab[1+j_bits][2*j_mr+1];
    VECND w2 = VECNOP(set_d)(Q->w2tab[j_bits][j_r]);
    VECND twowi = VECNOP(set_d)(vec1d_reduce_pm1n_to_pmhn(-2*wi, Q->p));
    ulong i = 0; do {
        VECND a, b, c, d, u, v;
        a = VECNOP(load)(X0+i);
        b = VECNOP(load)(X1+i);
        c = VECNOP(load)(X2+i);
        d = VECNOP(load)(X3+i);
        c = VECNOP(mulmod)(c, w2, n, ninv);
        d = VECNOP(mulmod)(d, w2, n, ninv);
        u = VECNOP(add)(a, b);
        v = VECNOP(sub)(a, b);
        u = VECNOP(add)(u, u);
        u = VECNOP(reduce_to_pm1n)(u, n, ninv);
        v = VECNOP(mulmod)(v, twowi, n, ninv);
        u = VECNOP(sub)(u, c);
        v = VECNOP(sub)(v, d);
        VECNOP(store)(X0+i, u);
        VECNOP(store)(X1+i, v);
    } while (i += N, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/*
k = 2, n = 2, z = 2, f = true
[            2                2]
[         2//w            -2//w]
[1//2*r + 1//2   -1//2*r + 1//2]

{x0, x1, x2} = {2*(x0 + x1), 2*w^-1*(x0 - x1), (x0+x1)/2 + (x0-x1)*i/2}
*/
static void radix_4_moth_inv_trunc_block_2_2_1(
    const sd_fft_lctx_t Q,
    ulong j, ulong j_bits,
    double* X0, double* X1, double* X2, double* X3)
{
    ulong j_mr = n_pow2(j_bits) - 1 - j;
    VECND n    = VECNOP(set_d)(Q->p);
    VECND ninv = VECNOP(set_d)(Q->pinv);
    double W = (j == 0) ? -1.0 : Q->w2tab[1+j_bits][2*j_mr+1];
    VECND c1 = VECNOP(set_d)(vec1d_reduce_pm1n_to_pmhn(-2*W, Q->p));  /* 2/w */
    VECND c2 = VECNOP(set_d)(vec1d_fnmadd(0.5, Q->p, 0.5));           /* 1/2 */
    VECND c3 = VECNOP(set_d)(Q->w2tab[1][0]);                         /* r */
    ulong i = 0; do {
        VECND u, v, s, t;
        u = VECNOP(load)(X0+i);
        v = VECNOP(load)(X1+i);
        s = VECNOP(add)(u, v);
        t = VECNOP(sub)(u, v);
        u = VECNOP(add)(s, s);
        u = VECNOP(reduce_to_pm1n)(u, n, ninv);
        v = VECNOP(mulmod)(t, c1, n, ninv);
        t = VECNOP(mulmod)(t, c3, n, ninv);
        s = VECNOP(add)(s, t);
        s = VECNOP(mulmod)(s, c2, n, ninv);
        VECNOP(store)(X0+i, u);
        VECNOP(store)(X1+i, v);
        VECNOP(store)(X2+i, s);
    } while (i += N, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/*
k = 2, n = 2, z = 2, f = 0

{x0, x1} = {2*(x0 + x1), 2*w^-1*(x0 - x1)}
*/
static void radix_4_moth_inv_trunc_block_2_2_0(
    const sd_fft_lctx_t Q,
    ulong j, ulong j_bits,
    double* X0, double* X1, double* X2, double* X3)
{
    ulong j_mr = n_pow2(j_bits) - 1 - j;
    VECND n    = VECNOP(set_d)(Q->p);
    VECND ninv = VECNOP(set_d)(Q->pinv);
    double W = (j == 0) ? -1.0 : Q->w2tab[1+j_bits][2*j_mr+1];
    VECND c0 = VECNOP(set_d)(2);
    VECND c1 = VECNOP(set_d)(vec1d_reduce_pm1n_to_pmhn(-2*W, Q->p));
    ulong i = 0; do {
        VECND u, v, s, t;
        u = VECNOP(load)(X0+i);
        v = VECNOP(load)(X1+i);
        s = VECNOP(add)(u, v);
        t = VECNOP(sub)(u, v);
        u = VECNOP(mulmod)(s, c0, n, ninv);
        v = VECNOP(mulmod)(t, c1, n, ninv);
        VECNOP(store)(X0+i, u);
        VECNOP(store)(X1+i, v);
    } while (i += N, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}


/*
k = 2, n = 1, z = 4, f = true
[4        -w   -w^2        -w^3]
[1   -1//2*w      0   -1//2*w^3]
*/
static void radix_4_moth_inv_trunc_block_1_4_1(
    const sd_fft_lctx_t Q,
    ulong j, ulong j_bits,
    double* X0, double* X1, double* X2, double* X3)
{
    ulong j_r = j & (n_pow2(j_bits)/2 - 1);
    double W2 = Q->w2tab[j_bits][j_r];
    double W  = FLINT_UNLIKELY(j == 0) ? 1.0 : Q->w2tab[1+j_bits][2*j_r];
    VECMD n    = VECMOP(set_d)(Q->p);
    VECMD ninv = VECMOP(set_d)(Q->pinv);
    VECMD f2 = VECMOP(set_d)(2);
    VECMD w2 = VECMOP(set_d)(W2);
    double ha = vec1d_fnmadd(0.5, Q->p, 0.5);
    double haW = vec1d_mulmod(W, ha, Q->p, Q->pinv);
    VECMD wo2 = VECMOP(set_d)(vec1d_reduce_pm1n_to_pmhn(haW, Q->p));
    ulong i = 0; do {
        VECMD a, b, c, d, u;
        a = VECMOP(load)(X0+i);
        a = VECMOP(reduce_to_pm1n)(a, n, ninv);
        b = VECMOP(load)(X1+i);
        c = VECMOP(load)(X2+i);
        d = VECMOP(load)(X3+i);
        c = VECMOP(nmulmod)(c, w2, n, ninv);
        d = VECMOP(mulmod)(d, w2, n, ninv);
        b = VECMOP(add)(b, d);
        b = VECMOP(mulmod)(b, wo2, n, ninv);
        u = VECMOP(fnmadd)(f2, a, b);
        b = VECMOP(sub)(a, b);
        a = VECMOP(reduce_to_pm1n)(VECMOP(add)(u, u), n, ninv);
        a = VECMOP(sub)(c, a);
        VECMOP(store)(X0+i, a);
        VECMOP(store)(X1+i, b);
    } while (i += M, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/*
k = 2, n = 1, z = 4, f = false
[4   -w   -w^2   -w^3]
*/
static void radix_4_moth_inv_trunc_block_1_4_0(
    const sd_fft_lctx_t Q,
    ulong j, ulong j_bits,
    double* X0, double* X1, double* X2, double* X3)
{
    ulong j_r = j & (n_pow2(j_bits)/2 - 1);
    VECND n    = VECNOP(set_d)(Q->p);
    VECND ninv = VECNOP(set_d)(Q->pinv);
    VECND f1 = VECNOP(set_d)(4.0);
    VECND w2 = VECNOP(set_d)(Q->w2tab[0+j_bits][j_r]);
    VECND w  = VECNOP(set_d)(FLINT_UNLIKELY(j == 0) ? 1.0 : Q->w2tab[1+j_bits][2*j_r]);
    ulong i = 0; do {
        VECND a, b, c, d;
        a = VECNOP(load)(X0+i);
        b = VECNOP(load)(X1+i);
        c = VECNOP(load)(X2+i);
        d = VECNOP(load)(X3+i);
        a = VECNOP(mul)(a, f1);
        d = VECNOP(mulmod)(d, w, n, ninv);
        a = VECNOP(reduce_to_pm1n)(a, n, ninv);
        b = VECNOP(mulmod)(b, w, n, ninv);
        a = VECNOP(sub)(a, b);
        c = VECNOP(add)(c, d);
        c = VECNOP(mulmod)(c, w2, n, ninv);
        a = VECNOP(sub)(a, c);
        VECNOP(store)(X0+i, a);
    } while (i += N, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/*
k = 2, n = 1, z = 1, f = true
[4]
[1]
*/
static void radix_4_moth_inv_trunc_block_1_1_1(
    const sd_fft_lctx_t Q,
    ulong j, ulong j_bits,
    double* X0, double* X1, double* X2, double* X3)
{
    VECND n    = VECNOP(set_d)(Q->p);
    VECND ninv = VECNOP(set_d)(Q->pinv);
    VECND f = VECNOP(set_d)(4.0);
    ulong i = 0; do {
        VECND a, b;
        a = VECNOP(load)(X0+i);
        b = VECNOP(mul)(f, a);
        VECNOP(store)(X0+i, VECNOP(reduce_to_pm1n)(b, n, ninv));
        VECNOP(store)(X1+i, VECNOP(reduce_to_pm1n)(a, n, ninv));
    } while (i += N, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/*
k = 2, n = 1, z = 1, f = false
[4]
*/
static void radix_4_moth_inv_trunc_block_1_1_0(
    const sd_fft_lctx_t Q,
    ulong j, ulong j_bits,
    double* X0, double* X1, double* X2, double* X3)
{
    VECND n    = VECNOP(set_d)(Q->p);
    VECND ninv = VECNOP(set_d)(Q->pinv);
    VECND f = VECNOP(set_d)(4.0);
    ulong i = 0; do {
        VECND a, b;
        a = VECNOP(load)(X0+i);
        b = VECNOP(mul)(f, a);
        VECNOP(store)(X0+i, VECNOP(reduce_to_pm1n)(b, n, ninv));
    } while (i += N, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/*
k = 2, n = 0, z = 4, f = true
[1//4   1//4*w   1//4*w^2   1//4*w^3]
*/
static void radix_4_moth_inv_trunc_block_0_4_1(
    const sd_fft_lctx_t Q,
    ulong j, ulong j_bits,
    double* X0, double* X1, double* X2, double* X3)
{
    ulong j_r = j & (n_pow2(j_bits)/2 - 1);
    VECND n    = VECNOP(set_d)(Q->p);
    VECND ninv = VECNOP(set_d)(Q->pinv);
    VECND one4th = VECNOP(set_d)(vec1d_fnmadd(0.25, Q->p, 0.25));
    VECND w2 = VECNOP(set_d)(Q->w2tab[j_bits][j_r]);
    VECND w  = VECNOP(set_d)(FLINT_UNLIKELY(j == 0) ? 1.0 : Q->w2tab[1+j_bits][2*j_r]);
    ulong i = 0; do {
        VECND a, b, c, d;
        a = VECNOP(load)(X0+i);
        b = VECNOP(load)(X1+i);
        c = VECNOP(load)(X2+i);
        d = VECNOP(load)(X3+i);
        b = VECNOP(mulmod)(b, w, n, ninv);
        d = VECNOP(mulmod)(d, w, n, ninv);
        a = VECNOP(add)(a, b);
        c = VECNOP(add)(c, d);
        c = VECNOP(mulmod)(c, w2, n, ninv);
        a = VECNOP(add)(a, c);
        VECNOP(store)(X0+i, VECNOP(mulmod)(a, one4th, n, ninv));
    } while (i += N, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/************************ the recursive stuff ********************************/

void sd_ifft_main_block(
    const sd_fft_lctx_t Q,
    ulong I, /* starting index */
    ulong S, /* stride */
    ulong k, /* BLK_SZ transforms each of length 2^k */
    ulong j)
{
    ulong j_bits, j_mr;

    if (k > 2)
    {
        ulong k1 = k/2;
        ulong k2 = k - k1;

        /* row ffts */
        ulong l1 = n_pow2(k1);
        ulong b = 0; do {
            sd_ifft_main_block(Q, I + (b<<k2)*S, S, k2, (j<<k1) + b);
        } while (b++, b < l1);

        /* column ffts */
        ulong l2 = n_pow2(k2);
        ulong a = 0; do {
            sd_ifft_main_block(Q, I + a*S, S<<k2, k1, j);
        } while (a++, a < l2);

        return;
    }

    SET_J_BITS_AND_J_MR(j_bits, j_mr, j);

    if (k == 2)
    {
        double* X0 = sd_fft_lctx_blk_index(Q, I + S*0);
        double* X1 = sd_fft_lctx_blk_index(Q, I + S*1);
        double* X2 = sd_fft_lctx_blk_index(Q, I + S*2);
        double* X3 = sd_fft_lctx_blk_index(Q, I + S*3);
        if (FLINT_UNLIKELY(j_bits == 0))
        {
            _RADIX_4_REVERSE_PARAM_J_IS_Z(VECND, Q)
            ulong i = 0; do {
                _RADIX_4_REVERSE_MOTH_J_IS_Z(VECND, X0+i, X1+i, X2+i, X3+i);
            } while(i += N, i < BLK_SZ);
        }
        else
        {
            _RADIX_4_REVERSE_PARAM_J_IS_NZ(VECND, Q, j_mr, j_bits)
            ulong i = 0; do {
                _RADIX_4_REVERSE_MOTH_J_IS_NZ(VECND, X0+i, X1+i, X2+i, X3+i);
            } while(i += N, i < BLK_SZ);
        }
    }
    else if (k == 1)
    {
        double* X0 = sd_fft_lctx_blk_index(Q, I + S*0);
        double* X1 = sd_fft_lctx_blk_index(Q, I + S*1);
        if (FLINT_UNLIKELY(j_bits == 0))
        {
            _RADIX_2_REVERSE_PARAM_J_IS_Z(VECND, Q)
            ulong i = 0; do {
                _RADIX_2_REVERSE_MOTH_J_IS_Z(VECND, X0+i, X1+i);
            } while (i += N, i < BLK_SZ);
        }
        else
        {
            _RADIX_2_REVERSE_PARAM_J_IS_NZ(VECND, Q, j_mr, j_bits)
            ulong i = 0; do {
                _RADIX_2_REVERSE_MOTH_J_IS_NZ(VECND, X0+i, X1+i);
            } while (i += N, i < BLK_SZ);
        }
    }
}

void sd_ifft_main(
    const sd_fft_lctx_t Q,
    ulong I, /* starting index */
    ulong S, /* stride */
    ulong k, /* transform length BLK_SZ*2^k */
    ulong j)
{
    if (k > 2)
    {
        ulong k1 = k/2;
        ulong k2 = k - k1;

        ulong l1 = n_pow2(k1);
        ulong b = 0; do {
            sd_ifft_main(Q, I + b*(S<<k2), S, k2, (j<<k1) + b);
        } while (b++, b < l1);

        ulong l2 = n_pow2(k2);
        ulong a = 0; do {
            sd_ifft_main_block(Q, I + a*S, S<<k2, k1, j);
        } while (a++, a < l2);

        return;
    }

    if (k == 2)
    {
        /* k1 = 2; k2 = 0 */
        sd_ifft_base_1(Q, I + S*0, 4*j+0);
        sd_ifft_base_0(Q, I + S*1, 4*j+1);
        sd_ifft_base_0(Q, I + S*2, 4*j+2);
        sd_ifft_base_0(Q, I + S*3, 4*j+3);
        sd_ifft_main_block(Q, I, S, 2, j);
    }
    else if (k == 1)
    {
        /* k1 = 1; k2 = 0 */
        sd_ifft_base_1(Q, I + S*0, 2*j+0);
        sd_ifft_base_0(Q, I + S*1, 2*j+1);
        sd_ifft_main_block(Q, I, S, 1, j);
    }
    else
    {
        sd_ifft_base_1(Q, I, j);
    }
}

void sd_ifft_trunc_block(
    const sd_fft_lctx_t Q,
    ulong I, /* starting index */
    ulong S, /* stride */
    ulong k, /* BLK_SZ transforms each of length 2^k */
    ulong j,
    ulong z, /* actual trunc is z */
    ulong n, /* actual trunc is n */
    int f)
{
    FLINT_ASSERT(f == 0 || f == 1);
    FLINT_ASSERT(n <= z);
    FLINT_ASSERT(1 <= z && z <= n_pow2(k));
    FLINT_ASSERT(1 <= n+f && n+f <= n_pow2(k));

    if (!f && z == n && n == n_pow2(k))
    {
        sd_ifft_main_block(Q, I, S, k, j);
        return;
    }

    if (k == 2)
    {
#define IT(nn, zz, ff) CAT4(radix_4_moth_inv_trunc_block, nn, zz, ff)
#define LOOKUP_IT(nn, zz, ff) tab[(ulong)(ff) + 2*((zz)-1 + 4*(nn))]
        static void (*tab[5*4*2])(const sd_fft_lctx_t, ulong, ulong, double*, double*, double*, double*) =
            {IT(0,1,0),IT(0,1,1), IT(0,2,0),IT(0,2,1), IT(0,3,0),IT(0,3,1), IT(0,4,0),IT(0,4,1),
             IT(1,1,0),IT(1,1,1), IT(1,2,0),IT(1,2,1), IT(1,3,0),IT(1,3,1), IT(1,4,0),IT(1,4,1),
             IT(2,1,0),IT(2,1,1), IT(2,2,0),IT(2,2,1), IT(2,3,0),IT(2,3,1), IT(2,4,0),IT(2,4,1),
             IT(3,1,0),IT(3,1,1), IT(3,2,0),IT(3,2,1), IT(3,3,0),IT(3,3,1), IT(3,4,0),IT(3,4,1),
             IT(4,1,0),IT(4,1,1), IT(4,2,0),IT(4,2,1), IT(4,3,0),IT(4,3,1), IT(4,4,0),IT(4,4,1)};

        static void (*fxn)(const sd_fft_lctx_t, ulong, ulong, double*, double*, double*, double*);
        fxn = LOOKUP_IT(n,z,f);
        if (FLINT_LIKELY(fxn != NULL))
        {
            fxn(Q, j, n_nbits(j), sd_fft_lctx_blk_index(Q, I + S*0),
                                  sd_fft_lctx_blk_index(Q, I + S*1),
                                  sd_fft_lctx_blk_index(Q, I + S*2),
                                  sd_fft_lctx_blk_index(Q, I + S*3));
            return;
        }
#undef LOOKUP_IT
#undef IT
    }

    if (k > 1)
    {
        ulong k1 = k/2;
        ulong k2 = k - k1;

        ulong l2 = n_pow2(k2);
        ulong n1 = n >> k2;
        ulong n2 = n & (l2 - 1);
        ulong z1 = z >> k2;
        ulong z2 = z & (l2 - 1);
        int fp = n2 + f > 0;
        ulong z2p = n_min(l2, z);
        ulong m = n_min(n2, z2);
        ulong mp = n_max(n2, z2);

        /* complete rows */
        for (ulong b = 0; b < n1; b++)
            sd_ifft_main_block(Q, I + b*(S << k2), S, k2, (j << k1) + b);

        /* rightmost columns */
        for (ulong a = n2; a < z2p; a++)
            sd_ifft_trunc_block(Q, I + a*S, S << k2, k1, j, z1 + (a < mp), n1, fp);

        /* last partial row */
        if (fp)
            sd_ifft_trunc_block(Q, I + n1*(S << k2), S, k2, (j << k1) + n1, z2p, n2, f);

        /* leftmost columns */
        for (ulong a = 0; a < n2; a++)
            sd_ifft_trunc_block(Q, I + a*S, S << k2, k1, j, z1 + (a < m), n1 + 1, 0);

        return;
    }

    if (k == 1)
    {
#define IT(nn, zz, ff) CAT4(radix_2_moth_inv_trunc_block, nn, zz, ff)
#define LOOKUP_IT(nn, zz, ff) tab[(ulong)(ff) + 2*((zz)-1 + 2*(nn))]
        static void (*tab[3*2*2*2])(const sd_fft_lctx_t, ulong, double*, double*) =
            {IT(0,1,0),IT(0,1,1), IT(0,2,0),IT(0,2,1),
             IT(1,1,0),IT(1,1,1), IT(1,2,0),IT(1,2,1),
             IT(2,1,0),IT(2,1,1), IT(2,2,0),IT(2,2,1)};

        LOOKUP_IT(n, z, f)(Q, j, sd_fft_lctx_blk_index(Q, I + S*0),
                                 sd_fft_lctx_blk_index(Q, I + S*1));
        return;
#undef LOOKUP_IT
#undef IT
    }
}


void sd_ifft_trunc(
    const sd_fft_lctx_t Q,
    ulong I, // starting index
    ulong S, // stride
    ulong k, // transform length 2^(k + LG_BLK_SZ)
    ulong j,
    ulong z,   // actual trunc is z*BLK_SZ
    ulong n,   // actual trunc is n*BLK_SZ
    int f)
{
    FLINT_ASSERT(f == 0 || f == 1);
    FLINT_ASSERT(n <= z);
    FLINT_ASSERT(1 <= BLK_SZ*z && BLK_SZ*z <= n_pow2(k+LG_BLK_SZ));
    FLINT_ASSERT(1 <= BLK_SZ*n+f && BLK_SZ*n+f <= n_pow2(k+LG_BLK_SZ));

    if (k > 2)
    {
        ulong k1 = k/2;
        ulong k2 = k - k1;

        ulong l2 = n_pow2(k2);
        ulong n1 = n >> k2;
        ulong n2 = n & (l2 - 1);
        ulong z1 = z >> k2;
        ulong z2 = z & (l2 - 1);
        int fp = n2 + f > 0;
        ulong z2p = n_min(l2, z);
        ulong m = n_min(n2, z2);
        ulong mp = n_max(n2, z2);

        /* complete rows */
        for (ulong b = 0; b < n1; b++)
            sd_ifft_main(Q, I + b*(S << k2), S, k2, (j << k1) + b);

        /* rightmost columns */
        for (ulong a = n2; a < z2p; a++)
            sd_ifft_trunc_block(Q, I + a*S, S << k2, k1, j, z1 + (a < mp), n1, fp);

        /* last partial row */
        if (fp)
            sd_ifft_trunc(Q, I + n1*(S << k2), S, k2, (j << k1) + n1, z2p, n2, f);

        /* leftmost columns */
        for (ulong a = 0; a < n2; a++)
            sd_ifft_trunc_block(Q, I + a*S, S << k2, k1, j, z1 + (a < m), n1 + 1, 0);

        return;
    }

    if (k == 2)
    {
                   sd_ifft_base_1(Q, I + S*0, 4*j+0);
        if (n > 1) sd_ifft_base_0(Q, I + S*1, 4*j+1);
        if (n > 2) sd_ifft_base_0(Q, I + S*2, 4*j+2);
        if (n > 3) sd_ifft_base_0(Q, I + S*3, 4*j+3);
        sd_ifft_trunc_block(Q, I, S, 2, j, z, n, f);
        if (f) sd_ifft_trunc(Q, I + S*n, S, 0, 4*j+n, 1, 0, f);
    }
    else if (k == 1)
    {
                   sd_ifft_base_1(Q, I + S*0, 2*j+0);
        if (n > 1) sd_ifft_base_0(Q, I + S*1, 2*j+1);
        sd_ifft_trunc_block(Q, I, S, 1, j, z, n, f);
        if (f) sd_ifft_trunc(Q, I + S*n, S, 0, 2*j+n, 1, 0, f);
    }
    else
    {
        FLINT_ASSERT(!f);
        sd_ifft_base_1(Q, I, j);
    }
}
