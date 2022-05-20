/*
    Copyright (C) 2022 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fft_small.h"
#include "machine_vectors.h"

/************************** inverse butterfly *********************************
    2*a0 =      (b0 + b1)
    2*a1 = w^-1*(b0 - b1)
    W := -w^-1
*/
#define RADIX_2_REVERSE_PARAM(V, Q, j) \
    V W = V##_set_d((UNLIKELY((j) == 0)) ? -Q->w2s[0] : Q->w2s[(j)^(n_next_pow2m1(j)>>1)]); \
    V n    = V##_set_d(Q->p); \
    V ninv = V##_set_d(Q->pinv);

#define RADIX_2_REVERSE_MOTH(V, X0, X1) \
{ \
    V x0, x1, y0, y1; \
    x0 = V##_load(X0); \
    x1 = V##_load(X1); \
    y0 = V##_add(x0, x1); \
    y1 = V##_sub(x1, x0); \
    y0 = V##_reduce_to_pm1n(y0, n, ninv); \
    y1 = V##_mulmod2(y1, W, n, ninv); \
    V##_store(X0, y0); \
    V##_store(X1, y1); \
}

/****************** inverse butterfly with truncation ************************/

/* the legal functions must be implemented */
#define DEFINE_IT(nn, zz, ff) \
static void CAT4(radix_2_moth_inv_trunc_block, nn, zz, ff)( \
    const sd_fft_ctx_t Q, \
    ulong j, \
    double* X0, double* X1) \
{ \
    int l = 2; \
    flint_printf("function l = %d, n = %d, z = %d, f = %d", l, nn, zz, ff); \
    if (1 <= zz && zz <= 2 && nn <= zz && 1 <= nn+ff && nn+ff <= l) \
    { \
        flint_printf(" is not implemented\n"); \
    } \
    else \
    { \
        flint_printf(" does not exist and should not be called\n"); \
    } \
    fflush(stdout); \
    flint_abort(); \
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
    const sd_fft_ctx_t Q,
    ulong j,
    double* X0, double* X1)
{
    vec4d n    = vec4d_set_d(Q->p);
    vec4d ninv = vec4d_set_d(Q->pinv);
    vec4d w = vec4d_set_d(Q->w2s[j]);
    vec4d c = vec4d_set_d(2);
    ulong i = 0; do {
        vec4d u0, u1;
        u0 = vec4d_load(X0 + i);
        u1 = vec4d_load(X1 + i);
        u0 = vec4d_reduce_to_pm1n(u0, n, ninv);
        u1 = vec4d_mulmod2(u1, w, n, ninv);
        vec4d_store(X0 + i, vec4d_fmsub(c, u0, u1));
        vec4d_store(X1 + i, vec4d_sub(u0, u1));
    } while (i += 4, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/* {x0} = {2*x0 - w*x1} */
static void radix_2_moth_inv_trunc_block_1_2_0(
    const sd_fft_ctx_t Q,
    ulong j,
    double* X0, double* X1)
{
    vec4d n    = vec4d_set_d(Q->p);
    vec4d ninv = vec4d_set_d(Q->pinv);
    vec4d w = vec4d_set_d(Q->w2s[j]);
    vec4d c = vec4d_set_d(2);
    ulong i = 0; do {
        vec4d u0, u1;
        u0 = vec4d_load(X0 + i);
        u1 = vec4d_load(X1 + i);
        u0 = vec4d_reduce_to_pm1n(u0, n, ninv);
        u1 = vec4d_mulmod2(u1, w, n, ninv);
        vec4d_store(X0 + i, vec4d_fmsub(c, u0, u1));
    } while (i += 4, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/* {x0, x1} = {2*x0, x0} */
static void radix_2_moth_inv_trunc_block_1_1_1(
    const sd_fft_ctx_t Q,
    ulong j,
    double* X0, double* X1)
{
    vec4d n    = vec4d_set_d(Q->p);
    vec4d ninv = vec4d_set_d(Q->pinv);
    ulong i = 0; do {
        vec4d u0;
        u0 = vec4d_load(X0 + i);
        u0 = vec4d_reduce_to_pm1n(u0, n, ninv);
        vec4d_store(X0 + i, vec4d_add(u0, u0));
        vec4d_store(X1 + i, u0);
    } while (i += 4, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/* {x0} = {2*x0} */
static void radix_2_moth_inv_trunc_block_1_1_0(
    const sd_fft_ctx_t Q,
    ulong j,
    double* X0, double* X1)
{
    vec4d n    = vec4d_set_d(Q->p);
    vec4d ninv = vec4d_set_d(Q->pinv);
    ulong i = 0; do {
        vec4d u0;
        u0 = vec4d_load(X0 + i);
        u0 = vec4d_reduce_to_pm1n(u0, n, ninv);
        vec4d_store(X0 + i, vec4d_add(u0, u0));
    } while (i += 4, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/* {x0} = {(x0 + w*x1)/2} */
static void radix_2_moth_inv_trunc_block_0_2_1(
    const sd_fft_ctx_t Q,
    ulong j,
    double* X0, double* X1)
{
    vec4d n    = vec4d_set_d(Q->p);
    vec4d ninv = vec4d_set_d(Q->pinv);
    vec4d w = vec4d_set_d(Q->w2s[j]);
    vec4d c = vec4d_set_d(vec1d_fnmadd(0.5, Q->p, 0.5));
    ulong i = 0; do {
        vec4d u0, u1;
        u0 = vec4d_load(X0 + i);
        u1 = vec4d_load(X1 + i);
        u1 = vec4d_mulmod2(u1, w, n, ninv);
        u0 = vec4d_mulmod2(vec4d_add(u0, u1), c, n, ninv);
        vec4d_store(X0 + i, u0);
    } while (i += 4, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/* {x0} = {(x0)/2} */
static void radix_2_moth_inv_trunc_block_0_1_1(
    const sd_fft_ctx_t Q,
    ulong j,
    double* X0, double* X1)
{
    vec4d n    = vec4d_set_d(Q->p);
    vec4d ninv = vec4d_set_d(Q->pinv);
    vec4d c = vec4d_set_d(vec1d_fnmadd(0.5, Q->p, 0.5));
    ulong i = 0; do {
        vec4d u0;
        u0 = vec4d_load(X0 + i);
        u0 = vec4d_mulmod2(u0, c, n, ninv);
        vec4d_store(X0 + i, u0);
    } while (i += 4, i < BLK_SZ);
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
#define RADIX_4_REVERSE_PARAM(V, Q, j, jm, j_can_be_0) \
    FLINT_ASSERT((j) == 0 || (jm) == ((j)^(n_next_pow2m1(j)>>1))); \
    V W  = V##_set_d(((j_can_be_0) && UNLIKELY((j) == 0)) ? -Q->w2s[0] : Q->w2s[2*(jm)+1]); \
    V W2 = V##_set_d(((j_can_be_0) && UNLIKELY((j) == 0)) ? -Q->w2s[0] : Q->w2s[(jm)]); \
    V IW = V##_set_d(((j_can_be_0) && UNLIKELY((j) == 0)) ?  Q->w2s[1] : Q->w2s[2*(jm)]); \
    V n    = V##_set_d(Q->p); \
    V ninv = V##_set_d(Q->pinv); \

#define RADIX_4_REVERSE_PARAM_SURE(V, Q, j, jm, j_is_0) \
    FLINT_ASSERT((j) == 0 || (jm) == ((j)^(n_next_pow2m1(j)>>1))); \
    V W  = V##_set_d((j_is_0) ? -Q->w2s[0] : Q->w2s[2*(jm)+1]); \
    V W2 = V##_set_d((j_is_0) ? -Q->w2s[0] : Q->w2s[(jm)]); \
    V IW = V##_set_d((j_is_0) ?  Q->w2s[1] : Q->w2s[2*(jm)]); \
    V n    = V##_set_d(Q->p); \
    V ninv = V##_set_d(Q->pinv); \

#define RADIX_4_REVERSE_MOTH(V, X0, X1, X2, X3) \
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
    y2 = V##_mulmod2(y2, W, n, ninv); \
    y3 = V##_mulmod2(y3, IW, n, ninv); \
    x0 = V##_add(y0, y1); \
    x1 = V##_sub(y3, y2); \
    V##_store(X1, x1); \
    x2 = V##_sub(y1, y0); \
    x3 = V##_add(y3, y2); \
    x0 = V##_reduce_to_pm1n(x0, n, ninv); \
    x2 = V##_mulmod2(x2, W2, n, ninv); \
    x3 = V##_mulmod2(x3, W2, n, ninv); \
    V##_store(X0, x0); \
    V##_store(X2, x2); \
    V##_store(X3, x3); \
}

/* second parameter is the bool j == 0 */

/*
    These functions are disabled because the ifft expects input in
    slightly-worse-than-bit-reversed order as in basecase_4.
*/
#if 0
/* length 1 */
FLINT_FORCE_INLINE void
sd_ifft_basecase_0_0(const sd_fft_ctx_t Q, double* X, ulong j, ulong jm)
{
    FLINT_ASSERT(0 == (j == 0));
}

FLINT_FORCE_INLINE void
sd_ifft_basecase_0_1(const sd_fft_ctx_t Q, double* X, ulong j, ulong jm)
{
    FLINT_ASSERT(1 == (j == 0));
}

/* length 2 */
FLINT_FORCE_INLINE void
sd_ifft_basecase_1_0(const sd_fft_ctx_t Q, double* X, ulong j, ulong jm)
{
    RADIX_2_REVERSE_PARAM(vec1d, Q, j)
    RADIX_2_REVERSE_MOTH(vec1d, X+0, X+1);
}

FLINT_FORCE_INLINE void
sd_ifft_basecase_1_1(const sd_fft_ctx_t Q, double* X, ulong j, ulong jm)
{
    RADIX_2_REVERSE_PARAM(vec1d, Q, j)
    RADIX_2_REVERSE_MOTH(vec1d, X+0, X+1);
}

/* length 4 */
FLINT_FORCE_INLINE void
sd_ifft_basecase_2_0(const sd_fft_ctx_t Q, double* X, ulong j, ulong jm)
{
    RADIX_4_REVERSE_PARAM_SURE(vec1d, Q, j, jm, 0)
    RADIX_4_REVERSE_MOTH(vec1d, X+0, X+1, X+2, X+3);
}

FLINT_FORCE_INLINE void
sd_ifft_basecase_2_1(const sd_fft_ctx_t Q, double* X, ulong j, ulong jm)
{
    RADIX_4_REVERSE_PARAM_SURE(vec1d, Q, j, jm, 1)
    RADIX_4_REVERSE_MOTH(vec1d, X+0, X+1, X+2, X+3);
}
#endif

/* length 16 */
#define DEFINE_IT(j_is_0) \
FLINT_FORCE_INLINE void sd_ifft_basecase_4_##j_is_0(\
    const sd_fft_ctx_t Q, double* X, ulong j, ulong jm) \
{ \
    vec4d n    = vec4d_set_d(Q->p); \
    vec4d ninv = vec4d_set_d(Q->pinv); \
    vec4d W, W2, IW, u, v; \
    vec4d x0, x1, x2, x3, y0, y1, y2, y3; \
    FLINT_ASSERT(j_is_0 == (j == 0)); \
    FLINT_ASSERT(j == 0 || jm == (j^(n_next_pow2m1(j)>>1))); \
 \
    x0 = vec4d_load(X + 0); \
    x1 = vec4d_load(X + 4); \
    x2 = vec4d_load(X + 8); \
    x3 = vec4d_load(X + 12); \
 \
    if (j_is_0) \
    { \
        W  = vec4d_set_d4(-Q->w2s[0], Q->w2s[3], Q->w2s[7], Q->w2s[5]); \
        IW = vec4d_set_d4( Q->w2s[1], Q->w2s[2], Q->w2s[6], Q->w2s[4]); \
        W2 = vec4d_set_d4(-Q->w2s[0], Q->w2s[1], Q->w2s[3], Q->w2s[2]); \
    } \
    else \
    { \
        W2 = vec4d_load_aligned(Q->w2s + 4*jm); /* a b c d */ \
        W2 = vec4d_permute_3_2_1_0(W2);         /* d c b a */ \
        u = vec4d_load_aligned(Q->w2s + 8*jm);  /* 0 1 2 3 */ \
        v = vec4d_load_aligned(Q->w2s + 8*jm+4);/* 4 5 6 7 */ \
        W  = vec4d_unpackhi(u, v);              /* 1 5 3 7 */ \
        IW = vec4d_unpacklo(u, v);              /* 0 4 2 6 */ \
        W  = vec4d_permute_3_1_2_0(W);          /* 7 5 3 1 */ \
        IW = vec4d_permute_3_1_2_0(IW);         /* 6 4 2 0 */ \
    } \
 \
    y0 = vec4d_add(x0, x1); \
    y1 = vec4d_add(x2, x3); \
    y2 = vec4d_sub(x0, x1); \
    y3 = vec4d_sub(x3, x2); \
    y2 = vec4d_mulmod2(y2, W, n, ninv); \
    y3 = vec4d_mulmod2(y3, IW, n, ninv); \
    x0 = vec4d_add(y0, y1); \
    x1 = vec4d_sub(y3, y2); \
    x2 = vec4d_sub(y1, y0); \
    x3 = vec4d_add(y3, y2); \
    x0 = vec4d_reduce_to_pm1n(x0, n, ninv); \
    x2 = vec4d_mulmod2(x2, W2, n, ninv); \
    x3 = vec4d_mulmod2(x3, W2, n, ninv); \
 \
    VEC4D_TRANSPOSE(x0, x1, x2, x3, x0, x1, x2, x3); \
 \
    W  = vec4d_set_d(j_is_0 ? -Q->w2s[0] : Q->w2s[2*jm+1]); \
    IW = vec4d_set_d(j_is_0 ?  Q->w2s[1] : Q->w2s[2*jm+0]); \
    W2 = vec4d_set_d(j_is_0 ? -Q->w2s[0] : Q->w2s[jm]); \
    y0 = vec4d_add(x0, x1); \
    y1 = vec4d_add(x2, x3); \
    y2 = vec4d_sub(x0, x1); \
    y3 = vec4d_sub(x3, x2); \
    y2 = vec4d_mulmod2(y2, W, n, ninv); \
    y3 = vec4d_mulmod2(y3, IW, n, ninv); \
    x0 = vec4d_add(y0, y1); \
    x1 = vec4d_sub(y3, y2); \
    x2 = vec4d_sub(y1, y0); \
    x3 = vec4d_add(y3, y2); \
    x0 = vec4d_reduce_to_pm1n(x0, n, ninv); \
    x2 = vec4d_mulmod2(x2, W2, n, ninv); \
    x3 = vec4d_mulmod2(x3, W2, n, ninv); \
    vec4d_store(X+0, x0); \
    vec4d_store(X+4, x1); \
    vec4d_store(X+8, x2); \
    vec4d_store(X+12, x3); \
}

DEFINE_IT(0)
DEFINE_IT(1)
#undef DEFINE_IT


/* use with N = M-2 and M >= 6 */
#define EXTEND_BASECASE(N, M, j_is_0) \
void CAT3(sd_ifft_basecase, M, j_is_0)(const sd_fft_ctx_t Q, double* X, ulong j, ulong jm) \
{ \
    ulong l = n_pow2(M - 2); \
    FLINT_ASSERT(j_is_0 == (j == 0)); \
    FLINT_ASSERT(j == 0 || jm == (j^(n_next_pow2m1(j)>>1))); \
    CAT3(sd_ifft_basecase, N, j_is_0)(Q, X+0*l, 4*j+0, j_is_0 ? 0 : 4*jm+3); \
    CAT3(sd_ifft_basecase, N, 0     )(Q, X+1*l, 4*j+1, j_is_0 ? 1 : 4*jm+2); \
    CAT3(sd_ifft_basecase, N, 0     )(Q, X+2*l, 4*j+2, j_is_0 ? 3 : 4*jm+1); \
    CAT3(sd_ifft_basecase, N, 0     )(Q, X+3*l, 4*j+3, j_is_0 ? 2 : 4*jm+0); \
    { \
        RADIX_4_REVERSE_PARAM_SURE(vec8d, Q, j, jm, j_is_0) \
        ulong i = 0; do { \
            RADIX_4_REVERSE_MOTH(vec8d, X+0*l+i, X+1*l+i, X+2*l+i, X+3*l+i) \
        } while (i += 8, i < l); \
        FLINT_ASSERT(i == l); \
    } \
}

EXTEND_BASECASE(4, 6, 0)
EXTEND_BASECASE(4, 6, 1)
EXTEND_BASECASE(6, 8, 0)
EXTEND_BASECASE(6, 8, 1)
#undef EXTEND_BASECASE

/* parameter 1: j can be zero */
void sd_ifft_base_1(const sd_fft_ctx_t Q, double* d, ulong I, ulong j)
{
    ulong jm = j^(n_next_pow2m1(j)>>1);
    double* x = sd_fft_ctx_blk_index(Q, d, I);
    if (j == 0)
        sd_ifft_basecase_8_1(Q, x, j, jm);
    else
        sd_ifft_basecase_8_0(Q, x, j, jm);
}

/* parameter 0: j cannot be zero */
void sd_ifft_base_0(const sd_fft_ctx_t Q, double* d, ulong I, ulong j)
{
    ulong jm = j^(n_next_pow2m1(j)>>1);
    double* x = sd_fft_ctx_blk_index(Q, d, I);
    FLINT_ASSERT(j != 0);
    sd_ifft_basecase_8_0(Q, x, j, jm);
}

/***************** inverse butterfy with truncation **************************/

/* the legal function are opt-in but print an annoying message if missing */
#define DEFINE_IT(nn, zz, ff) \
static int CAT4(radix_4_moth_inv_trunc_block, nn, zz, ff)( \
    const sd_fft_ctx_t Q, \
    ulong j, ulong jm, \
    double* X0, double* X1, double* X2, double* X3) \
{ \
    int l = 4; \
    flint_printf("function l = %d, n = %d, z = %d, f = %d", l, nn, zz, ff); \
    if (1 <= zz && zz <= 4 && nn <= zz && 1 <= nn+ff && nn+ff <= l) \
    { \
        flint_printf(" is not implemented\n"); \
        fflush(stdout); \
    } \
    else \
    { \
        flint_printf(" does not exist and should not be called\n"); \
        fflush(stdout); \
        flint_abort(); \
    } \
    return 0; \
}

DEFINE_IT(0,1,0)
DEFINE_IT(0,1,1)
DEFINE_IT(0,2,0)
DEFINE_IT(0,2,1)
DEFINE_IT(0,3,0)
DEFINE_IT(0,3,1)
DEFINE_IT(0,4,0)
DEFINE_IT(1,2,0)
DEFINE_IT(1,2,1)
DEFINE_IT(1,3,0)
DEFINE_IT(1,3,1)
DEFINE_IT(2,1,0)
DEFINE_IT(2,1,1)
DEFINE_IT(2,3,0)
DEFINE_IT(2,3,1)
DEFINE_IT(3,1,0)
DEFINE_IT(3,1,1)
DEFINE_IT(3,2,0)
DEFINE_IT(3,2,1)
DEFINE_IT(4,1,0)
DEFINE_IT(4,1,1)
DEFINE_IT(4,2,0)
DEFINE_IT(4,2,1)
DEFINE_IT(4,3,0)
DEFINE_IT(4,3,1)
DEFINE_IT(4,4,0)
DEFINE_IT(4,4,1)
#undef DEFINE_IT

/*
k = 2, n = 3, z = 4, f = true
[      -r + 1           r + 1         2   r*w^3]
[        2//w           -2//w         0    -w^2]
[(r + 1)//w^2   (-r + 1)//w^2   -2//w^2    -r*w]
[          -r               r         1   r*w^3]
*/
static int radix_4_moth_inv_trunc_block_3_4_1(
    const sd_fft_ctx_t Q,
    ulong j, ulong jm,
    double* X0, double* X1, double* X2, double* X3)
{
    vec4d n    = vec4d_set_d(Q->p);
    vec4d ninv = vec4d_set_d(Q->pinv);
    double W  = UNLIKELY(j == 0) ? -Q->w2s[0] : Q->w2s[2*jm+1];
    double W2 = UNLIKELY(j == 0) ? -Q->w2s[0] : Q->w2s[jm];
    double twoW = vec1d_reduce_pm1n_to_pmhn(-2*W, Q->p);
    vec4d f0 = vec4d_set_d(Q->w2s[1]);              /* r */
    vec4d f1 = vec4d_set_d(twoW);                   /* 2*w^-1 */
    vec4d f2 = vec4d_set_d(2);
    vec4d f3 = vec4d_set_d(W2);                     /* -w^-2 */
    vec4d fr = vec4d_set_d(-Q->w2s[2*j+1]);         /* -r*w */
    vec4d fq = vec4d_set_d(-Q->w2s[j]);             /* -w^2 */
    vec4d fp_ = vec4d_mulmod2(fr, fq, n, ninv);
    vec4d fp = vec4d_reduce_pm1n_to_pmhn(fp_, n);   /* r*w^3 */
    ulong i = 0; do {
        vec4d a, b, c, d, u, v, p, q, r;
        a = vec4d_load(X0+i);
        b = vec4d_load(X1+i);
        c = vec4d_load(X2+i);
        d = vec4d_load(X3+i);
        u = vec4d_add(a, b);
        v = vec4d_sub(a, b);
        p = vec4d_mulmod2(d, fp, n, ninv);
        q = vec4d_mulmod2(d, fq, n, ninv);
        r = vec4d_mulmod2(d, fr, n, ninv);
        c = vec4d_reduce_to_pm1n(c, n, ninv);
        u = vec4d_reduce_to_pm1n(u, n, ninv);
        b = vec4d_mulmod2(v, f1, n, ninv);
        v = vec4d_mulmod2(v, f0, n, ninv);
        d = vec4d_sub(c, v);
        c = vec4d_fmsub(f2, c, v);
        a = vec4d_add(c, u);
        c = vec4d_sub(c, u);
        c = vec4d_mulmod2(c, f3, n, ninv);
        vec4d_store(X0+i, vec4d_add(a, p));
        vec4d_store(X1+i, vec4d_add(b, q));
        vec4d_store(X2+i, vec4d_add(c, r));
        vec4d_store(X3+i, vec4d_add(d, p));
    } while (i += 4, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
    return 1;
}

/*
k = 2, n = 3, z = 4, f = false
[      -r + 1           r + 1         2   r*w^3]
[        2//w           -2//w         0    -w^2]
[(r + 1)//w^2   (-r + 1)//w^2   -2//w^2    -r*w]
*/
static int radix_4_moth_inv_trunc_block_3_4_0(
    const sd_fft_ctx_t Q,
    ulong j, ulong jm,
    double* X0, double* X1, double* X2, double* X3)
{
    vec4d n    = vec4d_set_d(Q->p);
    vec4d ninv = vec4d_set_d(Q->pinv);
    double W  = UNLIKELY(j == 0) ? -Q->w2s[0] : Q->w2s[2*jm+1];
    double W2 = UNLIKELY(j == 0) ? -Q->w2s[0] : Q->w2s[jm];
    double twoW = vec1d_reduce_pm1n_to_pmhn(-2*W, Q->p);
    vec4d f0 = vec4d_set_d(Q->w2s[1]);              /* r */
    vec4d f1 = vec4d_set_d(twoW);                   /* 2*w^-1 */
    vec4d f2 = vec4d_set_d(2);
    vec4d f3 = vec4d_set_d(W2);                     /* -w^-2 */
    vec4d fr = vec4d_set_d(Q->w2s[2*j+1]);          /* r*w */
    vec4d fq = vec4d_set_d(Q->w2s[j]);              /* w^2 */
    vec4d fp_ = vec4d_mulmod2(fr, fq, n, ninv);
    vec4d fp = vec4d_reduce_pm1n_to_pmhn(fp_, n);   /* r*w^3 */
    ulong i = 0; do {
        vec4d a, b, c, d, u, v, p, q, r;
        a = vec4d_load(X0+i);
        b = vec4d_load(X1+i);
        c = vec4d_load(X2+i);
        d = vec4d_load(X3+i);
        u = vec4d_add(a, b);
        v = vec4d_sub(a, b);
        p = vec4d_mulmod2(d, fp, n, ninv);
        q = vec4d_mulmod2(d, fq, n, ninv);
        r = vec4d_mulmod2(d, fr, n, ninv);
        c = vec4d_reduce_to_pm1n(c, n, ninv);
        u = vec4d_reduce_to_pm1n(u, n, ninv);
        b = vec4d_mulmod2(v, f1, n, ninv);
        v = vec4d_mulmod2(v, f0, n, ninv);
        c = vec4d_fmsub(f2, c, v);
        a = vec4d_add(c, u);
        c = vec4d_sub(c, u);
        c = vec4d_mulmod2(c, f3, n, ninv);
        a = vec4d_add(a, p);
        b = vec4d_sub(b, q);
        c = vec4d_sub(c, r);
        vec4d_store(X0+i, a);
        vec4d_store(X1+i, b);
        vec4d_store(X2+i, c);
    } while (i += 4, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
    return 1;
}

/*
k = 2, n = 3, z = 3, f = true
[      -r + 1           r + 1         2]
[        2//w           -2//w         0]
[(r + 1)//w^2   (-r + 1)//w^2   -2//w^2]
[          -r               r         1]

    {x0, x1, x3, x4} = {        -r*(x0 - x1) + (x0 + x1) + 2*x2,
                        2*w^-1*(x0 - x1),
                         -w^-2*(-r*(x0 - x1) - (x0 + x1) + 2*x2),
                                -r*(x0 - x1)             +   x2  }
*/
static int radix_4_moth_inv_trunc_block_3_3_1(
    const sd_fft_ctx_t Q,
    ulong j, ulong jm,
    double* X0, double* X1, double* X2, double* X3)
{
    vec8d n    = vec8d_set_d(Q->p);
    vec8d ninv = vec8d_set_d(Q->pinv);
    double W  = UNLIKELY(j == 0) ? -Q->w2s[0] : Q->w2s[2*jm+1];
    double W2 = UNLIKELY(j == 0) ? -Q->w2s[0] : Q->w2s[jm];
    double twoW = vec1d_reduce_pm1n_to_pmhn(-2*W, Q->p);
    vec8d f0 = vec8d_set_d(Q->w2s[1]);  /* r */
    vec8d f1 = vec8d_set_d(twoW);       /* 2*w^-1 */
    vec8d f2 = vec8d_set_d(2);
    vec8d f3 = vec8d_set_d(W2);         /* -w^-2 */
    ulong i = 0; do {
        vec8d a, b, c, u, v;
        a = vec8d_load(X0+i);
        b = vec8d_load(X1+i);
        c = vec8d_load(X2+i);
        v = vec8d_sub(a, b);
        vec8d_store(X1+i, vec8d_mulmod2(v, f1, n, ninv));
        c = vec8d_reduce_to_pm1n(c, n, ninv);
        v = vec8d_mulmod2(v, f0, n, ninv);
        vec8d_store(X3+i, vec8d_sub(c, v));
        u = vec8d_reduce_to_pm1n(vec8d_add(a, b), n, ninv);
        c = vec8d_fmsub(f2, c, v);
        a = vec8d_add(c, u);
        c = vec8d_sub(c, u);
        c = vec8d_mulmod2(c, f3, n, ninv);
        vec8d_store(X0+i, a);
        vec8d_store(X2+i, c);
    } while (i += 8, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
    return 1;
}

/*
k = 2, n = 3, z = 3, f = false
[      -r + 1           r + 1         2]
[        2//w           -2//w         0]
[(r + 1)//w^2   (-r + 1)//w^2   -2//w^2]

    {x0, x1, x3} = {        -r*(x0 - x1) + (x0 + x1) + 2*x2,
                    2*w^-1*(x0 - x1),
                     -w^-2*(-r*(x0 - x1) - (x0 + x1) + 2*x2)}

                 = {        2*x2 - r*v0 + u0,
                    2*w^-1*v0,
                     -w^-2*(2*x2 - r*v0 - u0)}
*/
static int radix_4_moth_inv_trunc_block_3_3_0(
    const sd_fft_ctx_t Q,
    ulong j, ulong jm,
    double* X0, double* X1, double* X2, double* X3)
{
    vec8d n    = vec8d_set_d(Q->p);
    vec8d ninv = vec8d_set_d(Q->pinv);
    double W  = UNLIKELY(j == 0) ? -Q->w2s[0] : Q->w2s[2*jm+1];
    double W2 = UNLIKELY(j == 0) ? -Q->w2s[0] : Q->w2s[jm];
    double twoW = vec1d_reduce_pm1n_to_pmhn(-2*W, Q->p);
    vec8d f0 = vec8d_set_d(Q->w2s[1]);  /* r */
    vec8d f1 = vec8d_set_d(twoW);       /* 2*w^-1 */
    vec8d f2 = vec8d_set_d(2);
    vec8d f3 = vec8d_set_d(W2);         /* -w^-2 */
    ulong i = 0; do {
        vec8d a, b, c, u, v;
        a = vec8d_load(X0+i);
        b = vec8d_load(X1+i);
        c = vec8d_load(X2+i);
        v = vec8d_sub(a, b);
        vec8d_store(X1+i, vec8d_mulmod2(v, f1, n, ninv));
        c = vec8d_reduce_to_pm1n(c, n, ninv);
        v = vec8d_mulmod2(v, f0, n, ninv);
        u = vec8d_add(a, b);
        u = vec8d_reduce_to_pm1n(u, n, ninv);
        c = vec8d_fmsub(f2, c, v);
        a = vec8d_add(c, u);
        c = vec8d_sub(c, u);
        c = vec8d_mulmod2(c, f3, n, ninv);
        vec8d_store(X0+i, a);
        vec8d_store(X2+i, c);
    } while (i += 8, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
    return 1;
}

/*
k = 2, n = 2, z = 4, f = true
[            2                2        -w^2             0]
[         2//w            -2//w           0          -w^2]
[1//2*r + 1//2   -1//2*r + 1//2   -1//2*w^2   -1//2*r*w^3]
*/
static int radix_4_moth_inv_trunc_block_2_4_1(
    const sd_fft_ctx_t Q,
    ulong j, ulong jm,
    double* X0, double* X1, double* X2, double* X3)
{
    vec4d n    = vec4d_set_d(Q->p);
    vec4d ninv = vec4d_set_d(Q->pinv);
    double W = UNLIKELY(j == 0) ? -Q->w2s[0] : Q->w2s[2*jm+1];
    double twoW = vec1d_reduce_pm1n_to_pmhn(-2*W, Q->p);
    double rw3 = vec1d_mulmod2(Q->w2s[j], Q->w2s[2*j+1], Q->p, Q->pinv);
    vec4d f0 = vec4d_set_d(2);
    vec4d f1 = vec4d_set_d(twoW);                                   /* 2*w^-1 */
    vec4d f2 = vec4d_set_d(vec1d_fnmadd(0.5, Q->p, 0.5));           /* 1/2 */
    vec4d f3 = vec4d_set_d(Q->w2s[1]);                              /* r */
    vec4d f4 = vec4d_set_d(Q->w2s[j]);                              /* w^2 */
    vec4d f5 = vec4d_set_d(vec1d_reduce_pm1n_to_pmhn(rw3, Q->p));   /* r*w^3 */
    ulong i = 0; do {
        vec4d a, b, u, v, s, t, g, h, p, q, r;
        u = vec4d_load(X0+i);
        v = vec4d_load(X1+i);
        a = vec4d_load(X2+i);
        b = vec4d_load(X3+i);
        p = vec4d_mulmod2(a, f4, n, ninv);
        q = vec4d_mulmod2(b, f4, n, ninv);
        r = vec4d_mulmod2(b, f5, n, ninv);
        s = vec4d_add(u, v);
        s = vec4d_reduce_to_pm1n(s, n, ninv);
        t = vec4d_sub(u, v);
        g = vec4d_mulmod2(s, f0, n, ninv);
        h = vec4d_mulmod2(t, f1, n, ninv);
        t = vec4d_mulmod2(t, f3, n, ninv);
        vec4d_store(X0+i, vec4d_sub(g, p));
        vec4d_store(X1+i, vec4d_sub(h, q));
        u = vec4d_add(s, t);
        v = vec4d_add(p, r);
        u = vec4d_sub(u, v);
        u = vec4d_mulmod2(u, f2, n, ninv);
        vec4d_store(X2+i, u);
    } while (i += 4, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
    return 1;
}

/*
k = 2, n = 2, z = 4, f = false
[   2       2   -w^2      0]
[2//w   -2//w      0   -w^2]
*/
static int radix_4_moth_inv_trunc_block_2_4_0(
    const sd_fft_ctx_t Q,
    ulong j, ulong jm,
    double* X0, double* X1, double* X2, double* X3)
{
    vec8d n    = vec8d_set_d(Q->p);
    vec8d ninv = vec8d_set_d(Q->pinv);
    double W = UNLIKELY(j == 0) ? -Q->w2s[0] : Q->w2s[2*jm+1];
    vec8d w2 = vec8d_set_d(Q->w2s[j]);
    vec8d twowi = vec8d_set_d(vec1d_reduce_pm1n_to_pmhn(-2*W, Q->p));
    ulong i = 0; do {
        vec8d a, b, c, d, u, v;
        a = vec8d_load(X0+i);
        b = vec8d_load(X1+i);
        c = vec8d_load(X2+i);
        d = vec8d_load(X3+i);
        c = vec8d_mulmod2(c, w2, n, ninv);
        d = vec8d_mulmod2(d, w2, n, ninv);
        u = vec8d_add(a, b);
        v = vec8d_sub(a, b);
        u = vec8d_add(u, u);
        u = vec8d_reduce_to_pm1n(u, n, ninv);
        v = vec8d_mulmod2(v, twowi, n, ninv);
        u = vec8d_sub(u, c);
        v = vec8d_sub(v, d);
        vec8d_store(X0+i, u);
        vec8d_store(X1+i, v);
    } while (i += 8, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
    return 1;
}

/*
k = 2, n = 2, z = 2, f = true
[            2                2]
[         2//w            -2//w]
[1//2*r + 1//2   -1//2*r + 1//2]

{x0, x1, x2} = {2*(x0 + x1), 2*w^-1*(x0 - x1), (x0+x1)/2 + (x0-x1)*i/2}
*/
static int radix_4_moth_inv_trunc_block_2_2_1(
    const sd_fft_ctx_t Q,
    ulong j, ulong jm,
    double* X0, double* X1, double* X2, double* X3)
{
    vec8d n    = vec8d_set_d(Q->p);
    vec8d ninv = vec8d_set_d(Q->pinv);
    double W = (j == 0) ? -Q->w2s[0] : Q->w2s[2*jm+1];
    double ha = vec1d_fnmadd(0.5, Q->p, 0.5);
    double har = vec1d_mulmod2(Q->w2s[1], ha, Q->p, Q->pinv);
    vec8d c1 = vec8d_set_d(vec1d_reduce_pm1n_to_pmhn(-2*W, Q->p));  /* 2/w */
    vec8d c2 = vec8d_set_d(ha);                                     /* 1/2 */
    vec8d c3 = vec8d_set_d(vec1d_reduce_pm1n_to_pmhn(har, Q->p));   /* r/2 */
    ulong i = 0; do {
        vec8d u, v, s, t;
        u = vec8d_load(X0+i);
        v = vec8d_load(X1+i);
        s = vec8d_add(u, v);
        t = vec8d_sub(u, v);
        u = vec8d_add(s, s);
        u = vec8d_reduce_to_pm1n(u, n, ninv);
        v = vec8d_mulmod2(t, c1, n, ninv);
        s = vec8d_mulmod2(s, c2, n, ninv);
        t = vec8d_mulmod2(t, c3, n, ninv);
        s = vec8d_add(s, t);
        vec8d_store(X0+i, u);
        vec8d_store(X1+i, v);
        vec8d_store(X2+i, s);
    } while (i += 8, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
    return 1;
}

/*
k = 2, n = 2, z = 2, f = 0

{x0, x1} = {2*(x0 + x1), 2*w^-1*(x0 - x1)}
*/
static int radix_4_moth_inv_trunc_block_2_2_0(
    const sd_fft_ctx_t Q,
    ulong j, ulong jm,
    double* X0, double* X1, double* X2, double* X3)
{
    vec8d n    = vec8d_set_d(Q->p);
    vec8d ninv = vec8d_set_d(Q->pinv);
    double W = (j == 0) ? -Q->w2s[0] : Q->w2s[2*jm+1];
    vec8d c0 = vec8d_set_d(2);
    vec8d c1 = vec8d_set_d(vec1d_reduce_pm1n_to_pmhn(-2*W, Q->p));
    ulong i = 0; do {
        vec8d u, v, s, t;
        u = vec8d_load(X0+i);
        v = vec8d_load(X1+i);
        s = vec8d_add(u, v);
        t = vec8d_sub(u, v);
        u = vec8d_mulmod2(s, c0, n, ninv);
        v = vec8d_mulmod2(t, c1, n, ninv);
        vec8d_store(X0+i, u);
        vec8d_store(X1+i, v);
    } while (i += 8, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
    return 1;
}


/*
k = 2, n = 1, z = 4, f = true
[4        -w   -w^2        -w^3]
[1   -1//2*w      0   -1//2*w^3]
*/
static int radix_4_moth_inv_trunc_block_1_4_1(
    const sd_fft_ctx_t Q,
    ulong j, ulong jm,
    double* X0, double* X1, double* X2, double* X3)
{
    vec4d n    = vec4d_set_d(Q->p);
    vec4d ninv = vec4d_set_d(Q->pinv);
    vec4d f2 = vec4d_set_d(2);
    vec4d w2 = vec4d_set_d(Q->w2s[j]);
    double ha = vec1d_fnmadd(0.5, Q->p, 0.5);
    vec4d wo2 = vec4d_set_d(vec1d_reduce_pm1n_to_pmhn(vec1d_mulmod2(Q->w2s[2*j], ha, Q->p, Q->pinv), Q->p));
    ulong i = 0; do {
        vec4d a, b, c, d, u;
        a = vec4d_load(X0+i);
        a = vec4d_reduce_to_pm1n(a, n, ninv);
        b = vec4d_load(X1+i);
        c = vec4d_load(X2+i);
        d = vec4d_load(X3+i);
        c = vec4d_mulmod2(c, w2, n, ninv);
        d = vec4d_mulmod2(d, w2, n, ninv);
        b = vec4d_add(b, d);
        b = vec4d_mulmod2(b, wo2, n, ninv);
        u = vec4d_fmsub(f2, a, b);
        b = vec4d_sub(a, b);
        a = vec4d_reduce_to_pm1n(vec4d_add(u, u), n, ninv);
        a = vec4d_sub(a, c);
        vec4d_store(X0+i, a);
        vec4d_store(X1+i, b);
    } while (i += 4, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
    return 1;
}

/*
k = 2, n = 1, z = 4, f = false
[4   -w   -w^2   -w^3]
*/
static int radix_4_moth_inv_trunc_block_1_4_0(
    const sd_fft_ctx_t Q,
    ulong j, ulong jm,
    double* X0, double* X1, double* X2, double* X3)
{
    vec8d n    = vec8d_set_d(Q->p);
    vec8d ninv = vec8d_set_d(Q->pinv);
    vec8d f1 = vec8d_set_d(4.0);
    vec8d w2 = vec8d_set_d(Q->w2s[j]);
    vec8d w  = vec8d_set_d(Q->w2s[2*j]);
    ulong i = 0; do {
        vec8d a, b, c, d;
        a = vec8d_load(X0+i);
        b = vec8d_load(X1+i);
        c = vec8d_load(X2+i);
        d = vec8d_load(X3+i);
        a = vec8d_mul(a, f1);
        d = vec8d_mulmod2(d, w, n, ninv);
        a = vec8d_reduce_to_pm1n(a, n, ninv);
        b = vec8d_mulmod2(b, w, n, ninv);
        a = vec8d_sub(a, b);
        c = vec8d_add(c, d);
        c = vec8d_mulmod2(c, w2, n, ninv);
        a = vec8d_sub(a, c);
        vec8d_store(X0+i, a);
    } while (i += 8, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
    return 1;
}

/*
k = 2, n = 1, z = 1, f = true
[4]
[1]
*/
static int radix_4_moth_inv_trunc_block_1_1_1(
    const sd_fft_ctx_t Q,
    ulong j, ulong jm,
    double* X0, double* X1, double* X2, double* X3)
{
    vec8d n    = vec8d_set_d(Q->p);
    vec8d ninv = vec8d_set_d(Q->pinv);
    vec8d f = vec8d_set_d(4.0);
    ulong i = 0; do {
        vec8d a, b;
        a = vec8d_load(X0+i);
        b = vec8d_mul(f, a);
        vec8d_store(X0+i, vec8d_reduce_to_pm1n(b, n, ninv));
        vec8d_store(X1+i, vec8d_reduce_to_pm1n(a, n, ninv));
    } while (i += 8, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
    return 1;
}

/*
k = 2, n = 1, z = 1, f = false
[4]
*/
static int radix_4_moth_inv_trunc_block_1_1_0(
    const sd_fft_ctx_t Q,
    ulong j, ulong jm,
    double* X0, double* X1, double* X2, double* X3)
{
    vec8d n    = vec8d_set_d(Q->p);
    vec8d ninv = vec8d_set_d(Q->pinv);
    vec8d f = vec8d_set_d(4.0);
    ulong i = 0; do {
        vec8d a, b;
        a = vec8d_load(X0+i);
        b = vec8d_mul(f, a);
        vec8d_store(X0+i, vec8d_reduce_to_pm1n(b, n, ninv));
    } while (i += 8, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
    return 1;
}

/*
k = 2, n = 0, z = 4, f = true
[1//4   1//4*w   1//4*w^2   1//4*w^3]
*/
static int radix_4_moth_inv_trunc_block_0_4_1(
    const sd_fft_ctx_t Q,
    ulong j, ulong jm,
    double* X0, double* X1, double* X2, double* X3)
{
    vec8d n    = vec8d_set_d(Q->p);
    vec8d ninv = vec8d_set_d(Q->pinv);
    vec8d one4th = vec8d_set_d(vec1d_fnmadd(0.25, Q->p, 0.25));
    vec8d w = vec8d_set_d(Q->w2s[2*j]);
    vec8d w2 = vec8d_set_d(Q->w2s[j]);
    ulong i = 0; do {
        vec8d a, b, c, d;
        a = vec8d_load(X0+i);
        b = vec8d_load(X1+i);
        c = vec8d_load(X2+i);
        d = vec8d_load(X3+i);
        b = vec8d_mulmod2(b, w, n, ninv);
        d = vec8d_mulmod2(d, w, n, ninv);
        a = vec8d_add(a, b);
        c = vec8d_add(c, d);
        c = vec8d_mulmod2(c, w2, n, ninv);
        a = vec8d_add(a, c);
        vec8d_store(X0+i, vec8d_mulmod2(a, one4th, n, ninv));
    } while (i += 8, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
    return 1;
}

/************************ the recursive stuff ********************************/

void sd_ifft_main_block(
    const sd_fft_ctx_t Q,
    double* d,
    ulong I, /* starting index */
    ulong S, /* stride */
    ulong k, /* BLK_SZ transforms each of length 2^k */
    ulong j)
{
    if (k > 2)
    {
        ulong k1 = k/2;
        ulong k2 = k - k1;

        /* row ffts */
        ulong l1 = n_pow2(k1);
        ulong b = 0; do {
            sd_ifft_main_block(Q, d, I + (b<<k2)*S, S, k2, (j<<k1) + b);
        } while (b++, b < l1);

        /* column ffts */
        ulong l2 = n_pow2(k2);
        ulong a = 0; do {
            sd_ifft_main_block(Q, d, I + a*S, S<<k2, k1, j);
        } while (a++, a < l2);

        return;
    }

    ulong jm = j^(n_next_pow2m1(j)>>1);

    if (k == 2)
    {
        double* X0 = sd_fft_ctx_blk_index(Q, d, I + S*0);
        double* X1 = sd_fft_ctx_blk_index(Q, d, I + S*1);
        double* X2 = sd_fft_ctx_blk_index(Q, d, I + S*2);
        double* X3 = sd_fft_ctx_blk_index(Q, d, I + S*3);
        RADIX_4_REVERSE_PARAM(vec8d, Q, j, jm, 1)
        ulong i = 0; do {
            RADIX_4_REVERSE_MOTH(vec8d, X0+i, X1+i, X2+i, X3+i);
        } while(i += 8, i < BLK_SZ);
    }
    else if (k == 1)
    {
        double* X0 = sd_fft_ctx_blk_index(Q, d, I + S*0);
        double* X1 = sd_fft_ctx_blk_index(Q, d, I + S*1);
        RADIX_2_REVERSE_PARAM(vec8d, Q, j)
        ulong i = 0; do {
            RADIX_2_REVERSE_MOTH(vec8d, X0+i, X1+i);
        } while (i += 8, i < BLK_SZ);
    }
}

void sd_ifft_main(
    const sd_fft_ctx_t Q,
    double* d,
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
            sd_ifft_main(Q, d, I + b*(S<<k2), S, k2, (j<<k1) + b);
        } while (b++, b < l1);

        ulong l2 = n_pow2(k2);
        ulong a = 0; do {
            sd_ifft_main_block(Q, d, I + a*S, S<<k2, k1, j);
        } while (a++, a < l2);

        return;
    }

    if (k == 2)
    {
        /* k1 = 2; k2 = 0 */
        sd_ifft_base_1(Q, d, I + S*0, 4*j+0);
        sd_ifft_base_0(Q, d, I + S*1, 4*j+1);
        sd_ifft_base_0(Q, d, I + S*2, 4*j+2);
        sd_ifft_base_0(Q, d, I + S*3, 4*j+3);
        sd_ifft_main_block(Q, d, I, S, 2, j);
    }
    else if (k == 1)
    {
        /* k1 = 1; k2 = 0 */
        sd_ifft_base_1(Q, d, I + S*0, 2*j+0);
        sd_ifft_base_0(Q, d, I + S*1, 2*j+1);
        sd_ifft_main_block(Q, d, I, S, 1, j);
    }
    else
    {
        sd_ifft_base_1(Q, d, I, j);
    }
}

void sd_ifft_trunc_block(
    const sd_fft_ctx_t Q,
    double* d,
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
        sd_ifft_main_block(Q, d, I, S, k, j);
        return;
    }

    if (k == 2)
    {
#define IT(nn, zz, ff) CAT4(radix_4_moth_inv_trunc_block, nn, zz, ff)
#define LOOKUP_IT(nn, zz, ff) tab[(ulong)(ff) + 2*((zz)-1 + 4*(nn))]
        static int (*tab[5*4*2])(const sd_fft_ctx_t, ulong, ulong, double*, double*, double*, double*) =
            {IT(0,1,0),IT(0,1,1), IT(0,2,0),IT(0,2,1), IT(0,3,0),IT(0,3,1), IT(0,4,0),IT(0,4,1),
             IT(1,1,0),IT(1,1,1), IT(1,2,0),IT(1,2,1), IT(1,3,0),IT(1,3,1), IT(1,4,0),IT(1,4,1),
             IT(2,1,0),IT(2,1,1), IT(2,2,0),IT(2,2,1), IT(2,3,0),IT(2,3,1), IT(2,4,0),IT(2,4,1),
             IT(3,1,0),IT(3,1,1), IT(3,2,0),IT(3,2,1), IT(3,3,0),IT(3,3,1), IT(3,4,0),IT(3,4,1),
             IT(4,1,0),IT(4,1,1), IT(4,2,0),IT(4,2,1), IT(4,3,0),IT(4,3,1), IT(4,4,0),IT(4,4,1)};

        ulong jm = j^(n_next_pow2m1(j)>>1);
        if (LOOKUP_IT(n,z,f)(Q, j, jm, sd_fft_ctx_blk_index(Q, d, I + S*0),
                                       sd_fft_ctx_blk_index(Q, d, I + S*1),
                                       sd_fft_ctx_blk_index(Q, d, I + S*2),
                                       sd_fft_ctx_blk_index(Q, d, I + S*3)))
        {
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
            sd_ifft_main_block(Q, d, I + b*(S << k2), S, k2, (j << k1) + b);

        /* rightmost columns */
        for (ulong a = n2; a < z2p; a++)
            sd_ifft_trunc_block(Q, d, I + a*S, S << k2, k1, j, z1 + (a < mp), n1, fp);

        /* last partial row */
        if (fp)
            sd_ifft_trunc_block(Q, d, I + n1*(S << k2), S, k2, (j << k1) + n1, z2p, n2, f);

        /* leftmost columns */
        for (ulong a = 0; a < n2; a++)
            sd_ifft_trunc_block(Q, d, I + a*S, S << k2, k1, j, z1 + (a < m), n1 + 1, 0);

        return;
    }

    if (k == 1)
    {
#define IT(nn, zz, ff) CAT4(radix_2_moth_inv_trunc_block, nn, zz, ff)
#define LOOKUP_IT(nn, zz, ff) tab[(ulong)(ff) + 2*((zz)-1 + 2*(nn))]
        static void (*tab[3*2*2])(const sd_fft_ctx_t, ulong, double*, double*) =
            {IT(0,1,0),IT(0,1,1), IT(0,2,0),IT(0,2,1),
             IT(1,1,0),IT(1,1,1), IT(1,2,0),IT(1,2,1),
             IT(2,1,0),IT(2,1,1), IT(2,2,0),IT(2,2,1)};

        LOOKUP_IT(n,z,f)(Q, j, sd_fft_ctx_blk_index(Q, d, I + S*0),
                               sd_fft_ctx_blk_index(Q, d, I + S*1));
        return;
#undef LOOKUP_IT
#undef IT
    }
}


void sd_ifft_trunc(
    const sd_fft_ctx_t Q,
    double* d,
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
            sd_ifft_main(Q, d, I + b*(S << k2), S, k2, (j << k1) + b);

        /* rightmost columns */
        for (ulong a = n2; a < z2p; a++)
            sd_ifft_trunc_block(Q, d, I + a*S, S << k2, k1, j, z1 + (a < mp), n1, fp);

        /* last partial row */
        if (fp)
            sd_ifft_trunc(Q, d, I + n1*(S << k2), S, k2, (j << k1) + n1, z2p, n2, f);

        /* leftmost columns */
        for (ulong a = 0; a < n2; a++)
            sd_ifft_trunc_block(Q, d, I + a*S, S << k2, k1, j, z1 + (a < m), n1 + 1, 0);

        return;
    }

    if (k == 2)
    {
                   sd_ifft_base_1(Q, d, I + S*0, 4*j+0);
        if (n > 1) sd_ifft_base_0(Q, d, I + S*1, 4*j+1);
        if (n > 2) sd_ifft_base_0(Q, d, I + S*2, 4*j+2);
        if (n > 3) sd_ifft_base_0(Q, d, I + S*3, 4*j+3);
        sd_ifft_trunc_block(Q, d, I, S, 2, j, z, n, f);
        if (f) sd_ifft_trunc(Q, d, I + S*n, S, 0, 4*j+n, 1, 0, f);
    }
    else if (k == 1)
    {
                   sd_ifft_base_1(Q, d, I + S*0, 2*j+0);
        if (n > 1) sd_ifft_base_0(Q, d, I + S*1, 2*j+1);
        sd_ifft_trunc_block(Q, d, I, S, 1, j, z, n, f);
        if (f) sd_ifft_trunc(Q, d, I + S*n, S, 0, 2*j+n, 1, 0, f);
    }
    else
    {
        FLINT_ASSERT(!f);
        sd_ifft_base_1(Q, d, I, j);
    }
}

