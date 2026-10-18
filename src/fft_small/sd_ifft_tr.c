/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fft_small.h"
#include "butterflies.h"

/*
    This module contains the function `sd_ifft_trunc_tr`, the transpose of the
    truncated inverse FFT `sd_ifft_trunc`.

    `sd_ifft_trunc(Q, d, L, trunc)` is a linear map on the first
    `BLK_SZ*ceil(trunc/BLK_SZ)` entries of `d` (the first `2^L` entries when
    `L <= LG_BLK_SZ`); `sd_ifft_trunc_tr` computes its transpose, in place, so
    that for all `a` and `v`

        <sd_ifft_trunc(a), v> == <a, sd_ifft_trunc_tr(v)>.

    Together with `sd_fft_trunc_tr` this gives the transpose of a truncated
    convolution: the transpose of `x -> ifft(fft(x) * fft(y))` with respect to
    `x` is `v -> fft_tr(fft(y) * ifft_tr(v))`.

    ********************************************************************************

    As in `sd_fft_tr.c`, this file is the mechanical transposition of
    `sd_ifft.c`: the order of the steps is reversed and each step is replaced
    by its transpose. See the comment at the top of `sd_fft_tr.c` for the
    general recipe. Two things are specific to this file:

    - the transpose of the decimation-in-frequency butterfly
      `(a0, a1) -> (a0 + a1, w^-1*(a0 - a1))` used by `sd_ifft.c` is the
      decimation-in-time butterfly `(a0, a1) -> (a0 + w^-1*a1, a0 - w^-1*a1)`,
      i.e. the shape used by `sd_fft.c`, with the inverse twiddles of
      `sd_ifft.c`. So the twiddle lookups here are those of `sd_ifft.c`
      (`SET_J_BITS_AND_J_MR`), and the butterflies are those of `sd_fft.c`
      (see butterflies.h), with `w = -W` given by multiplying with
      `nmulmod_fast`, hence with the same reduction placement and bounds;

    - the truncated inverse butterflies of `sd_ifft.c` are not restrictions of
      one formula, they are 13 + 6 separately derived maps, each documented by
      its matrix in a comment there. Each of them is transposed here
      separately, and the matrix in each comment below is the transpose of the
      matrix in the corresponding comment of `sd_ifft.c`. The naming and the
      `(n, z, f)` parameters are unchanged, so
      `radix_4_moth_tr_trunc_block_{n}_{z}_{f}` reads `n + f` blocks and writes
      `z` blocks, i.e. the other way round from
      `radix_4_moth_inv_trunc_block_{n}_{z}_{f}`, and the dispatch tables and
      the recursion index bookkeeping are identical to those of `sd_ifft.c`.

    Throughout, `r` denotes a square root of `-1`, `w` denotes the twiddle
    `sd_fft_ctx_w(Q, 2*j)` of the two children of node `j`, so that `w^2` is
    the twiddle `sd_fft_ctx_w2(Q, j)` of node `j` itself and `r*w` is the
    twiddle of the second child, exactly as in `sd_ifft.c`.

    `sd_ifft_trunc_tr` performs the same number of operations as
    `sd_ifft_trunc`, and costs the same. Its input and output ranges are those
    of `sd_fft_trunc`, whose butterfly shape it shares: entries in `(-3n, 3n)`
    are accepted and produced. Like `sd_fft_trunc` it
    reads the coefficient side of the transform, where entries are usually
    reduced residues in `[0, n)`.
*/

#define N 8
#define VECND vec8d
#define VECNOP(op) vec8d_##op

#define M 4
#define VECMD vec4d
#define VECMOP(op) vec4d_##op

/***************** transposed inverse butterfly *******************************
    The transpose of
        2*a0 =      (b0 + b1)
        2*a1 = w^-1*(b0 - b1)
    is
        b0 = a0 + w^-1*a1
        b1 = a0 - w^-1*a1
    with, as in sd_ifft.c, W := -w^-1: the DIT butterfly of sd_fft.c with
    twiddle -W, see butterflies.h.
*/
#define RADIX_2_REVERSE_TR_PARAM_J_IS_Z(V, Q) \
    V n    = CAT(V, set_d)(Q->p); \
    V ninv = CAT(V, set_d)(Q->pinv);

#define RADIX_2_REVERSE_TR_MOTH_J_IS_Z(V, X0, X1) \
    SD_DIT2_MOTH_J_IS_Z(V, X0, X1)

#define RADIX_2_REVERSE_TR_PARAM_J_IS_NZ(V, Q, j_mr, j_bits) \
    V W = CAT(V, set_d)(Q->w2tab[j_bits][j_mr]); \
    V n    = CAT(V, set_d)(Q->p); \
    V ninv = CAT(V, set_d)(Q->pinv);

#define RADIX_2_REVERSE_TR_MOTH_J_IS_NZ(V, X0, X1) \
    SD_DIT2_MOTH_J_IS_NZ(V, nmulmod_fast, W, X0, X1)

/***************** transposed inverse butterfly *******************************
    The transpose of
        4*a0 =            (b0 + b1) +        (b2 + b3)
        4*a1 =       w^-1*(b0 - b1) - i*w^-1*(b2 - b3)
        4*a2 = w^-2*(     (b0 + b1) -        (b2 + b3))
        4*a3 = w^-2*(w^-1*(b0 - b1) + i*w^-1*(b2 - b3))
    is
        b0 = (a0 - W2*a2) - W *(a1 - W2*a3)
        b1 = (a0 - W2*a2) + W *(a1 - W2*a3)
        b2 = (a0 + W2*a2) - IW*(a1 + W2*a3)
        b3 = (a0 + W2*a2) + IW*(a1 + W2*a3)
    with, as in sd_ifft.c, W := -w^-1, W2 := -w^-2, IW := i*w^-1: the radix-4
    butterfly of sd_fft.c with twiddles (-W2, -W, -IW).
*/
#define RADIX_4_REVERSE_TR_PARAM_J_IS_Z(V, Q) \
    V IW = CAT(V, set_d)(Q->w2tab[0][1]); \
    V n    = CAT(V, set_d)(Q->p); \
    V ninv = CAT(V, set_d)(Q->pinv); \

#define RADIX_4_REVERSE_TR_MOTH_J_IS_Z(V, X0, X1, X2, X3) \
    SD_DIT4_MOTH_J_IS_Z(V, nmulmod_fast, IW, X0, X1, X2, X3)

#define RADIX_4_REVERSE_TR_PARAM_J_IS_NZ(V, Q, j_mr, j_bits) \
    V W  = CAT(V, set_d)(Q->w2tab[1+j_bits][2*j_mr+1]); \
    V W2 = CAT(V, set_d)(Q->w2tab[0+j_bits][j_mr]); \
    V IW = CAT(V, set_d)(Q->w2tab[1+j_bits][2*j_mr+0]); \
    V n    = CAT(V, set_d)(Q->p); \
    V ninv = CAT(V, set_d)(Q->pinv); \

#define RADIX_4_REVERSE_TR_MOTH_J_IS_NZ(V, X0, X1, X2, X3) \
    SD_DIT4_MOTH_J_IS_NZ(V, nmulmod_fast, W2, W, IW, X0, X1, X2, X3)

#define _RADIX_4_REVERSE_TR_PARAM_J_IS_Z(...)  RADIX_4_REVERSE_TR_PARAM_J_IS_Z(__VA_ARGS__)
#define _RADIX_4_REVERSE_TR_MOTH_J_IS_Z(...)   RADIX_4_REVERSE_TR_MOTH_J_IS_Z(__VA_ARGS__)
#define _RADIX_4_REVERSE_TR_PARAM_J_IS_NZ(...) RADIX_4_REVERSE_TR_PARAM_J_IS_NZ(__VA_ARGS__)
#define _RADIX_4_REVERSE_TR_MOTH_J_IS_NZ(...)  RADIX_4_REVERSE_TR_MOTH_J_IS_NZ(__VA_ARGS__)

/* transposes of the LENGTH{2,4,8}INV macros of sd_ifft.c: the DIT register
   basecases of sd_fft.c with negated twiddles. For j != 0 the 3-layer one is
   the LENGTH8_ANY_J of sd_fft.c, which performs the same operations as
   composing LENGTH4INVTR_ANY_J and LENGTH2INVTR_ANY_J (as done below for
   j == 0), in another order. */

#define LENGTH2INVTR_ANY_J(T, ...) SD_DIT_LENGTH2_ANY_J(T, nmulmod_fast, __VA_ARGS__)

#define LENGTH2INVTR_ZERO_J(T, ...) SD_DIT_LENGTH2_ZERO_J(T, __VA_ARGS__)

#define LENGTH4INVTR_ANY_J(T, ...) SD_DIT_LENGTH4_ANY_J(T, nmulmod_fast, __VA_ARGS__)

#define LENGTH4INVTR_ZERO_J(T, ...) SD_DIT_LENGTH4_ZERO_J(T, nmulmod_fast, __VA_ARGS__)

#define LENGTH8INVTR_ANY_J(T, ...) SD_DIT_LENGTH8_ANY_J(T, nmulmod_fast, __VA_ARGS__)

#define LENGTH8INVTR_ZERO_J(T, x0, x1, x2, x3, x4, x5, x6, x7, n, ninv, ww1, www2, www3) \
{ \
    T A0 = x0, A1 = x1, A2 = x2, A3 = x3, A4 = x4, A5 = x5, A6 = x6, A7 = x7; \
    LENGTH4INVTR_ZERO_J(T, A0,A2,A4,A6, n,ninv, ww1) \
    LENGTH4INVTR_ZERO_J(T, A1,A3,A5,A7, n,ninv, ww1) \
    LENGTH2INVTR_ZERO_J(T, A0,A1, n,ninv) \
    LENGTH2INVTR_ANY_J(T, A2,A3, n,ninv, ww1) \
    LENGTH2INVTR_ANY_J(T, A4,A5, n,ninv, www2) \
    LENGTH2INVTR_ANY_J(T, A6,A7, n,ninv, www3) \
    x0 = A0; \
    x1 = A1; \
    x2 = A2; \
    x3 = A3; \
    x4 = A4; \
    x5 = A5; \
    x6 = A6; \
    x7 = A7; \
}

/****** transposed basecase inverse transform of size dividing BLK_SZ *********/

static void sd_ifft_tr_basecase_0_1(const sd_fft_ctx_t FLINT_UNUSED(Q), double* FLINT_UNUSED(X))
{
}


static void sd_ifft_tr_basecase_1_1(const sd_fft_ctx_t Q, double* X)
{
    LENGTH2INVTR_ZERO_J(vec1d, X[0], X[1], Q->p, Q->pinv);
}


static void sd_ifft_tr_basecase_2_1(const sd_fft_ctx_t Q, double* X)
{
    LENGTH4INVTR_ZERO_J(vec1d, X[0], X[1], X[2], X[3], Q->p, Q->pinv, Q->w2tab[1][0]);
}


static void sd_ifft_tr_basecase_3_1(const sd_fft_ctx_t Q, double* X)
{
    LENGTH8INVTR_ZERO_J(vec1d, X[0], X[1], X[2], X[3], X[4], X[5], X[6], X[7], Q->p, Q->pinv,
                       Q->w2tab[1][0], Q->w2tab[2][1], Q->w2tab[2][0]);
}


static void sd_ifft_tr_basecase_4_1(const sd_fft_ctx_t Q, double* X)
{
    vec4d n    = vec4d_set_d(Q->p);
    vec4d ninv = vec4d_set_d(Q->pinv);
    vec4d W, W2, IW;

    vec4d x0 = vec4d_load(X+4*0);
    vec4d x1 = vec4d_load(X+4*1);
    vec4d x2 = vec4d_load(X+4*2);
    vec4d x3 = vec4d_load(X+4*3);

    IW = vec4d_set_d(Q->w2tab[0][1]);
    LENGTH4INVTR_ZERO_J(vec4d, x0,x1,x2,x3, n,ninv, IW);
    VEC4D_TRANSPOSE(x0,x1,x2,x3, x0,x1,x2,x3);

    W  = vec4d_set_d4(-Q->w2tab[0][0], Q->w2tab[0][3], Q->w2tab[0][7], Q->w2tab[0][5]);
    IW = vec4d_set_d4( Q->w2tab[0][1], Q->w2tab[0][2], Q->w2tab[0][6], Q->w2tab[0][4]);
    W2 = vec4d_set_d4(-Q->w2tab[0][0], Q->w2tab[0][1], Q->w2tab[0][3], Q->w2tab[0][2]);
    LENGTH4INVTR_ANY_J(vec4d, x0,x1,x2,x3, n,ninv, W2, W, IW);

    vec4d_store(X+4*0, x0);
    vec4d_store(X+4*1, x1);
    vec4d_store(X+4*2, x2);
    vec4d_store(X+4*3, x3);
}

static void sd_ifft_tr_basecase_4_0(const sd_fft_ctx_t Q, double* X, ulong j_mr, ulong j_bits)
{
    vec4d n    = vec4d_set_d(Q->p);
    vec4d ninv = vec4d_set_d(Q->pinv);
    vec4d W, W2, IW, u, v;

    vec4d x0 = vec4d_load(X+4*0);
    vec4d x1 = vec4d_load(X+4*1);
    vec4d x2 = vec4d_load(X+4*2);
    vec4d x3 = vec4d_load(X+4*3);

    W  = vec4d_set_d(Q->w2tab[1+j_bits][2*j_mr+1]);
    IW = vec4d_set_d(Q->w2tab[1+j_bits][2*j_mr+0]);
    W2 = vec4d_set_d(Q->w2tab[0+j_bits][1*j_mr+0]);
    LENGTH4INVTR_ANY_J(vec4d, x0,x1,x2,x3, n,ninv, W2, W, IW);
    VEC4D_TRANSPOSE(x0,x1,x2,x3, x0,x1,x2,x3);

    W2 = vec4d_load_aligned(&Q->w2tab[2+j_bits][4*j_mr+0]); /* a b c d */
    W2 = vec4d_permute_3_2_1_0(W2);                         /* d c b a */
    u  = vec4d_load_aligned(&Q->w2tab[3+j_bits][8*j_mr+0]); /* 0 1 2 3 */
    v  = vec4d_load_aligned(&Q->w2tab[3+j_bits][8*j_mr+4]); /* 4 5 6 7 */
    W  = vec4d_unpackhi_permute_3_1_2_0(u, v);              /* 7 5 3 1 */
    IW = vec4d_unpacklo_permute_3_1_2_0(u, v);              /* 6 4 2 0 */
    LENGTH4INVTR_ANY_J(vec4d, x0,x1,x2,x3, n,ninv, W2, W, IW);

    vec4d_store(X+4*0, x0);
    vec4d_store(X+4*1, x1);
    vec4d_store(X+4*2, x2);
    vec4d_store(X+4*3, x3);
}


static void sd_ifft_tr_basecase_5_1(const sd_fft_ctx_t Q, double* X)
{
    vec4d n    = vec4d_set_d(Q->p);
    vec4d ninv = vec4d_set_d(Q->pinv);
    vec4d w0, ww0, ww1, www2, www3;

    vec4d x0 = vec4d_load(X+4*0);
    vec4d x1 = vec4d_load(X+4*1);
    vec4d x2 = vec4d_load(X+4*2);
    vec4d x3 = vec4d_load(X+4*3);
    vec4d x4 = vec4d_load(X+4*4);
    vec4d x5 = vec4d_load(X+4*5);
    vec4d x6 = vec4d_load(X+4*6);
    vec4d x7 = vec4d_load(X+4*7);

    ww1  = vec4d_set_d(Q->w2tab[1][0]);
    www2 = vec4d_set_d(Q->w2tab[2][1]);
    www3 = vec4d_set_d(Q->w2tab[2][0]);
    LENGTH8INVTR_ZERO_J(vec4d, x0,x1,x2,x3,x4,x5,x6,x7, n,ninv, ww1, www2,www3);

    VEC4D_TRANSPOSE(x0,x1,x2,x3, x0,x1,x2,x3);
    VEC4D_TRANSPOSE(x4,x5,x6,x7, x4,x5,x6,x7);

    /* j = 0, 1, 2, 3  then {j,2j+0,2j+1}^-1 in each column */
    w0  = vec4d_set_d4(-Q->w2tab[0][0], Q->w2tab[1][0], Q->w2tab[2][1], Q->w2tab[2][0]);
    ww0 = vec4d_set_d4(-Q->w2tab[0][0], Q->w2tab[2][1], Q->w2tab[3][3], Q->w2tab[3][1]);
    ww1 = vec4d_set_d4( Q->w2tab[1][0], Q->w2tab[2][0], Q->w2tab[3][2], Q->w2tab[3][0]);
    LENGTH4INVTR_ANY_J(vec4d, x0,x1,x2,x3, n,ninv, w0, ww0, ww1);

    /* j = 4, 5, 6, 7  then {j,2j+0,2j+1}^-1 in each column */
    w0  = vec4d_set_d4(Q->w2tab[3][3], Q->w2tab[3][2], Q->w2tab[3][1], Q->w2tab[3][0]);
    ww0 = vec4d_set_d4(Q->w2tab[4][7], Q->w2tab[4][5], Q->w2tab[4][3], Q->w2tab[4][1]);
    ww1 = vec4d_set_d4(Q->w2tab[4][6], Q->w2tab[4][4], Q->w2tab[4][2], Q->w2tab[4][0]);
    LENGTH4INVTR_ANY_J(vec4d, x4,x5,x6,x7, n,ninv, w0, ww0, ww1);

    vec4d_store(X+4*0, x0);
    vec4d_store(X+4*1, x1);
    vec4d_store(X+4*2, x2);
    vec4d_store(X+4*3, x3);
    vec4d_store(X+4*4, x4);
    vec4d_store(X+4*5, x5);
    vec4d_store(X+4*6, x6);
    vec4d_store(X+4*7, x7);
}

static void sd_ifft_tr_basecase_5_0(const sd_fft_ctx_t Q, double* X, ulong j_mr, ulong j_bits)
{
    vec4d n    = vec4d_set_d(Q->p);
    vec4d ninv = vec4d_set_d(Q->pinv);
    vec4d u, v, w0, ww0, ww1, www0, www1, www2, www3;

    vec4d x0 = vec4d_load(X+4*0);
    vec4d x1 = vec4d_load(X+4*1);
    vec4d x2 = vec4d_load(X+4*2);
    vec4d x3 = vec4d_load(X+4*3);
    vec4d x4 = vec4d_load(X+4*4);
    vec4d x5 = vec4d_load(X+4*5);
    vec4d x6 = vec4d_load(X+4*6);
    vec4d x7 = vec4d_load(X+4*7);

    w0   = vec4d_set_d(Q->w2tab[0+j_bits][1*j_mr+0]);
    ww0  = vec4d_set_d(Q->w2tab[1+j_bits][2*j_mr+1]);
    ww1  = vec4d_set_d(Q->w2tab[1+j_bits][2*j_mr+0]);
    www0 = vec4d_set_d(Q->w2tab[2+j_bits][4*j_mr+3]);
    www1 = vec4d_set_d(Q->w2tab[2+j_bits][4*j_mr+2]);
    www2 = vec4d_set_d(Q->w2tab[2+j_bits][4*j_mr+1]);
    www3 = vec4d_set_d(Q->w2tab[2+j_bits][4*j_mr+0]);
    LENGTH8INVTR_ANY_J(vec4d, x0,x1,x2,x3,x4,x5,x6,x7, n,ninv, w0, ww0,ww1, www0,www1,www2,www3);

    VEC4D_TRANSPOSE(x0,x1,x2,x3, x0,x1,x2,x3);
    VEC4D_TRANSPOSE(x4,x5,x6,x7, x4,x5,x6,x7);

    w0 = vec4d_load_aligned(&Q->w2tab[3+j_bits][8*j_mr+4]);
    w0 = vec4d_permute_3_2_1_0(w0);
    u  = vec4d_load_aligned(&Q->w2tab[4+j_bits][16*j_mr+8]);
    v  = vec4d_load_aligned(&Q->w2tab[4+j_bits][16*j_mr+12]);
    ww0 = vec4d_unpackhi_permute_3_1_2_0(u, v);
    ww1 = vec4d_unpacklo_permute_3_1_2_0(u, v);
    LENGTH4INVTR_ANY_J(vec4d, x0,x1,x2,x3, n,ninv, w0,ww0,ww1);

    w0 = vec4d_load_aligned(&Q->w2tab[3+j_bits][8*j_mr+0]);
    w0 = vec4d_permute_3_2_1_0(w0);
    u  = vec4d_load_aligned(&Q->w2tab[4+j_bits][16*j_mr+0]);
    v  = vec4d_load_aligned(&Q->w2tab[4+j_bits][16*j_mr+4]);
    ww0 = vec4d_unpackhi_permute_3_1_2_0(u, v);
    ww1 = vec4d_unpacklo_permute_3_1_2_0(u, v);
    LENGTH4INVTR_ANY_J(vec4d, x4,x5,x6,x7, n,ninv, w0,ww0,ww1);

    vec4d_store(X+4*0, x0);
    vec4d_store(X+4*1, x1);
    vec4d_store(X+4*2, x2);
    vec4d_store(X+4*3, x3);
    vec4d_store(X+4*4, x4);
    vec4d_store(X+4*5, x5);
    vec4d_store(X+4*6, x6);
    vec4d_store(X+4*7, x7);
}


/* use with n = m-2 and m >= 6 */
#define EXTEND_BASECASE(n, m) \
static void CAT3(sd_ifft_tr_basecase, m, 1)(const sd_fft_ctx_t Q, double* X) \
{ \
    ulong l = n_pow2(m - 2); \
    { \
        _RADIX_4_REVERSE_TR_PARAM_J_IS_Z(VECND, Q) \
        ulong i = 0; do { \
            _RADIX_4_REVERSE_TR_MOTH_J_IS_Z(VECND, X+0*l+i, X+1*l+i, X+2*l+i, X+3*l+i) \
        } while (i += N, i < l); \
        FLINT_ASSERT(i == l); \
    } \
    CAT3(sd_ifft_tr_basecase, n, 1)(Q, X+0*l); \
    CAT3(sd_ifft_tr_basecase, n, 0)(Q, X+1*l, 0, 1); \
    CAT3(sd_ifft_tr_basecase, n, 0)(Q, X+2*l, 1, 2); \
    CAT3(sd_ifft_tr_basecase, n, 0)(Q, X+3*l, 0, 2); \
} \
static void CAT3(sd_ifft_tr_basecase, m, 0)(const sd_fft_ctx_t Q, double* X, ulong j_mr, ulong j_bits) \
{ \
    ulong l = n_pow2(m - 2); \
    FLINT_ASSERT(j_bits != 0); \
    { \
        _RADIX_4_REVERSE_TR_PARAM_J_IS_NZ(VECND, Q, j_mr, j_bits) \
        ulong i = 0; do { \
            _RADIX_4_REVERSE_TR_MOTH_J_IS_NZ(VECND, X+0*l+i, X+1*l+i, X+2*l+i, X+3*l+i) \
        } while (i += N, i < l); \
        FLINT_ASSERT(i == l); \
    } \
    CAT3(sd_ifft_tr_basecase, n, 0)(Q, X+0*l, 4*j_mr+3, 2+j_bits); \
    CAT3(sd_ifft_tr_basecase, n, 0)(Q, X+1*l, 4*j_mr+2, 2+j_bits); \
    CAT3(sd_ifft_tr_basecase, n, 0)(Q, X+2*l, 4*j_mr+1, 2+j_bits); \
    CAT3(sd_ifft_tr_basecase, n, 0)(Q, X+3*l, 4*j_mr+0, 2+j_bits); \
}

EXTEND_BASECASE(4, 6)
EXTEND_BASECASE(5, 7)
EXTEND_BASECASE(6, 8)
EXTEND_BASECASE(7, 9)
#undef EXTEND_BASECASE

/* parameter 1: j can be zero */
static void sd_ifft_tr_base_8_1(const sd_fft_ctx_t Q, double* x, ulong j)
{
    ulong j_bits, j_mr;

    SET_J_BITS_AND_J_MR(j_bits, j_mr, j);

    if (j == 0)
        sd_ifft_tr_basecase_8_1(Q, x);
    else
        sd_ifft_tr_basecase_8_0(Q, x, j_mr, j_bits);
}

/* parameter 0: j cannot be zero */
static void sd_ifft_tr_base_8_0(const sd_fft_ctx_t Q, double* x, ulong j)
{
    ulong j_bits, j_mr;

    FLINT_ASSERT(j != 0);

    SET_J_BITS_AND_J_MR(j_bits, j_mr, j);

    sd_ifft_tr_basecase_8_0(Q, x, j_mr, j_bits);
}

static void sd_ifft_tr_base_9_1(const sd_fft_ctx_t Q, double* x, ulong j)
{
    ulong j_bits, j_mr;

    SET_J_BITS_AND_J_MR(j_bits, j_mr, j);

    if (j == 0)
        sd_ifft_tr_basecase_9_1(Q, x);
    else
        sd_ifft_tr_basecase_9_0(Q, x, j_mr, j_bits);
}


/******** transposed inverse butterfly with truncation, radix 2 ***************/

/* the legal functions must be implemented */
#define DEFINE_IT(nn, zz, ff) \
static void CAT4(radix_2_moth_tr_trunc_block, nn, zz, ff)( \
    const sd_fft_ctx_t FLINT_UNUSED(Q), \
    ulong FLINT_UNUSED(j), \
    double* FLINT_UNUSED(X0), double* FLINT_UNUSED(X1)) \
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


/* transpose of {x0, x1} = {2*x0 - w*x1, x0 - w*x1}, i.e.
   {x0, x1} = {2*x0 + x1, -w*(x0 + x1)} */
static void radix_2_moth_tr_trunc_block_1_2_1(
    const sd_fft_ctx_t Q,
    ulong j,
    double* X0, double* X1)
{
    VECND n    = VECNOP(set_d)(Q->p);
    VECND ninv = VECNOP(set_d)(Q->pinv);
    VECND w = VECNOP(set_d)(sd_fft_ctx_w2(Q, j));
    VECND c = VECNOP(set_d)(2);
    ulong i = 0; do {
        VECND a, b, u, v;
        a = VECNOP(load)(X0 + i);
        /* both operands are inputs in (-3n, 3n): 2*a + b would reach 9n,
           beyond 2^53 for the largest primes, so reduce a first */
        a = VECNOP(reduce_to_pm1n)(a, n, ninv);
        b = VECNOP(load)(X1 + i);
        u = VECNOP(fmadd)(c, a, b);
        v = VECNOP(add)(a, b);
        u = VECNOP(reduce_to_pm1n)(u, n, ninv);
        v = VECNOP(nmulmod_fast)(v, w, n, ninv);
        VECNOP(store)(X0 + i, u);
        VECNOP(store)(X1 + i, v);
    } while (i += N, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/* transpose of {x0} = {2*x0 - w*x1}, i.e. {x0, x1} = {2*x0, -w*x0} */
static void radix_2_moth_tr_trunc_block_1_2_0(
    const sd_fft_ctx_t Q,
    ulong j,
    double* X0, double* X1)
{
    VECND n    = VECNOP(set_d)(Q->p);
    VECND ninv = VECNOP(set_d)(Q->pinv);
    VECND w = VECNOP(set_d)(sd_fft_ctx_w2(Q, j));
    ulong i = 0; do {
        VECND a, u, v;
        a = VECNOP(load)(X0 + i);
        a = VECNOP(reduce_to_pm1n)(a, n, ninv);
        u = VECNOP(add)(a, a);
        v = VECNOP(nmulmod_fast)(a, w, n, ninv);
        VECNOP(store)(X0 + i, u);
        VECNOP(store)(X1 + i, v);
    } while (i += N, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/* transpose of {x0, x1} = {2*x0, x0}, i.e. {x0} = {2*x0 + x1} */
static void radix_2_moth_tr_trunc_block_1_1_1(
    const sd_fft_ctx_t Q,
    ulong FLINT_UNUSED(j),
    double* X0, double* X1)
{
    VECND n    = VECNOP(set_d)(Q->p);
    VECND ninv = VECNOP(set_d)(Q->pinv);
    VECND c = VECNOP(set_d)(2);
    ulong i = 0; do {
        VECND a, b, u;
        a = VECNOP(load)(X0 + i);
        /* as in 1_2_1: 2*a + b with both in (-3n, 3n) would reach 9n */
        a = VECNOP(reduce_to_pm1n)(a, n, ninv);
        b = VECNOP(load)(X1 + i);
        u = VECNOP(fmadd)(c, a, b);
        u = VECNOP(reduce_to_pm1n)(u, n, ninv);
        VECNOP(store)(X0 + i, u);
    } while (i += N, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/* transpose of {x0} = {2*x0}, i.e. {x0} = {2*x0} */
static void radix_2_moth_tr_trunc_block_1_1_0(
    const sd_fft_ctx_t Q,
    ulong FLINT_UNUSED(j),
    double* X0, double* FLINT_UNUSED(X1))
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

/* transpose of {x0} = {(x0 + w*x1)/2}, i.e. {x0, x1} = {x0/2, w*x0/2} */
static void radix_2_moth_tr_trunc_block_0_2_1(
    const sd_fft_ctx_t Q,
    ulong j,
    double* X0, double* X1)
{
    VECND n    = VECNOP(set_d)(Q->p);
    VECND ninv = VECNOP(set_d)(Q->pinv);
    VECND w = VECNOP(set_d)(sd_fft_ctx_w2(Q, j));
    VECND c = VECNOP(set_d)(vec1d_fnmadd(0.5, Q->p, 0.5));
    ulong i = 0; do {
        VECND a, b;
        a = VECNOP(load)(X0 + i);
        a = VECNOP(mulmod_fast)(a, c, n, ninv);
        b = VECNOP(mulmod_fast)(a, w, n, ninv);
        VECNOP(store)(X0 + i, a);
        VECNOP(store)(X1 + i, b);
    } while (i += N, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/* transpose of {x0} = {(x0)/2}, i.e. {x0} = {x0/2} */
static void radix_2_moth_tr_trunc_block_0_1_1(
    const sd_fft_ctx_t Q,
    ulong FLINT_UNUSED(j),
    double* X0, double* FLINT_UNUSED(X1))
{
    VECND n    = VECNOP(set_d)(Q->p);
    VECND ninv = VECNOP(set_d)(Q->pinv);
    VECND c = VECNOP(set_d)(vec1d_fnmadd(0.5, Q->p, 0.5));
    ulong i = 0; do {
        VECND a;
        a = VECNOP(load)(X0 + i);
        a = VECNOP(mulmod_fast)(a, c, n, ninv);
        VECNOP(store)(X0 + i, a);
    } while (i += N, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}


/******** transposed inverse butterfly with truncation, radix 4 ***************/

/* the legal function are opt-in, and illegal/unimplemented are NULL */

#define radix_4_moth_tr_trunc_block_0_1_0 NULL
#define radix_4_moth_tr_trunc_block_0_1_1 NULL
#define radix_4_moth_tr_trunc_block_0_2_0 NULL
#define radix_4_moth_tr_trunc_block_0_2_1 NULL
#define radix_4_moth_tr_trunc_block_0_3_0 NULL
#define radix_4_moth_tr_trunc_block_0_3_1 NULL
#define radix_4_moth_tr_trunc_block_0_4_0 NULL
#define radix_4_moth_tr_trunc_block_1_2_0 NULL
#define radix_4_moth_tr_trunc_block_1_2_1 NULL
#define radix_4_moth_tr_trunc_block_1_3_0 NULL
#define radix_4_moth_tr_trunc_block_1_3_1 NULL
#define radix_4_moth_tr_trunc_block_2_1_0 NULL
#define radix_4_moth_tr_trunc_block_2_1_1 NULL
#define radix_4_moth_tr_trunc_block_2_3_0 NULL
#define radix_4_moth_tr_trunc_block_2_3_1 NULL
#define radix_4_moth_tr_trunc_block_3_1_0 NULL
#define radix_4_moth_tr_trunc_block_3_1_1 NULL
#define radix_4_moth_tr_trunc_block_3_2_0 NULL
#define radix_4_moth_tr_trunc_block_3_2_1 NULL
#define radix_4_moth_tr_trunc_block_4_1_0 NULL
#define radix_4_moth_tr_trunc_block_4_1_1 NULL
#define radix_4_moth_tr_trunc_block_4_2_0 NULL
#define radix_4_moth_tr_trunc_block_4_2_1 NULL
#define radix_4_moth_tr_trunc_block_4_3_0 NULL
#define radix_4_moth_tr_trunc_block_4_3_1 NULL
#define radix_4_moth_tr_trunc_block_4_4_0 NULL
#define radix_4_moth_tr_trunc_block_4_4_1 NULL

/*
Transpose of radix_4_moth_inv_trunc_block_3_4_1, whose matrix is
k = 2, n = 3, z = 4, f = true
[      -r + 1           r + 1         2   r*w^3]
[        2//w           -2//w         0    -w^2]
[(r + 1)//w^2   (-r + 1)//w^2   -2//w^2    -r*w]
[          -r               r         1   r*w^3]

    so that, with u = x0 + x3, s = x2/w^2, t = u - s and A = x0 + s,

    {x0, x1, x2, x3} = {A - (r*t - 2//w*x1),
                        A + (r*t - 2//w*x1),
                        2*(x0 - s) + x3,
                        r*w^3*t - w^2*x1}
*/
static void radix_4_moth_tr_trunc_block_3_4_1(
    const sd_fft_ctx_t Q,
    ulong j, ulong j_bits,
    double* X0, double* X1, double* X2, double* X3)
{
    VECMD n    = VECMOP(set_d)(Q->p);
    VECMD ninv = VECMOP(set_d)(Q->pinv);
    ulong j_mr = n_pow2(j_bits) - 1 - j;
    ulong j_r = j & (n_pow2(j_bits)/2 - 1);
    double W  = FLINT_UNLIKELY(j == 0) ? -1.0 : Q->w2tab[1+j_bits][2*j_mr+1];
    double W2 = FLINT_UNLIKELY(j == 0) ? -1.0 : Q->w2tab[0+j_bits][j_mr];
    double ww2 = Q->w2tab[0+j_bits][j_r];
    double rww = FLINT_UNLIKELY(j == 0) ? Q->w2tab[0][1] : Q->w2tab[1+j_bits][2*j_r+1];
    double rww3 = vec1d_reduce_pm1n_to_pmhn(vec1d_mulmod_fast(rww, ww2, Q->p, Q->pinv), Q->p);
    VECMD f0 = VECMOP(set_d)(Q->w2tab[0][1]);                                     /* r */
    VECMD f1 = VECMOP(set_d)(vec1d_reduce_pm1n_to_pmhn(-2*W, Q->p));              /* 2*w^-1 */
    VECMD f3 = VECMOP(set_d)(-W2);                                                /* w^-2 */
    VECMD fq = VECMOP(set_d)(ww2);                                                /* w^2 */
    VECMD fp = VECMOP(set_d)(rww3);                                               /* r*w^3 */
    ulong i = 0; do {
        VECMD a, b, c, d, s, t, u, v, A, B, C;
        a = VECMOP(load)(X0+i);
        a = VECMOP(reduce_to_pm1n)(a, n, ninv);
        b = VECMOP(load)(X1+i);
        c = VECMOP(load)(X2+i);
        d = VECMOP(load)(X3+i);
        s = VECMOP(mulmod_fast)(c, f3, n, ninv);
        u = VECMOP(reduce_to_pm1n)(VECMOP(add)(a, d), n, ninv);
        t = VECMOP(sub)(u, s);
        A = VECMOP(reduce_to_pm1n)(VECMOP(add)(a, s), n, ninv);
        B = VECMOP(mulmod_fast)(t, f0, n, ninv);
        C = VECMOP(mulmod_fast)(b, f1, n, ninv);
        B = VECMOP(sub)(B, C);
        v = VECMOP(sub)(a, s);
        v = VECMOP(add)(VECMOP(add)(v, v), d);
        v = VECMOP(reduce_to_pm1n)(v, n, ninv);
        VECMOP(store)(X0+i, VECMOP(sub)(A, B));
        VECMOP(store)(X1+i, VECMOP(add)(A, B));
        VECMOP(store)(X2+i, v);
        VECMOP(store)(X3+i, VECMOP(sub)(VECMOP(mulmod_fast)(t, fp, n, ninv),
                                        VECMOP(mulmod_fast)(b, fq, n, ninv)));
    } while (i += M, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/*
Transpose of radix_4_moth_inv_trunc_block_3_4_0, whose matrix is
k = 2, n = 3, z = 4, f = false
[      -r + 1           r + 1         2   r*w^3]
[        2//w           -2//w         0    -w^2]
[(r + 1)//w^2   (-r + 1)//w^2   -2//w^2    -r*w]

    i.e. the previous one with x3 = 0.
*/
static void radix_4_moth_tr_trunc_block_3_4_0(
    const sd_fft_ctx_t Q,
    ulong j, ulong j_bits,
    double* X0, double* X1, double* X2, double* X3)
{
    VECMD n    = VECMOP(set_d)(Q->p);
    VECMD ninv = VECMOP(set_d)(Q->pinv);
    ulong j_mr = n_pow2(j_bits) - 1 - j;
    ulong j_r = j & (n_pow2(j_bits)/2 - 1);
    double W  = FLINT_UNLIKELY(j == 0) ? -1.0 : Q->w2tab[1+j_bits][2*j_mr+1];
    double W2 = FLINT_UNLIKELY(j == 0) ? -1.0 : Q->w2tab[0+j_bits][j_mr];
    double ww2 = Q->w2tab[0+j_bits][j_r];
    double rww = FLINT_UNLIKELY(j == 0) ? Q->w2tab[0][1] : Q->w2tab[1+j_bits][2*j_r+1];
    double rww3 = vec1d_reduce_pm1n_to_pmhn(vec1d_mulmod_fast(rww, ww2, Q->p, Q->pinv), Q->p);
    VECMD f0 = VECMOP(set_d)(Q->w2tab[0][1]);                                     /* r */
    VECMD f1 = VECMOP(set_d)(vec1d_reduce_pm1n_to_pmhn(-2*W, Q->p));              /* 2*w^-1 */
    VECMD f3 = VECMOP(set_d)(-W2);                                                /* w^-2 */
    VECMD fq = VECMOP(set_d)(ww2);                                                /* w^2 */
    VECMD fp = VECMOP(set_d)(rww3);                                               /* r*w^3 */
    ulong i = 0; do {
        VECMD a, b, c, s, t, v, A, B, C;
        a = VECMOP(load)(X0+i);
        b = VECMOP(load)(X1+i);
        c = VECMOP(load)(X2+i);
        s = VECMOP(mulmod_fast)(c, f3, n, ninv);
        a = VECMOP(reduce_to_pm1n)(a, n, ninv);
        t = VECMOP(sub)(a, s);
        A = VECMOP(reduce_to_pm1n)(VECMOP(add)(a, s), n, ninv);
        B = VECMOP(mulmod_fast)(t, f0, n, ninv);
        C = VECMOP(mulmod_fast)(b, f1, n, ninv);
        B = VECMOP(sub)(B, C);
        v = VECMOP(reduce_to_pm1n)(VECMOP(add)(t, t), n, ninv);
        VECMOP(store)(X0+i, VECMOP(sub)(A, B));
        VECMOP(store)(X1+i, VECMOP(add)(A, B));
        VECMOP(store)(X2+i, v);
        VECMOP(store)(X3+i, VECMOP(sub)(VECMOP(mulmod_fast)(t, fp, n, ninv),
                                        VECMOP(mulmod_fast)(b, fq, n, ninv)));
    } while (i += M, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/*
Transpose of radix_4_moth_inv_trunc_block_3_3_1, whose matrix is
k = 2, n = 3, z = 3, f = true
[      -r + 1           r + 1         2]
[        2//w           -2//w         0]
[(r + 1)//w^2   (-r + 1)//w^2   -2//w^2]
[          -r               r         1]

    i.e. the transpose of the 3_4_1 matrix without its last row.
*/
static void radix_4_moth_tr_trunc_block_3_3_1(
    const sd_fft_ctx_t Q,
    ulong j, ulong j_bits,
    double* X0, double* X1, double* X2, double* X3)
{
    VECND n    = VECNOP(set_d)(Q->p);
    VECND ninv = VECNOP(set_d)(Q->pinv);
    ulong j_mr = n_pow2(j_bits) - 1 - j;
    double W  = FLINT_UNLIKELY(j == 0) ? -1.0 : Q->w2tab[1+j_bits][2*j_mr+1];
    double W2 = FLINT_UNLIKELY(j == 0) ? -1.0 : Q->w2tab[0+j_bits][j_mr];
    VECND f0 = VECNOP(set_d)(Q->w2tab[0][1]);                                     /* r */
    VECND f1 = VECNOP(set_d)(vec1d_reduce_pm1n_to_pmhn(-2*W, Q->p));              /* 2*w^-1 */
    VECND f3 = VECNOP(set_d)(-W2);                                                /* w^-2 */
    ulong i = 0; do {
        VECND a, b, c, d, s, t, u, v, A, B, C;
        a = VECNOP(load)(X0+i);
        a = VECNOP(reduce_to_pm1n)(a, n, ninv);
        b = VECNOP(load)(X1+i);
        c = VECNOP(load)(X2+i);
        d = VECNOP(load)(X3+i);
        s = VECNOP(mulmod_fast)(c, f3, n, ninv);
        u = VECNOP(reduce_to_pm1n)(VECNOP(add)(a, d), n, ninv);
        t = VECNOP(sub)(u, s);
        A = VECNOP(reduce_to_pm1n)(VECNOP(add)(a, s), n, ninv);
        B = VECNOP(mulmod_fast)(t, f0, n, ninv);
        C = VECNOP(mulmod_fast)(b, f1, n, ninv);
        B = VECNOP(sub)(B, C);
        v = VECNOP(sub)(a, s);
        v = VECNOP(add)(VECNOP(add)(v, v), d);
        v = VECNOP(reduce_to_pm1n)(v, n, ninv);
        VECNOP(store)(X0+i, VECNOP(sub)(A, B));
        VECNOP(store)(X1+i, VECNOP(add)(A, B));
        VECNOP(store)(X2+i, v);
    } while (i += N, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/*
Transpose of radix_4_moth_inv_trunc_block_3_3_0, whose matrix is
k = 2, n = 3, z = 3, f = false
[      -r + 1           r + 1         2]
[        2//w           -2//w         0]
[(r + 1)//w^2   (-r + 1)//w^2   -2//w^2]
*/
static void radix_4_moth_tr_trunc_block_3_3_0(
    const sd_fft_ctx_t Q,
    ulong j, ulong j_bits,
    double* X0, double* X1, double* X2, double* FLINT_UNUSED(X3))
{
    VECND n    = VECNOP(set_d)(Q->p);
    VECND ninv = VECNOP(set_d)(Q->pinv);
    ulong j_mr = n_pow2(j_bits) - 1 - j;
    double W  = FLINT_UNLIKELY(j == 0) ? -1.0 : Q->w2tab[1+j_bits][2*j_mr+1];
    double W2 = FLINT_UNLIKELY(j == 0) ? -1.0 : Q->w2tab[0+j_bits][j_mr];
    VECND f0 = VECNOP(set_d)(Q->w2tab[0][1]);                                     /* r */
    VECND f1 = VECNOP(set_d)(vec1d_reduce_pm1n_to_pmhn(-2*W, Q->p));              /* 2*w^-1 */
    VECND f3 = VECNOP(set_d)(-W2);                                                /* w^-2 */
    ulong i = 0; do {
        VECND a, b, c, s, t, A, B, C;
        a = VECNOP(load)(X0+i);
        b = VECNOP(load)(X1+i);
        c = VECNOP(load)(X2+i);
        s = VECNOP(mulmod_fast)(c, f3, n, ninv);
        a = VECNOP(reduce_to_pm1n)(a, n, ninv);
        t = VECNOP(sub)(a, s);
        A = VECNOP(reduce_to_pm1n)(VECNOP(add)(a, s), n, ninv);
        B = VECNOP(mulmod_fast)(t, f0, n, ninv);
        C = VECNOP(mulmod_fast)(b, f1, n, ninv);
        B = VECNOP(sub)(B, C);
        VECNOP(store)(X0+i, VECNOP(sub)(A, B));
        VECNOP(store)(X1+i, VECNOP(add)(A, B));
        VECNOP(store)(X2+i, VECNOP(reduce_to_pm1n)(VECNOP(add)(t, t), n, ninv));
    } while (i += N, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/*
Transpose of radix_4_moth_inv_trunc_block_2_4_1, whose matrix is
k = 2, n = 2, z = 4, f = true
[            2                2        -w^2             0]
[         2//w            -2//w           0          -w^2]
[1//2*r + 1//2   -1//2*r + 1//2   -1//2*w^2   -1//2*r*w^3]

    so that, with h = x2/2,

    {x0, x1, x2, x3} = {2*x0 + h + (2//w*x1 + r*h),
                        2*x0 + h - (2//w*x1 + r*h),
                        -w^2*(x0 + h),
                        -w^2*x1 - r*w^3*h}
*/
static void radix_4_moth_tr_trunc_block_2_4_1(
    const sd_fft_ctx_t Q,
    ulong j, ulong j_bits,
    double* X0, double* X1, double* X2, double* X3)
{
    VECMD n    = VECMOP(set_d)(Q->p);
    VECMD ninv = VECMOP(set_d)(Q->pinv);
    ulong j_mr = n_pow2(j_bits) - 1 - j;
    ulong j_r = j & (n_pow2(j_bits)/2 - 1);
    double W  = FLINT_UNLIKELY(j == 0) ? -1.0 : Q->w2tab[1+j_bits][2*j_mr+1];
    double ww2 = Q->w2tab[0+j_bits][j_r];
    double rww = FLINT_UNLIKELY(j == 0) ? Q->w2tab[0][1] : Q->w2tab[1+j_bits][2*j_r+1];
    double rww3 = vec1d_reduce_pm1n_to_pmhn(vec1d_mulmod_fast(rww, ww2, Q->p, Q->pinv), Q->p);
    VECMD f1 = VECMOP(set_d)(vec1d_reduce_pm1n_to_pmhn(-2*W, Q->p));              /* 2*w^-1 */
    VECMD f2 = VECMOP(set_d)(vec1d_fnmadd(0.5, Q->p, 0.5));                       /* 1/2 */
    VECMD f3 = VECMOP(set_d)(Q->w2tab[0][1]);                                     /* r */
    VECMD f4 = VECMOP(set_d)(ww2);                                                /* w^2 */
    VECMD f5 = VECMOP(set_d)(rww3);                                               /* r*w^3 */
    ulong i = 0; do {
        VECMD a, b, c, h, p, R, C;
        a = VECMOP(load)(X0+i);
        b = VECMOP(load)(X1+i);
        c = VECMOP(load)(X2+i);
        a = VECMOP(reduce_to_pm1n)(a, n, ninv);
        h = VECMOP(mulmod_fast)(c, f2, n, ninv);
        R = VECMOP(mulmod_fast)(h, f3, n, ninv);
        C = VECMOP(mulmod_fast)(b, f1, n, ninv);
        C = VECMOP(add)(C, R);
        p = VECMOP(add)(VECMOP(add)(a, a), h);
        p = VECMOP(reduce_to_pm1n)(p, n, ninv);
        VECMOP(store)(X0+i, VECMOP(add)(p, C));
        VECMOP(store)(X1+i, VECMOP(sub)(p, C));
        VECMOP(store)(X2+i, VECMOP(nmulmod_fast)(VECMOP(add)(a, h), f4, n, ninv));
        VECMOP(store)(X3+i, VECMOP(sub)(VECMOP(nmulmod_fast)(b, f4, n, ninv),
                                        VECMOP(mulmod_fast)(h, f5, n, ninv)));
    } while (i += M, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/*
Transpose of radix_4_moth_inv_trunc_block_2_4_0, whose matrix is
k = 2, n = 2, z = 4, f = false
[   2       2   -w^2      0]
[2//w   -2//w      0   -w^2]

    {x0, x1, x2, x3} = {2*x0 + 2//w*x1, 2*x0 - 2//w*x1, -w^2*x0, -w^2*x1}
*/
static void radix_4_moth_tr_trunc_block_2_4_0(
    const sd_fft_ctx_t Q,
    ulong j, ulong j_bits,
    double* X0, double* X1, double* X2, double* X3)
{
    VECND n    = VECNOP(set_d)(Q->p);
    VECND ninv = VECNOP(set_d)(Q->pinv);
    ulong j_mr = n_pow2(j_bits) - 1 - j;
    ulong j_r = j & (n_pow2(j_bits)/2 - 1);
    double W  = FLINT_UNLIKELY(j == 0) ? -1.0 : Q->w2tab[1+j_bits][2*j_mr+1];
    double ww2 = Q->w2tab[0+j_bits][j_r];
    VECND w2 = VECNOP(set_d)(ww2);                                                /* w^2 */
    VECND twowi = VECNOP(set_d)(vec1d_reduce_pm1n_to_pmhn(-2*W, Q->p));           /* 2*w^-1 */
    ulong i = 0; do {
        VECND a, b, u, v;
        a = VECNOP(load)(X0+i);
        b = VECNOP(load)(X1+i);
        u = VECNOP(reduce_to_pm1n)(VECNOP(add)(a, a), n, ninv);
        v = VECNOP(mulmod_fast)(b, twowi, n, ninv);
        VECNOP(store)(X2+i, VECNOP(nmulmod_fast)(a, w2, n, ninv));
        VECNOP(store)(X3+i, VECNOP(nmulmod_fast)(b, w2, n, ninv));
        VECNOP(store)(X0+i, VECNOP(add)(u, v));
        VECNOP(store)(X1+i, VECNOP(sub)(u, v));
    } while (i += N, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/*
Transpose of radix_4_moth_inv_trunc_block_2_2_1, whose matrix is
k = 2, n = 2, z = 2, f = true
[            2                2]
[         2//w            -2//w]
[1//2*r + 1//2   -1//2*r + 1//2]

    so that, with h = x2/2,
    {x0, x1} = {2*x0 + h + (2//w*x1 + r*h), 2*x0 + h - (2//w*x1 + r*h)}
*/
static void radix_4_moth_tr_trunc_block_2_2_1(
    const sd_fft_ctx_t Q,
    ulong j, ulong j_bits,
    double* X0, double* X1, double* X2, double* FLINT_UNUSED(X3))
{
    VECND n    = VECNOP(set_d)(Q->p);
    VECND ninv = VECNOP(set_d)(Q->pinv);
    ulong j_mr = n_pow2(j_bits) - 1 - j;
    double W  = FLINT_UNLIKELY(j == 0) ? -1.0 : Q->w2tab[1+j_bits][2*j_mr+1];
    VECND c1 = VECNOP(set_d)(vec1d_reduce_pm1n_to_pmhn(-2*W, Q->p));              /* 2/w */
    VECND c2 = VECNOP(set_d)(vec1d_fnmadd(0.5, Q->p, 0.5));                       /* 1/2 */
    VECND c3 = VECNOP(set_d)(Q->w2tab[1][0]);                                     /* r */
    ulong i = 0; do {
        VECND a, b, c, h, p, s;
        a = VECNOP(load)(X0+i);
        b = VECNOP(load)(X1+i);
        c = VECNOP(load)(X2+i);
        a = VECNOP(reduce_to_pm1n)(a, n, ninv);
        h = VECNOP(mulmod_fast)(c, c2, n, ninv);
        s = VECNOP(mulmod_fast)(h, c3, n, ninv);
        s = VECNOP(add)(s, VECNOP(mulmod_fast)(b, c1, n, ninv));
        p = VECNOP(add)(VECNOP(add)(a, a), h);
        p = VECNOP(reduce_to_pm1n)(p, n, ninv);
        VECNOP(store)(X0+i, VECNOP(add)(p, s));
        VECNOP(store)(X1+i, VECNOP(sub)(p, s));
    } while (i += N, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/*
Transpose of radix_4_moth_inv_trunc_block_2_2_0
k = 2, n = 2, z = 2, f = 0

{x0, x1} = {2*x0 + 2//w*x1, 2*x0 - 2//w*x1}
*/
static void radix_4_moth_tr_trunc_block_2_2_0(
    const sd_fft_ctx_t Q,
    ulong j, ulong j_bits,
    double* X0, double* X1, double* FLINT_UNUSED(X2), double* FLINT_UNUSED(X3))
{
    VECND n    = VECNOP(set_d)(Q->p);
    VECND ninv = VECNOP(set_d)(Q->pinv);
    ulong j_mr = n_pow2(j_bits) - 1 - j;
    double W  = FLINT_UNLIKELY(j == 0) ? -1.0 : Q->w2tab[1+j_bits][2*j_mr+1];
    VECND c1 = VECNOP(set_d)(vec1d_reduce_pm1n_to_pmhn(-2*W, Q->p));              /* 2/w */
    ulong i = 0; do {
        VECND a, b, u, v;
        a = VECNOP(load)(X0+i);
        b = VECNOP(load)(X1+i);
        u = VECNOP(reduce_to_pm1n)(VECNOP(add)(a, a), n, ninv);
        v = VECNOP(mulmod_fast)(b, c1, n, ninv);
        VECNOP(store)(X0+i, VECNOP(add)(u, v));
        VECNOP(store)(X1+i, VECNOP(sub)(u, v));
    } while (i += N, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/*
Transpose of radix_4_moth_inv_trunc_block_1_4_1, whose matrix is
k = 2, n = 1, z = 4, f = true
[4        -w   -w^2        -w^3]
[1   -1//2*w      0   -1//2*w^3]

    so that, with g = x0 + x1/2,
    {x0, x1, x2, x3} = {4*x0 + x1, -w*g, -w^2*x0, w^2*(-w*g)}
*/
static void radix_4_moth_tr_trunc_block_1_4_1(
    const sd_fft_ctx_t Q,
    ulong j, ulong j_bits,
    double* X0, double* X1, double* X2, double* X3)
{
    VECMD n    = VECMOP(set_d)(Q->p);
    VECMD ninv = VECMOP(set_d)(Q->pinv);
    ulong j_r = j & (n_pow2(j_bits)/2 - 1);
    double ww2 = Q->w2tab[0+j_bits][j_r];
    double ww  = FLINT_UNLIKELY(j == 0) ? 1.0 : Q->w2tab[1+j_bits][2*j_r];
    VECMD f4 = VECMOP(set_d)(4);
    VECMD f2 = VECMOP(set_d)(vec1d_fnmadd(0.5, Q->p, 0.5));                       /* 1/2 */
    VECMD w2 = VECMOP(set_d)(ww2);                                                /* w^2 */
    VECMD w  = VECMOP(set_d)(ww);                                                 /* w */
    ulong i = 0; do {
        VECMD a, b, g, u, v;
        a = VECMOP(load)(X0+i);
        a = VECMOP(reduce_to_pm1n)(a, n, ninv);
        b = VECMOP(load)(X1+i);
        u = VECMOP(reduce_to_pm1n)(VECMOP(fmadd)(f4, a, b), n, ninv);
        g = VECMOP(add)(a, VECMOP(mulmod_fast)(b, f2, n, ninv));
        g = VECMOP(reduce_to_pm1n)(g, n, ninv);
        v = VECMOP(nmulmod_fast)(g, w, n, ninv);
        VECMOP(store)(X2+i, VECMOP(nmulmod_fast)(a, w2, n, ninv));
        VECMOP(store)(X3+i, VECMOP(mulmod_fast)(v, w2, n, ninv));
        VECMOP(store)(X0+i, u);
        VECMOP(store)(X1+i, v);
    } while (i += M, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/*
Transpose of radix_4_moth_inv_trunc_block_1_4_0, whose matrix is
k = 2, n = 1, z = 4, f = false
[4   -w   -w^2   -w^3]

    {x0, x1, x2, x3} = {4*x0, -w*x0, -w^2*x0, w^2*(-w*x0)}
*/
static void radix_4_moth_tr_trunc_block_1_4_0(
    const sd_fft_ctx_t Q,
    ulong j, ulong j_bits,
    double* X0, double* X1, double* X2, double* X3)
{
    VECND n    = VECNOP(set_d)(Q->p);
    VECND ninv = VECNOP(set_d)(Q->pinv);
    ulong j_r = j & (n_pow2(j_bits)/2 - 1);
    double ww2 = Q->w2tab[0+j_bits][j_r];
    double ww  = FLINT_UNLIKELY(j == 0) ? 1.0 : Q->w2tab[1+j_bits][2*j_r];
    VECND f1 = VECNOP(set_d)(4.0);
    VECND w2 = VECNOP(set_d)(ww2);
    VECND w  = VECNOP(set_d)(ww);
    ulong i = 0; do {
        VECND a, b;
        a = VECNOP(load)(X0+i);
        a = VECNOP(reduce_to_pm1n)(a, n, ninv);
        b = VECNOP(nmulmod_fast)(a, w, n, ninv);
        VECNOP(store)(X3+i, VECNOP(mulmod_fast)(b, w2, n, ninv));
        VECNOP(store)(X2+i, VECNOP(nmulmod_fast)(a, w2, n, ninv));
        VECNOP(store)(X1+i, b);
        VECNOP(store)(X0+i, VECNOP(reduce_to_pm1n)(VECNOP(mul)(a, f1), n, ninv));
    } while (i += N, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/*
Transpose of radix_4_moth_inv_trunc_block_1_1_1, whose matrix is
k = 2, n = 1, z = 1, f = true
[4]
[1]

    {x0} = {4*x0 + x1}
*/
static void radix_4_moth_tr_trunc_block_1_1_1(
    const sd_fft_ctx_t Q,
    ulong FLINT_UNUSED(j), ulong FLINT_UNUSED(j_bits),
    double* X0, double* X1, double* FLINT_UNUSED(X2), double* FLINT_UNUSED(X3))
{
    VECND n    = VECNOP(set_d)(Q->p);
    VECND ninv = VECNOP(set_d)(Q->pinv);
    VECND f = VECNOP(set_d)(4.0);
    ulong i = 0; do {
        VECND a, b;
        a = VECNOP(load)(X0+i);
        a = VECNOP(reduce_to_pm1n)(a, n, ninv);
        b = VECNOP(load)(X1+i);
        a = VECNOP(fmadd)(f, a, b);
        VECNOP(store)(X0+i, VECNOP(reduce_to_pm1n)(a, n, ninv));
    } while (i += N, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/*
Transpose of radix_4_moth_inv_trunc_block_1_1_0, whose matrix is
k = 2, n = 1, z = 1, f = false
[4]
*/
static void radix_4_moth_tr_trunc_block_1_1_0(
    const sd_fft_ctx_t Q,
    ulong FLINT_UNUSED(j), ulong FLINT_UNUSED(j_bits),
    double* X0, double* FLINT_UNUSED(X1), double* FLINT_UNUSED(X2), double* FLINT_UNUSED(X3))
{
    VECND n    = VECNOP(set_d)(Q->p);
    VECND ninv = VECNOP(set_d)(Q->pinv);
    VECND f = VECNOP(set_d)(4.0);
    ulong i = 0; do {
        VECND a, b;
        a = VECNOP(load)(X0+i);
        a = VECNOP(reduce_to_pm1n)(a, n, ninv);
        b = VECNOP(mul)(f, a);
        VECNOP(store)(X0+i, VECNOP(reduce_to_pm1n)(b, n, ninv));
    } while (i += N, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/*
Transpose of radix_4_moth_inv_trunc_block_0_4_1, whose matrix is
k = 2, n = 0, z = 4, f = true
[1//4   1//4*w   1//4*w^2   1//4*w^3]

    so that, with q = x0/4,
    {x0, x1, x2, x3} = {q, w*q, w^2*q, w^2*(w*q)}
*/
static void radix_4_moth_tr_trunc_block_0_4_1(
    const sd_fft_ctx_t Q,
    ulong j, ulong j_bits,
    double* X0, double* X1, double* X2, double* X3)
{
    VECND n    = VECNOP(set_d)(Q->p);
    VECND ninv = VECNOP(set_d)(Q->pinv);
    ulong j_r = j & (n_pow2(j_bits)/2 - 1);
    double ww2 = Q->w2tab[0+j_bits][j_r];
    double ww  = FLINT_UNLIKELY(j == 0) ? 1.0 : Q->w2tab[1+j_bits][2*j_r];
    VECND one4th = VECNOP(set_d)(vec1d_fnmadd(0.25, Q->p, 0.25));
    VECND w2 = VECNOP(set_d)(ww2);
    VECND w  = VECNOP(set_d)(ww);
    ulong i = 0; do {
        VECND a, b;
        a = VECNOP(load)(X0+i);
        a = VECNOP(mulmod_fast)(a, one4th, n, ninv);
        b = VECNOP(mulmod_fast)(a, w, n, ninv);
        VECNOP(store)(X0+i, a);
        VECNOP(store)(X1+i, b);
        VECNOP(store)(X2+i, VECNOP(mulmod_fast)(a, w2, n, ninv));
        VECNOP(store)(X3+i, VECNOP(mulmod_fast)(b, w2, n, ninv));
    } while (i += N, i < BLK_SZ);
    FLINT_ASSERT(i == BLK_SZ);
}

/************************ the recursive stuff ********************************/

static void sd_ifft_tr_no_trunc_block(
    const sd_fft_ctx_t Q,
    double* x,
    ulong S, /* stride */
    ulong k, /* BLK_SZ transforms each of length 2^k */
    ulong j)
{
    ulong j_bits, j_mr;

    if (k > 2)
    {
        ulong k1 = k/2;
        ulong k2 = k - k1;

        /* column ffts */
        ulong l2 = n_pow2(k2);
        ulong a = 0; do {
            sd_ifft_tr_no_trunc_block(Q, x + BLK_SZ*(a*S), S<<k2, k1, j);
        } while (a++, a < l2);

        /* row ffts */
        ulong l1 = n_pow2(k1);
        ulong b = 0; do {
            sd_ifft_tr_no_trunc_block(Q, x + BLK_SZ*((b<<k2)*S), S, k2, (j<<k1) + b);
        } while (b++, b < l1);

        return;
    }

    SET_J_BITS_AND_J_MR(j_bits, j_mr, j);

    if (k == 2)
    {
        double* X0 = x + BLK_SZ*(S*0);
        double* X1 = x + BLK_SZ*(S*1);
        double* X2 = x + BLK_SZ*(S*2);
        double* X3 = x + BLK_SZ*(S*3);
        if (FLINT_UNLIKELY(j_bits == 0))
        {
            RADIX_4_REVERSE_TR_PARAM_J_IS_Z(VECND, Q)
            ulong i = 0; do {
                RADIX_4_REVERSE_TR_MOTH_J_IS_Z(VECND, X0+i, X1+i, X2+i, X3+i);
            } while(i += N, i < BLK_SZ);
        }
        else
        {
            RADIX_4_REVERSE_TR_PARAM_J_IS_NZ(VECND, Q, j_mr, j_bits)
            ulong i = 0; do {
                RADIX_4_REVERSE_TR_MOTH_J_IS_NZ(VECND, X0+i, X1+i, X2+i, X3+i);
            } while(i += N, i < BLK_SZ);
        }
    }
    else if (k == 1)
    {
        double* X0 = x + BLK_SZ*(S*0);
        double* X1 = x + BLK_SZ*(S*1);
        if (FLINT_UNLIKELY(j_bits == 0))
        {
            RADIX_2_REVERSE_TR_PARAM_J_IS_Z(VECND, Q)
            ulong i = 0; do {
                RADIX_2_REVERSE_TR_MOTH_J_IS_Z(VECND, X0+i, X1+i);
            } while (i += N, i < BLK_SZ);
        }
        else
        {
            RADIX_2_REVERSE_TR_PARAM_J_IS_NZ(VECND, Q, j_mr, j_bits)
            ulong i = 0; do {
                RADIX_2_REVERSE_TR_MOTH_J_IS_NZ(VECND, X0+i, X1+i);
            } while (i += N, i < BLK_SZ);
        }
    }
}

static void sd_ifft_tr_no_trunc_internal(
    const sd_fft_ctx_t Q,
    double* x,
    ulong k, /* transform length BLK_SZ*2^k */
    ulong j)
{
    if (k > 2)
    {
        ulong k1 = k/2;
        ulong k2 = k - k1;

        ulong l2 = n_pow2(k2);
        ulong a = 0; do {
            sd_ifft_tr_no_trunc_block(Q, x + BLK_SZ*a, n_pow2(k2), k1, j);
        } while (a++, a < l2);

        ulong l1 = n_pow2(k1);
        ulong b = 0; do {
            sd_ifft_tr_no_trunc_internal(Q, x + BLK_SZ*(b<<k2), k2, (j<<k1) + b);
        } while (b++, b < l1);

        return;
    }

    if (k == 2)
    {
        /* k1 = 2; k2 = 0 */
        sd_ifft_tr_no_trunc_block(Q, x, 1, 2, j);
        sd_ifft_tr_base_8_1(Q, x + BLK_SZ*0, 4*j+0);
        sd_ifft_tr_base_8_0(Q, x + BLK_SZ*1, 4*j+1);
        sd_ifft_tr_base_8_0(Q, x + BLK_SZ*2, 4*j+2);
        sd_ifft_tr_base_8_0(Q, x + BLK_SZ*3, 4*j+3);
    }
    else if (k == 1)
    {
        sd_ifft_tr_base_9_1(Q, x, j);
    }
    else
    {
        sd_ifft_tr_base_8_1(Q, x, j);
    }
}

static void sd_ifft_tr_trunc_block(
    const sd_fft_ctx_t Q,
    double* x,  /* data + BLK_SZ*I */
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
        sd_ifft_tr_no_trunc_block(Q, x, S, k, j);
        return;
    }

    if (k == 2)
    {
#define IT(nn, zz, ff) CAT4(radix_4_moth_tr_trunc_block, nn, zz, ff)
#define LOOKUP_IT(nn, zz, ff) tab[(ulong)(ff) + 2*((zz)-1 + 4*(nn))]
        static void (*tab[5*4*2])(const sd_fft_ctx_t, ulong, ulong, double*, double*, double*, double*) =
            {IT(0,1,0),IT(0,1,1), IT(0,2,0),IT(0,2,1), IT(0,3,0),IT(0,3,1), IT(0,4,0),IT(0,4,1),
             IT(1,1,0),IT(1,1,1), IT(1,2,0),IT(1,2,1), IT(1,3,0),IT(1,3,1), IT(1,4,0),IT(1,4,1),
             IT(2,1,0),IT(2,1,1), IT(2,2,0),IT(2,2,1), IT(2,3,0),IT(2,3,1), IT(2,4,0),IT(2,4,1),
             IT(3,1,0),IT(3,1,1), IT(3,2,0),IT(3,2,1), IT(3,3,0),IT(3,3,1), IT(3,4,0),IT(3,4,1),
             IT(4,1,0),IT(4,1,1), IT(4,2,0),IT(4,2,1), IT(4,3,0),IT(4,3,1), IT(4,4,0),IT(4,4,1)};

        static void (*fxn)(const sd_fft_ctx_t, ulong, ulong, double*, double*, double*, double*);
        fxn = LOOKUP_IT(n,z,f);
        if (FLINT_LIKELY(fxn != NULL))
        {
            double* X0 = x + BLK_SZ*(S*0);
            double* X1 = x + BLK_SZ*(S*1);
            double* X2 = x + BLK_SZ*(S*2);
            double* X3 = x + BLK_SZ*(S*3);
            fxn(Q, j, n_nbits(j), X0, X1, X2, X3);
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

        /* leftmost columns */
        for (ulong a = 0; a < n2; a++)
            sd_ifft_tr_trunc_block(Q, x + BLK_SZ*(a*S), S << k2, k1, j, z1 + (a < m), n1 + 1, 0);

        /* last partial row */
        if (fp)
            sd_ifft_tr_trunc_block(Q, x + BLK_SZ*(n1*(S << k2)), S, k2, (j << k1) + n1, z2p, n2, f);

        /* rightmost columns */
        for (ulong a = n2; a < z2p; a++)
            sd_ifft_tr_trunc_block(Q, x + BLK_SZ*(a*S), S << k2, k1, j, z1 + (a < mp), n1, fp);

        /* complete rows */
        for (ulong b = 0; b < n1; b++)
            sd_ifft_tr_no_trunc_block(Q, x + BLK_SZ*(b*(S << k2)), S, k2, (j << k1) + b);

        return;
    }

    if (k == 1)
    {
#define IT(nn, zz, ff) CAT4(radix_2_moth_tr_trunc_block, nn, zz, ff)
#define LOOKUP_IT(nn, zz, ff) tab[(ulong)(ff) + 2*((zz)-1 + 2*(nn))]
        static void (*tab[3*2*2*2])(const sd_fft_ctx_t, ulong, double*, double*) =
            {IT(0,1,0),IT(0,1,1), IT(0,2,0),IT(0,2,1),
             IT(1,1,0),IT(1,1,1), IT(1,2,0),IT(1,2,1),
             IT(2,1,0),IT(2,1,1), IT(2,2,0),IT(2,2,1)};

        LOOKUP_IT(n, z, f)(Q, j, x + BLK_SZ*(S*0), x + BLK_SZ*(S*1));
        return;
#undef LOOKUP_IT
#undef IT
    }
}


static void sd_ifft_tr_trunc_internal(
    const sd_fft_ctx_t Q,
    double* x,  /* x = data + BLK_SZ*I  where I = starting index */
    ulong k,    /* transform length 2^(k + LG_BLK_SZ) */
    ulong j,
    ulong z,    /* actual trunc is z*BLK_SZ */
    ulong n,    /* actual trunc is n*BLK_SZ */
    int f)
{
    FLINT_ASSERT(f == 0 || f == 1);
    FLINT_ASSERT(n <= z);
    FLINT_ASSERT(1 <= BLK_SZ*z && BLK_SZ*z <= n_pow2(k+LG_BLK_SZ));
    FLINT_ASSERT(1 <= BLK_SZ*n+f && BLK_SZ*n+f <= n_pow2(k+LG_BLK_SZ));

    if (!f && z == n && n == n_pow2(k))
    {
        sd_ifft_tr_no_trunc_internal(Q, x, k, j);
        return;
    }

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

        /* leftmost columns */
        for (ulong a = 0; a < n2; a++)
            sd_ifft_tr_trunc_block(Q, x + BLK_SZ*a, n_pow2(k2), k1, j, z1 + (a < m), n1 + 1, 0);

        /* last partial row */
        if (fp)
            sd_ifft_tr_trunc_internal(Q, x + BLK_SZ*(n1 << k2), k2, (j << k1) + n1, z2p, n2, f);

        /* rightmost columns */
        for (ulong a = n2; a < z2p; a++)
            sd_ifft_tr_trunc_block(Q, x + BLK_SZ*a, n_pow2(k2), k1, j, z1 + (a < mp), n1, fp);

        /* complete rows */
        for (ulong b = 0; b < n1; b++)
            sd_ifft_tr_no_trunc_internal(Q, x + BLK_SZ*(b << k2), k2, (j << k1) + b);

        return;
    }

    if (k == 2)
    {
        if (f) sd_ifft_tr_trunc_internal(Q, x + BLK_SZ*n, 0, 4*j+n, 1, 0, f);
        sd_ifft_tr_trunc_block(Q, x, 1, 2, j, z, n, f);
                   sd_ifft_tr_base_8_1(Q, x + BLK_SZ*0, 4*j+0);
        if (n > 1) sd_ifft_tr_base_8_0(Q, x + BLK_SZ*1, 4*j+1);
        if (n > 2) sd_ifft_tr_base_8_0(Q, x + BLK_SZ*2, 4*j+2);
        if (n > 3) sd_ifft_tr_base_8_0(Q, x + BLK_SZ*3, 4*j+3);
    }
    else if (k == 1)
    {
        if (f) sd_ifft_tr_trunc_internal(Q, x + BLK_SZ*n, 0, 2*j+n, 1, 0, f);
        sd_ifft_tr_trunc_block(Q, x, 1, 1, j, z, n, f);
                   sd_ifft_tr_base_8_1(Q, x + BLK_SZ*0, 2*j+0);
        if (n > 1) sd_ifft_tr_base_8_0(Q, x + BLK_SZ*1, 2*j+1);
    }
    else
    {
        sd_ifft_tr_base_8_1(Q, x, j);
    }
}


/********************* interface functions ***********************/

/*
Transpose of the truncated inverse Fourier transform `sd_ifft_trunc`, in place
in `d`, with the same parameter, so that transposing a call means adding the
suffix. Both transforms are linear maps on the first
`BLK_SZ*ceil(trunc/BLK_SZ)` entries of `d` (on the first `2^L` entries when
`L <= LG_BLK_SZ`), and
    <sd_ifft_trunc(a), v> == <a, sd_ifft_trunc_tr(v)>.
The array `d` needs to have size at least `n_pow2(L)`.
*/
void sd_ifft_trunc_tr(
    sd_fft_ctx_t Q,
    double* d,
    ulong L,    /* convolution length 2^L */
    ulong trunc)
{
    FLINT_ASSERT(trunc <= n_pow2(L));

    sd_fft_ctx_fit_depth(Q, L);

    if (L > LG_BLK_SZ)
    {
        ulong new_trunc = n_cdiv(trunc, BLK_SZ);

        sd_ifft_tr_trunc_internal(Q, d, L - LG_BLK_SZ, 0, new_trunc, new_trunc, 0);
        return;
    }

    switch (L) {
        case 0: sd_ifft_tr_basecase_0_1(Q, d); break;
        case 1: sd_ifft_tr_basecase_1_1(Q, d); break;
        case 2: sd_ifft_tr_basecase_2_1(Q, d); break;
        case 3: sd_ifft_tr_basecase_3_1(Q, d); break;
        case 4: sd_ifft_tr_basecase_4_1(Q, d); break;
        case 5: sd_ifft_tr_basecase_5_1(Q, d); break;
        case 6: sd_ifft_tr_basecase_6_1(Q, d); break;
        case 7: sd_ifft_tr_basecase_7_1(Q, d); break;
        case 8: sd_ifft_tr_basecase_8_1(Q, d); break;
        default: FLINT_ASSERT(0);
    }
}
