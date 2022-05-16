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

/*
inverse butterfly:

    2*a0 =      (b0 + b1)
    2*a1 = w^-1*(b0 - b1)

    W  := -w^-1
*/
#define RADIX_2_REVERSE_PARAM(V, Q, j) \
    V W = V##_set_d((UNLIKELY((j) == 0)) ? -Q->w2s[0] : Q->w2s[(j)^(n_saturate_bits(j)>>1)]); \
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


/*
inverse butterfly:

    4*a0 =            (b0 + b1) +        (b2 + b3)
    4*a1 =       w^-1*(b0 - b1) - i*w^-1*(b2 - b3)
    4*a2 = w^-2*(     (b0 + b1) -        (b2 + b3))
    4*a3 = w^-2*(w^-1*(b0 - b1) + i*w^-1*(b2 - b3))

    W  := -w^-1
    W2 := -w^-2
    IW := i*w^-1
*/
#define RADIX_4_REVERSE_PARAM(V, Q, j, jm, j_can_be_0) \
    FLINT_ASSERT((j) == 0 || ((jm) == (j)^(n_saturate_bits(j)>>1))); \
    V W  = V##_set_d(((j_can_be_0) && UNLIKELY((j) == 0)) ? -Q->w2s[0] : Q->w2s[2*(jm)+1]); \
    V W2 = V##_set_d(((j_can_be_0) && UNLIKELY((j) == 0)) ? -Q->w2s[0] : Q->w2s[(jm)]); \
    V IW = V##_set_d(((j_can_be_0) && UNLIKELY((j) == 0)) ?  Q->w2s[1] : Q->w2s[2*(jm)]); \
    V n    = V##_set_d(Q->p); \
    V ninv = V##_set_d(Q->pinv); \

#define RADIX_4_REVERSE_PARAM_SURE(V, Q, j, jm, j_is_0) \
    FLINT_ASSERT((j) == 0 || ((jm) == (j)^(n_saturate_bits(j)>>1))); \
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




#define DEFINE_IT(nn, zz, ff) \
static void CAT4(radix_2_moth_inv_trunc_block, nn, zz, ff)( \
    const sd_fft_ctx_t Q, \
    double* X0, double* X1, \
    ulong j) \
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
DEFINE_IT(0,1,1)
DEFINE_IT(0,2,0)
DEFINE_IT(0,2,1)
DEFINE_IT(1,1,0)
DEFINE_IT(1,1,1)
DEFINE_IT(1,2,0)
DEFINE_IT(1,2,1)
DEFINE_IT(2,1,0)
DEFINE_IT(2,1,1)
DEFINE_IT(2,2,0)
DEFINE_IT(2,2,1)
#undef DEFINE_IT



#define DEFINE_IT(nn, zz, ff) \
static int CAT4(radix_4_moth_inv_trunc_block, nn, zz, ff)( \
    const sd_fft_ctx_t Q, \
    double* X0, double* X1, double* X2, double* X3, \
    ulong j, \
    ulong jm) \
{ \
    int l = 2; \
    flint_printf("function l = $d, n = %d, z = %d, f = %d", l, nn, zz, ff); \
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
DEFINE_IT(0,4,1)
DEFINE_IT(1,1,0)
DEFINE_IT(1,1,1)
DEFINE_IT(1,2,0)
DEFINE_IT(1,2,1)
DEFINE_IT(1,3,0)
DEFINE_IT(1,3,1)
DEFINE_IT(1,4,0)
DEFINE_IT(1,4,1)
DEFINE_IT(2,1,0)
DEFINE_IT(2,1,1)
DEFINE_IT(2,2,0)
DEFINE_IT(2,2,1)
DEFINE_IT(2,3,0)
DEFINE_IT(2,3,1)
DEFINE_IT(2,4,0)
DEFINE_IT(2,4,1)
DEFINE_IT(3,1,0)
DEFINE_IT(3,1,1)
DEFINE_IT(3,2,0)
DEFINE_IT(3,2,1)
DEFINE_IT(3,3,0)
DEFINE_IT(3,3,1)
DEFINE_IT(3,4,0)
DEFINE_IT(3,4,1)
DEFINE_IT(4,1,0)
DEFINE_IT(4,1,1)
DEFINE_IT(4,2,0)
DEFINE_IT(4,2,1)
DEFINE_IT(4,3,0)
DEFINE_IT(4,3,1)
DEFINE_IT(4,4,0)
DEFINE_IT(4,4,1)

#undef DEFINE_IT


void sd_ifft_basecase_8_0(const sd_fft_ctx_t Q, double* X, ulong j, ulong jm)
{
    flint_printf("not implemented");
    flint_abort();
}

void sd_ifft_basecase_8_1(const sd_fft_ctx_t Q, double* X, ulong j, ulong jm)
{
    flint_printf("not implemented");
    flint_abort();
}


/* parameter 1: j can be zero */
void sd_ifft_base_1(const sd_fft_ctx_t Q, ulong I, ulong j)
{
    ulong jm = j^(n_saturate_bits(j)>>1);
    if (j == 0)
        sd_ifft_basecase_8_1(Q, sd_fft_ctx_blk_index(Q, I), j, jm);
    else
        sd_ifft_basecase_8_0(Q, sd_fft_ctx_blk_index(Q, I), j, jm);
}

/* parameter 0: j cannot be zero */
void sd_ifft_base_0(const sd_fft_ctx_t Q, ulong I, ulong j)
{
    ulong jm = j^(n_saturate_bits(j)>>1);
    sd_ifft_basecase_8_0(Q, sd_fft_ctx_blk_index(Q, I), j, jm);
}

void sd_ifft_main_block(
    const sd_fft_ctx_t Q,
    ulong I, /* starting index */
    ulong S, /* stride */
    ulong k, /* BLK_SZ transforms each of length 2^k */
    ulong j)
{
    if (k > 2)
    {
        ulong k1 = k/2;
        ulong k2 = k - k1;

        // row ffts
        ulong l1 = n_pow2(k1);
        ulong b = 0; do {
            sd_ifft_main_block(Q, I + (b<<k2)*S, S, k2, (j<<k1) + b);
        } while (b++, b < l1);

        // column ffts
        ulong l2 = n_pow2(k2);
        ulong a = 0; do {
            sd_ifft_main_block(Q, I + a*S, S<<k2, k1, j);
        } while (a++, a < l2);

        return;
    }

    ulong jm = j^(n_saturate_bits(j)>>1);

    if (k == 2)
    {
        double* X0 = sd_fft_ctx_blk_index(Q, I + S*0);
        double* X1 = sd_fft_ctx_blk_index(Q, I + S*1);
        double* X2 = sd_fft_ctx_blk_index(Q, I + S*2);
        double* X3 = sd_fft_ctx_blk_index(Q, I + S*3);
        RADIX_4_REVERSE_PARAM(vec8d, Q, j, jm, 1)
        ulong i = 0; do {
            RADIX_4_REVERSE_MOTH(vec8d, X0+i, X1+i, X2+i, X3+i);
        } while(i += 8, i < BLK_SZ);
    }
    else if (k == 1)
    {
        double* X0 = sd_fft_ctx_blk_index(Q, I + S*0);
        double* X1 = sd_fft_ctx_blk_index(Q, I + S*1);
        RADIX_2_REVERSE_PARAM(vec8d, Q, j)
        ulong i = 0; do {
            RADIX_2_REVERSE_MOTH(vec8d, X0+i, X1+i);
        } while (i += 8, i < BLK_SZ);
    }
}

void sd_ifft_main(
    const sd_fft_ctx_t Q,
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
        sd_ifft_base_1(Q, I+S*0, 4*j+0);
        sd_ifft_base_0(Q, I+S*1, 4*j+1);
        sd_ifft_base_0(Q, I+S*2, 4*j+2);
        sd_ifft_base_0(Q, I+S*3, 4*j+3);
        sd_ifft_main_block(Q, I, S, 2, j);
    }
    else if (k == 1)
    {
        /* k1 = 1; k2 = 0 */
        sd_ifft_base_1(Q, I+S*0, 2*j+0);
        sd_ifft_base_0(Q, I+S*1, 2*j+1);
        sd_ifft_main_block(Q, I, S, 1, j);
    }
    else
    {
        sd_ifft_base_1(Q, I, j);
    }
}

void sd_ifft_trunc_block(
    const sd_fft_ctx_t Q,
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

        static int (*tab[5*4*2])(const sd_fft_ctx_t, double*, double*, double*, double*, ulong, ulong) =
            {IT(0,1,0),IT(0,1,1), IT(0,2,0),IT(0,2,1), IT(0,3,0),IT(0,3,1), IT(0,4,0),IT(0,4,1),
             IT(1,1,0),IT(1,1,1), IT(1,2,0),IT(1,2,1), IT(1,3,0),IT(1,3,1), IT(1,4,0),IT(1,4,1),
             IT(2,1,0),IT(2,1,1), IT(2,2,0),IT(2,2,1), IT(2,3,0),IT(2,3,1), IT(2,4,0),IT(2,4,1),
             IT(3,1,0),IT(3,1,1), IT(3,2,0),IT(3,2,1), IT(3,3,0),IT(3,3,1), IT(3,4,0),IT(3,4,1),
             IT(4,1,0),IT(4,1,1), IT(4,2,0),IT(4,2,1), IT(4,3,0),IT(4,3,1), IT(4,4,0),IT(4,4,1)};

        ulong jm = j^(n_saturate_bits(j)>>1);
        if (LOOKUP_IT(n,z,f)(Q, sd_fft_ctx_blk_index(Q, I+S*0),
                                sd_fft_ctx_blk_index(Q, I+S*1),
                                sd_fft_ctx_blk_index(Q, I+S*2),
                                sd_fft_ctx_blk_index(Q, I+S*3), j, jm))
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

        static void (*tab[3*2*2])(const sd_fft_ctx_t, double*, double*, ulong) =
            {IT(0,1,0),IT(0,1,1), IT(0,2,0),IT(0,2,1),
             IT(1,1,0),IT(1,1,1), IT(1,2,0),IT(1,2,1),
             IT(2,1,0),IT(2,1,1), IT(2,2,0),IT(2,2,1)};

        LOOKUP_IT(n,z,f)(Q, sd_fft_ctx_blk_index(Q, I+S*0),
                            sd_fft_ctx_blk_index(Q, I+S*1), j);

#undef LOOKUP_IT
#undef IT

        return;
    }
}


void sd_ifft_trunc(
    const sd_fft_ctx_t Q,
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
                   sd_ifft_base_1(Q, I+S*0, 4*j+0);
        if (n > 1) sd_ifft_base_0(Q, I+S*1, 4*j+1);
        if (n > 2) sd_ifft_base_0(Q, I+S*2, 4*j+2);
        if (n > 3) sd_ifft_base_0(Q, I+S*3, 4*j+3);
        sd_ifft_trunc_block(Q, I, S, 2, j, z, n, f);
        if (f) sd_ifft_trunc(Q, I + S*n, S, 0, 4*j+n, 1, 0, f);
        
    }
    else if (k == 1)
    {
                   sd_ifft_base_1(Q, I+S*0, 2*j+0);
        if (n > 1) sd_ifft_base_0(Q, I+S*1, 2*j+1);
        sd_ifft_trunc_block(Q, I, S, 1, j, z, n, f);
        if (f) sd_ifft_trunc(Q, I + S*n, S, 0, 2*j+n, 1, 0, f);
    }
    else
    {
        FLINT_ASSERT(!f);
        sd_ifft_base_1(Q, I, j);
    }
}

