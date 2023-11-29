/*
    Copyright (C) 2022 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdint.h>
#include <string.h>
#include "thread_support.h"
#include "ulong_extras.h"
#include "nmod.h"
#include "fft_small.h"
#include "crt_helpers.h"

void crt_data_init(crt_data_t C, ulong prime, ulong coeff_len, ulong nprimes)
{
    C->prime = prime;
    C->coeff_len = coeff_len;
    C->nprimes = nprimes;
    C->data = FLINT_ARRAY_ALLOC(nprimes*coeff_len + coeff_len + nprimes, ulong);
}

void crt_data_clear(crt_data_t C)
{
    flint_free(C->data);
}



/*
    need  ceil(64*bn/bits) <= prod_primes/2^(2*bits)
    i.e.  (64*bn+bits-1)/bits <= prod_primes/2^(2*bits)
           64*bn <= bits*prod_primes/2^(2*bits) - (bits-1)
*/
static ulong crt_data_find_bn_bound(const crt_data_t C, ulong bits)
{
    ulong bound = 0;
    ulong q = (2*bits)/FLINT_BITS;
    ulong r = (2*bits)%FLINT_BITS;
    ulong n = C->coeff_len;
    ulong i;
    ulong* x;
    TMP_INIT;

    TMP_START;
    x = TMP_ARRAY_ALLOC(n+1, ulong);

    x[n] = mpn_mul_1(x, crt_data_prod_primes(C), n, bits);

    if (q < n+1)
    {
        if (r > 0)
            mpn_rshift(x + q, x + q, n + 1 - q, r);

        if (!mpn_sub_1(x + q, x + q, n + 1 - q, bits - 1))
        {
            mpn_rshift(x + q, x + q, n + 1 - q, 6);
            bound = (x + q)[0];
            for (i = q + 1; i < n + 1; i++)
                if (x[i] != 0)
                    bound = -UWORD(1);
        }
    }

    TMP_END;
    return bound;
}

/*
    need  ceil(64*bn/bits) <= prod_primes/2^(2*bits)
    first try bits = (nbits(prod_primes) - nbits(bn))/2 then adjust
    also require bits > 64 for some applications below
*/
static ulong crt_data_find_bits(const crt_data_t C, ulong bn)
{
    ulong p_nbits = flint_mpn_nbits(crt_data_prod_primes(C), C->coeff_len);
    ulong bits = n_max(66, (p_nbits - n_nbits(bn))/2);

    if (bn > crt_data_find_bn_bound(C, bits))
    {
        do {
            bits -= 1;
        } while (bits > 65 && bn > crt_data_find_bn_bound(C, bits));
    }
    else
    {
        while (bn <= crt_data_find_bn_bound(C, bits + 1))
        {
            bits += 1;
        }
    }

    FLINT_ASSERT(bits > 64 && bn <= crt_data_find_bn_bound(C, bn));
    return bits;
}





#define aindex(i) (a[i])

#if 0
void slow_mpn_to_fft_easy(
    sd_fft_lctx_t Q,
    double* z,
    const uint32_t* a,
    ulong iq_stop_easy,
    ulong bits,
    const double* two_pow)
{
    vec8d p    = vec8d_set_d(Q->p);
    vec8d pinv = vec8d_set_d(Q->pinv);

    ulong s = BLK_SZ/8/32*bits;

    FLINT_ASSERT(64 < bits);

    __m256i tab[256];

    for (ulong iq = 0; iq < iq_stop_easy; iq++)
    {
        double* zI = sd_fft_ctx_blk_index(z, iq);
        const uint32_t* aa = a + iq*(BLK_SZ/32)*bits;
        ulong k = 0;

        /*
            transposing this block all at once first is faster than gathering
            it in the loop
        */
#if 1
        do {
            __m256i A, B, C, D, E, F, G, H, X0, X1, Y0, Y1, Z0, Z1, W0, W1, P, Q, R, S, T, U, V, W;

            A = _mm256_loadu_si256((const __m256i*) (aa + k + 0*s));
            B = _mm256_loadu_si256((const __m256i*) (aa + k + 1*s));
            C = _mm256_loadu_si256((const __m256i*) (aa + k + 2*s));
            D = _mm256_loadu_si256((const __m256i*) (aa + k + 3*s));
            E = _mm256_loadu_si256((const __m256i*) (aa + k + 4*s));
            F = _mm256_loadu_si256((const __m256i*) (aa + k + 5*s));
            G = _mm256_loadu_si256((const __m256i*) (aa + k + 6*s));
            H = _mm256_loadu_si256((const __m256i*) (aa + k + 7*s));

            X0 = _mm256_unpacklo_epi32(A, B); // a0 b0 | a1 b1 | a4 b4 | a5 b5
            Y0 = _mm256_unpacklo_epi32(C, D); // c0 d0 | c1 d1 | c4 d4 | c5 d5
            Z0 = _mm256_unpacklo_epi32(E, F); // e0 f0 | e1 f1 | e4 f4 | e5 f5
            W0 = _mm256_unpacklo_epi32(G, H); // g0 h0 | g1 h1 | g4 h4 | g5 h5
            X1 = _mm256_unpackhi_epi32(A, B); // a2 b2 | a3 b3 | a6 b6 | a7 b7
            Y1 = _mm256_unpackhi_epi32(C, D); // c2 d2 | c3 d3 | c6 d6 | c7 d7
            Z1 = _mm256_unpackhi_epi32(E, F); // e2 f2 | e3 f3 | e6 f6 | e7 f7
            W1 = _mm256_unpackhi_epi32(G, H); // g2 h2 | g3 h3 | g6 h6 | g7 h7

            P = _mm256_unpacklo_epi64(X0, Y0); // a0 b0 | c0 d0 | a4 b4 | c4 d4
            Q = _mm256_unpackhi_epi64(X0, Y0); // a1 b1 | c1 d1 | a5 b5 | c5 d5
            R = _mm256_unpacklo_epi64(Z0, W0); // e0 f0 | g0 h0 | e4 f4 | g4 h4
            S = _mm256_unpackhi_epi64(Z0, W0); // e1 f1 | g1 h1 | e5 f5 | g5 h5
            T = _mm256_unpacklo_epi64(X1, Y1); // a2 b2 | c2 d2 | a6 b6 | c6 d6
            U = _mm256_unpackhi_epi64(X1, Y1); // a3 b3 | c3 d3 | a7 b7 | c7 d7
            V = _mm256_unpacklo_epi64(Z1, W1); // e2 f2 | g2 h2 | e6 f6 | g6 h6
            W = _mm256_unpackhi_epi64(Z1, W1); // e3 f3 | g3 h3 | e7 f7 | g7 h7

            tab[k+0] = _mm256_permute2x128_si256(P, R, 0 + 16*2);
            tab[k+4] = _mm256_permute2x128_si256(P, R, 1 + 16*3);
            tab[k+1] = _mm256_permute2x128_si256(Q, S, 0 + 16*2);
            tab[k+5] = _mm256_permute2x128_si256(Q, S, 1 + 16*3);
            tab[k+2] = _mm256_permute2x128_si256(T, V, 0 + 16*2);
            tab[k+6] = _mm256_permute2x128_si256(T, V, 1 + 16*3);
            tab[k+3] = _mm256_permute2x128_si256(U, W, 0 + 16*2);
            tab[k+7] = _mm256_permute2x128_si256(U, W, 1 + 16*3);

        } while (k += 8, k + 8 <= s);
#endif

        while (k < s) {
            tab[k] = _mm256_set_epi32(aa[k + 7*s], aa[k + 6*s],
                                      aa[k + 5*s], aa[k + 4*s],
                                      aa[k + 3*s], aa[k + 2*s],
                                      aa[k + 1*s], aa[k + 0*s]);
            k++;
        }

        for (ulong ir = 0; ir < BLK_SZ/8; ir++)
        {
            ulong k = (ir*bits)/32;
            ulong j = (ir*bits)%32;

            __m256i AK = tab[k];/*_mm256_i32gather_epi32((const int*)aa+k, index, 4);*/
            __m256i AKJ = _mm256_srl_epi32(AK, _mm_set_epi32(0,0,0,j));
            vec8d ak = _vec8i32_convert_vec8d(AKJ);

            vec8d X = ak;
            k++;
            j = 32 - j;
            do {
                AK = tab[k];/*_mm256_i32gather_epi32((const int*)aa+k, index, 4);*/
                AKJ = AK;
                ak = _vec8i32_convert_vec8d(AKJ);
                X = vec8d_add(X, vec8d_mulmod(ak, vec8d_set_d(two_pow[j]), p, pinv));
                k++;
                j += 32;
            } while (j + 32 <= bits);

            if ((bits-j) != 0)
            {
                AK = tab[k];/*_mm256_i32gather_epi32((const int*)aa+k, index, 4);*/
                AKJ = _mm256_sll_epi32(AK, _mm_set_epi32(0,0,0,32-(bits-j)));
                ak = _vec8i32_convert_vec8d(AKJ);
                X = vec8d_add(X, vec8d_mulmod(ak, vec8d_set_d(two_pow[bits - 32]), p, pinv));
                k++;
            }

            X = vec8d_reduce_to_pm1n(X, p, pinv);

            /* _vec8i32_convert_vec8d make the Xs slightly out of order */
            zI[ir+0*BLK_SZ/8] = X.e1[0];
            zI[ir+1*BLK_SZ/8] = X.e1[1];
            zI[ir+4*BLK_SZ/8] = X.e1[2];
            zI[ir+5*BLK_SZ/8] = X.e1[3];
            zI[ir+2*BLK_SZ/8] = X.e2[0];
            zI[ir+3*BLK_SZ/8] = X.e2[1];
            zI[ir+6*BLK_SZ/8] = X.e2[2];
            zI[ir+7*BLK_SZ/8] = X.e2[3];
        }
    }
}

#else
void slow_mpn_to_fft_easy(
    sd_fft_lctx_t Q,
    double* z,
    const uint32_t* a,
    ulong iq_stop_easy,
    ulong bits,
    const double* two_pow)
{
    vec8d p    = vec8d_set_d(Q->p);
    vec8d pinv = vec8d_set_d(Q->pinv);

    for (ulong iq = 0; iq < iq_stop_easy; iq++)
    {
        double* zI = sd_fft_ctx_blk_index(z, iq);

        for (ulong ir = 0; ir < BLK_SZ/8; ir++)
        {
            ulong k = iq*(BLK_SZ/32)*bits + (ir*bits)/32;
            ulong j = (ir*bits)%32;

            vec8d ak = vec8d_set_d8(aindex(k+0*BLK_SZ/8/32*bits) >> j,
                                    aindex(k+1*BLK_SZ/8/32*bits) >> j,
                                    aindex(k+2*BLK_SZ/8/32*bits) >> j,
                                    aindex(k+3*BLK_SZ/8/32*bits) >> j,
                                    aindex(k+4*BLK_SZ/8/32*bits) >> j,
                                    aindex(k+5*BLK_SZ/8/32*bits) >> j,
                                    aindex(k+6*BLK_SZ/8/32*bits) >> j,
                                    aindex(k+7*BLK_SZ/8/32*bits) >> j);
            vec8d X = ak;
            k++;
            j = 32 - j;
            while (j + 32 <= bits)
            {
                ak = vec8d_set_d8(aindex(k+0*BLK_SZ/8/32*bits),
                                  aindex(k+1*BLK_SZ/8/32*bits),
                                  aindex(k+2*BLK_SZ/8/32*bits),
                                  aindex(k+3*BLK_SZ/8/32*bits),
                                  aindex(k+4*BLK_SZ/8/32*bits),
                                  aindex(k+5*BLK_SZ/8/32*bits),
                                  aindex(k+6*BLK_SZ/8/32*bits),
                                  aindex(k+7*BLK_SZ/8/32*bits));
                X = vec8d_add(X, vec8d_mulmod(ak, vec8d_set_d(two_pow[j]), p, pinv));
                k++;
                j += 32;
            }

            if ((bits-j) != 0)
            {
                ak = vec8d_set_d8(aindex(k+0*BLK_SZ/8/32*bits) << (32-(bits-j)),
                                  aindex(k+1*BLK_SZ/8/32*bits) << (32-(bits-j)),
                                  aindex(k+2*BLK_SZ/8/32*bits) << (32-(bits-j)),
                                  aindex(k+3*BLK_SZ/8/32*bits) << (32-(bits-j)),
                                  aindex(k+4*BLK_SZ/8/32*bits) << (32-(bits-j)),
                                  aindex(k+5*BLK_SZ/8/32*bits) << (32-(bits-j)),
                                  aindex(k+6*BLK_SZ/8/32*bits) << (32-(bits-j)),
                                  aindex(k+7*BLK_SZ/8/32*bits) << (32-(bits-j)));
                X = vec8d_add(X, vec8d_mulmod(ak, vec8d_set_d(two_pow[bits-32]), p, pinv));
            }

            X = vec8d_reduce_to_pm1n(X, p, pinv);

            zI[ir+0*BLK_SZ/8] = vec8d_get_index(X, 0);
            zI[ir+1*BLK_SZ/8] = vec8d_get_index(X, 1);
            zI[ir+2*BLK_SZ/8] = vec8d_get_index(X, 2);
            zI[ir+3*BLK_SZ/8] = vec8d_get_index(X, 3);
            zI[ir+4*BLK_SZ/8] = vec8d_get_index(X, 4);
            zI[ir+5*BLK_SZ/8] = vec8d_get_index(X, 5);
            zI[ir+6*BLK_SZ/8] = vec8d_get_index(X, 6);
            zI[ir+7*BLK_SZ/8] = vec8d_get_index(X, 7);
        }

    }
}
#endif

#undef aindex


void slow_mpn_to_fft(
    sd_fft_lctx_t Q,
    double* z, ulong ztrunc,
    const ulong* a_, ulong an_,
    ulong bits,
    const double* two_pow)
{
    const uint32_t* a = (const uint32_t*)(a_);
    ulong an = 2*an_;
    ulong iq, i_stop_easy, iq_stop_easy;

    /* the highest index read from two_pow is bits - 32 */
    FLINT_ASSERT(bits - 32 < MPN_CTX_TWO_POWER_TAB_SIZE);

    /* if i*bits + 32 < 32*an, then the index into a is always in bounds */
    i_stop_easy = n_min(ztrunc, (32*an - 33)/bits);
    iq_stop_easy = i_stop_easy/BLK_SZ;

    slow_mpn_to_fft_easy(Q, z, a, iq_stop_easy, bits, two_pow);

    /* now the hard ones */
    {
        vec1d p = vec1d_set_d(Q->p);
        vec1d pinv = vec1d_set_d(Q->pinv);
#define aindex(i) (((i) < an) ? a[i] : (uint32_t)(0))
        for (iq = iq_stop_easy; iq < ztrunc/BLK_SZ; iq++)
        {
            double* zI = sd_fft_ctx_blk_index(z, iq);

            for (ulong ir = 0; ir < BLK_SZ; ir++)
            {
                ulong k = ((iq*BLK_SZ+ir)*bits)/32;
                ulong j = ((iq*BLK_SZ+ir)*bits)%32;

                vec1d ak = vec1d_set_d(aindex(k) >> j);
                vec1d X = ak;
                k++;
                j = 32 - j;
                while (j + 32 <= bits)
                {
                    ak = vec1d_set_d(aindex(k));
                    X = vec1d_add(X, vec1d_mulmod(ak, two_pow[j], p, pinv));
                    k++;
                    j += 32;
                }

                if ((bits-j) != 0)
                {
                    ak = vec1d_set_d(aindex(k) << (32-(bits-j)));
                    X = vec1d_add(X, vec1d_mulmod(ak, two_pow[bits-32], p, pinv));
                }

                X = vec1d_reduce_to_pm1n(X, p, pinv);

                zI[ir] = X;
            }
        }
#undef aindex
    }
}





#define aindex(i) (((i) < an) ? a[i] : (uint32_t)(0))

#define N_CDIV(a, b) (((a) + (b) - 1) / (b))

#define DEFINE_IT(NP) \
FLINT_STATIC_NOINLINE void CAT(mpn_to_ffts_hard, NP)( \
    sd_fft_ctx_struct* Rffts, double* d, ulong dstride, \
    const uint32_t* a, ulong an, ulong atrunc, \
    const vec4d* two_pow, \
    ulong start_hard, ulong stop_hard, \
    ulong bits) \
{ \
    ulong np = NP; \
    ulong nvs = N_CDIV(NP, VEC_SZ); \
    vec4d X[N_CDIV(NP, VEC_SZ)]; \
    vec4d P[N_CDIV(NP, VEC_SZ)]; \
    vec4d PINV[N_CDIV(NP, VEC_SZ)]; \
 \
    for (ulong l = 0; l < nvs; l++) \
    { \
        P[l]    = vec4d_set_d4(Rffts[4*l+0].p, Rffts[4*l+1].p, Rffts[4*l+2].p, Rffts[4*l+3].p); \
        PINV[l] = vec4d_set_d4(Rffts[4*l+0].pinv, Rffts[4*l+1].pinv, Rffts[4*l+2].pinv, Rffts[4*l+3].pinv); \
    } \
 \
    for (ulong i = start_hard; i < stop_hard; i++) \
    { \
        ulong k = (i*bits)/32; \
        ulong j = (i*bits)%32; \
 \
        vec4d ak = vec4d_set_d((double)(aindex(k) >> j)); \
        for (ulong l = 0; l < nvs; l++) \
            X[l] = ak; \
        k++; \
        j = 32 - j; \
        while (j + 32 <= bits) \
        { \
            ak = vec4d_set_d((double)(aindex(k))); \
            for (ulong l = 0; l < nvs; l++) \
                X[l] = vec4d_add(X[l], vec4d_mulmod(ak, two_pow[j*nvs+l], P[l], PINV[l])); \
            k++; \
            j += 32; \
        } \
 \
        if ((bits-j) != 0) \
        { \
            ak = vec4d_set_d((double)(aindex(k) << (32-(bits-j)))); \
            for (ulong l = 0; l < nvs; l++) \
                X[l] = vec4d_add(X[l], vec4d_mulmod(ak, two_pow[(bits-32)*nvs+l], P[l], PINV[l])); \
        } \
 \
        for (ulong l = 0; l < nvs; l++) \
            X[l] = vec4d_reduce_to_pm1n(X[l], P[l], PINV[l]); \
 \
        for (ulong l = 0; l < np; l++) \
            sd_fft_ctx_set_index(d + l*dstride, i, vec4d_get_index(X[l/VEC_SZ], l%VEC_SZ)); \
    } \
 \
    for (ulong l = 0; l < np; l++) \
        for (ulong i = stop_hard; i < atrunc; i++) \
            sd_fft_ctx_set_index(d + l*dstride, i, 0.0); \
}

DEFINE_IT(4)
DEFINE_IT(5)
DEFINE_IT(6)
DEFINE_IT(7)
DEFINE_IT(8)
#undef DEFINE_IT
#undef aindex
#undef N_CDIV

#define CODE(ir) \
{ \
    ulong k = ((i+ir)*bits)/32; \
    ulong j = ((  ir)*bits)%32; \
 \
    vec4d ak = vec4d_set_d((double)(a[k] >> j)); \
    for (ulong l = 0; l < nvs; l++) \
        X[l] = ak; \
    k++; \
    j = 32 - j; \
    while (j + 32 <= bits) \
    { \
        ak = vec4d_set_d((double)(a[k])); \
        for (ulong l = 0; l < nvs; l++) \
            X[l] = vec4d_add(X[l], vec4d_mulmod(ak, two_pow[j*nvs+l], P[l], PINV[l])); \
        k++; \
        j += 32; \
    } \
 \
    if ((bits-j) != 0) \
    { \
        ak = vec4d_set_d((double)(a[k] << (32-(bits-j)))); \
        for (ulong l = 0; l < nvs; l++) \
            X[l] = vec4d_add(X[l], vec4d_mulmod(ak, two_pow[(bits-32)*nvs+l], P[l], PINV[l])); \
    } \
 \
    for (ulong l = 0; l < nvs; l++) \
        X[l] = vec4d_reduce_to_pm1n(X[l], P[l], PINV[l]); \
 \
    for (ulong l = 0; l < np; l++) \
        sd_fft_ctx_set_index(d + l*dstride, i+ir, vec4d_get_index(X[l/VEC_SZ], l%VEC_SZ)); \
}

#define N_CDIV(a, b) (((a) + (b) - 1) / (b))

/* The the l^th fft ctx Rffts[l] is expected to have data at d + l*dstride */
#define DEFINE_IT(NP, BITS) \
static void CAT3(mpn_to_ffts, NP, BITS)( \
    sd_fft_ctx_struct* Rffts, double* d, ulong dstride, \
    const ulong* a_, ulong an_, ulong atrunc, \
    const vec4d* two_pow, \
    ulong start_easy, ulong stop_easy, \
    ulong start_hard, ulong stop_hard) \
{ \
    ulong np = NP; \
    ulong bits = BITS; \
    ulong nvs = N_CDIV(NP, VEC_SZ); \
 \
    FLINT_ASSERT(bits >= FLINT_BITS); \
    FLINT_ASSERT(bits - 32 < MPN_CTX_TWO_POWER_TAB_SIZE); \
 \
    const uint32_t* a = (const uint32_t*)(a_); \
    ulong an = 2*an_; \
 \
    vec4d X[N_CDIV(NP, VEC_SZ)]; \
    vec4d P[N_CDIV(NP, VEC_SZ)]; \
    vec4d PINV[N_CDIV(NP, VEC_SZ)]; \
 \
    for (ulong l = 0; l < nvs; l++) \
    { \
        P[l]    = vec4d_set_d4(Rffts[4*l+0].p, Rffts[4*l+1].p, Rffts[4*l+2].p, Rffts[4*l+3].p); \
        PINV[l] = vec4d_set_d4(Rffts[4*l+0].pinv, Rffts[4*l+1].pinv, Rffts[4*l+2].pinv, Rffts[4*l+3].pinv); \
    } \
 \
    if ((bits % 8) == 0) \
    { \
        FLINT_ASSERT(start_easy % 4 == 0); \
        FLINT_ASSERT(stop_easy % 4 == 0); \
        for (ulong i = start_easy ; i < stop_easy; i += 4) \
        { \
            CODE(0);CODE(1);CODE(2);CODE(3); \
        } \
    } \
    else if ((bits % 4) == 0) \
    { \
        FLINT_ASSERT(start_easy % 8 == 0); \
        FLINT_ASSERT(stop_easy % 8 == 0); \
        for (ulong i = start_easy ; i < stop_easy; i += 8) \
        { \
            CODE(0);CODE(1);CODE(2);CODE(3); \
            CODE(4);CODE(5);CODE(6);CODE(7); \
        } \
    } \
    else if ((bits % 2) == 0) \
    { \
        FLINT_ASSERT(start_easy % 16 == 0); \
        FLINT_ASSERT(stop_easy % 16 == 0); \
        for (ulong i = start_easy ; i < stop_easy; i += 16) \
        { \
            CODE(0);CODE(1);CODE(2);CODE(3); \
            CODE(4);CODE(5);CODE(6);CODE(7); \
            CODE(8);CODE(9);CODE(10);CODE(11); \
            CODE(12);CODE(13);CODE(14);CODE(15); \
        } \
    } \
    else \
    { \
        FLINT_ASSERT(0); \
    } \
 \
    CAT(mpn_to_ffts_hard, NP)(Rffts, d, dstride, a, an, atrunc, two_pow, \
                              start_hard, stop_hard, bits); \
}

DEFINE_IT(4, 84)
DEFINE_IT(4, 88)
DEFINE_IT(4, 92)
DEFINE_IT(5,112)
DEFINE_IT(5,116)
DEFINE_IT(5,120)
DEFINE_IT(6,136)
DEFINE_IT(6,140)
DEFINE_IT(6,144)
DEFINE_IT(7,160)
DEFINE_IT(7,164)
DEFINE_IT(7,168)
DEFINE_IT(8,184)
DEFINE_IT(8,188)
DEFINE_IT(8,192)
#undef DEFINE_IT
#undef CODE
#undef aindex
#undef N_CDIV



#define DEFINE_IT(n, n_plus_1) \
FLINT_FORCE_INLINE void CAT(_add_to_answer_easy, n)(ulong z[], ulong r[], ulong zn, ulong toff, ulong tshift) \
{ \
    FLINT_ASSERT(zn > toff); \
    if (tshift == 0) \
    { \
        CAT(multi_add, n)(z + toff, r); \
    } \
    else \
    { \
        r[n] = r[n-1] >> (64-tshift); \
        for (ulong k = n; k >= 2; k--) \
            r[k-1] = (r[k-1] << (tshift)) | (r[k-2] >> (64-tshift)); \
        r[0] =  r[0] << (tshift); \
        CAT(multi_add, n_plus_1)(z + toff, r); \
    } \
} \
FLINT_FORCE_INLINE void CAT(_add_to_answer_hard, n)(ulong z[], ulong r[], ulong zn, ulong toff, ulong tshift) \
{ \
    FLINT_ASSERT(zn > toff); \
    if (tshift == 0) \
    { \
        if (zn - toff >= n) \
        { \
            CAT(multi_add, n)(z + toff, r); \
            return; \
        } \
    } \
    else \
    { \
        r[n] = r[n-1] >> (64-tshift); \
        for (ulong k = n; k >= 2; k--) \
            r[k-1] = (r[k-1] << (tshift)) | (r[k-2] >> (64-tshift)); \
        r[0] =  r[0] << (tshift); \
        if (zn - toff > n) \
        { \
            CAT(multi_add, n_plus_1)(z + toff, r); \
            return; \
        } \
    } \
    FLINT_ASSERT(zn - toff <= n); \
    mpn_add_n(z + toff, z + toff, r, zn - toff); \
}

DEFINE_IT(4, 5)
DEFINE_IT(5, 6)
DEFINE_IT(6, 7)
DEFINE_IT(7, 8)
#undef DEFINE_IT

typedef void (*from_ffts_func)(
    ulong* z, ulong zn, ulong zlen,
    sd_fft_ctx_struct* Rffts, double* d, ulong dstride,
    crt_data_struct* Rcrts,
    ulong bits,
    ulong start_easy, ulong stop_easy,
    ulong* overhang);

/*
    The "n" here is the limb count Rcrts[np-1].coeff_len, which is big enough
    to hold (product of primes)*(number of primes), so it can hold the
    intermediate dot products f[0]*x[0] + ... + f[np-1]*x[np-1]. The x[i] are
    single limb and the f[i] are of length "m". The number of primes is "np".

    The coefficient of X^i, 0 <= i < zlen needs to be reconstructed and added
    to the answer mpn (z, zn). This involves the limbs

       z[floor(i*bits/64)] ... z[floor(i*bits/64)+n]

    so is easy if floor(i*bits/64)+n < zn.

    The the l^th fft ctx Rffts[l] is expected to have data at d + l*dstride

    if overhang = NULL

        handle output coefficients from [start_easy, zlen)
        end_easy is still expected to be valid

    if overhang != NULL

        overhang has space for n words

        handle output coefficients from [start_easy, end_easy) where
        start_easy and stop_easy are divisible by BLK_SZ

        write to output words
        [start_easy*bits/64, stop_easy*bits/64) [overhang+0, overhang+n)
*/
#define DEFINE_IT(NP, N, M) \
static void CAT(_mpn_from_ffts, NP)( \
    ulong* z, ulong zn, ulong zlen, \
    sd_fft_ctx_struct* Rffts, double* d, ulong dstride, \
    crt_data_struct* Rcrts, \
    ulong bits, \
    ulong start_easy, ulong stop_easy, \
    ulong* overhang) \
{ \
    ulong np = NP; \
    ulong n = N; \
    ulong m = M; \
    ulong zn_start = start_easy*bits/64; \
    ulong zn_stop  = (overhang == NULL) ? zn : stop_easy*bits/64; \
 \
    FLINT_ASSERT(n == Rcrts[np-1].coeff_len); \
    FLINT_ASSERT(start_easy <= stop_easy); \
 \
    if (n == m + 1) \
    { \
        for (ulong l = 0; l < np; l++) { \
            FLINT_ASSERT(crt_data_co_prime(Rcrts + np - 1, l)[m] == 0); \
        } \
    } \
    else \
    { \
        FLINT_ASSERT(n == m); \
    } \
 \
    memset(z + zn_start, 0, (zn_stop - zn_start)*sizeof(ulong)); \
 \
    ulong Xs[BLK_SZ*NP]; \
 \
    if (overhang != NULL) \
    { \
        for (ulong i = 0; i < n; i++) \
            overhang[i] = 0; \
        if (start_easy >= stop_easy) \
            return; \
        stop_easy -= BLK_SZ; \
    } \
 \
    for (ulong i = start_easy; i < stop_easy; i += BLK_SZ) \
    { \
        _convert_block(Xs, Rffts, d, dstride, np, i/BLK_SZ); \
 \
        for (ulong j = 0; j < BLK_SZ; j += 1) \
        { \
            ulong r[N + 1]; \
            ulong t[N + 1]; \
            ulong l = 0; \
 \
            CAT3(_big_mul, N, M)(r, t, _crt_data_co_prime(Rcrts + np - 1, l, n), Xs[l*BLK_SZ + j]); \
            for (l++; l < np; l++) \
                CAT3(_big_addmul, N, M)(r, t, _crt_data_co_prime(Rcrts + np - 1, l, n), Xs[l*BLK_SZ + j]); \
 \
            CAT(_reduce_big_sum, N)(r, t, crt_data_prod_primes(Rcrts + np - 1)); \
 \
            ulong toff = ((i+j)*bits)/FLINT_BITS; \
            ulong tshift = ((i+j)*bits)%FLINT_BITS; \
 \
            FLINT_ASSERT(zn_stop > n + toff); \
 \
            CAT(_add_to_answer_easy, N)(z, r, zn_stop, toff, tshift); \
        } \
    } \
 \
    if (overhang != NULL) \
    { \
        ulong i = stop_easy; \
        _convert_block(Xs, Rffts, d, dstride, np, i/BLK_SZ); \
 \
        for (ulong j = 0; j < BLK_SZ; j += 1) \
        { \
            ulong r[N + 1]; \
            ulong t[N + 1]; \
            ulong l = 0; \
 \
            CAT3(_big_mul, N, M)(r, t, _crt_data_co_prime(Rcrts + np - 1, l, n), Xs[l*BLK_SZ + j]); \
            for (l++; l < np; l++) \
                CAT3(_big_addmul, N, M)(r, t, _crt_data_co_prime(Rcrts + np - 1, l, n), Xs[l*BLK_SZ + j]); \
 \
            CAT(_reduce_big_sum, N)(r, t, crt_data_prod_primes(Rcrts + np - 1)); \
 \
            ulong toff = ((i+j)*bits)/FLINT_BITS; \
            ulong tshift = ((i+j)*bits)%FLINT_BITS; \
 \
            if (n + toff < zn_stop) \
            { \
                CAT(_add_to_answer_easy, N)(z, r, zn_stop, toff, tshift); \
            } \
            else \
            { \
                if (tshift == 0) \
                { \
                    r[n] = 0; \
                } \
                else \
                { \
                    r[n] = r[n-1] >> (64-tshift); \
                    for (ulong k = n; k >= 2; k--) \
                        r[k-1] = (r[k-1] << (tshift)) | (r[k-2] >> (64-tshift)); \
                    r[0] =  r[0] << (tshift); \
                } \
                /* add zn_stop - toff words to the answer */ \
                /* and n + 1 + toff - zn_stop words to the overhang */ \
                unsigned char cf = 0; \
                ulong k = 0; \
                for (; k < zn_stop - toff; k++) \
                    cf = _addcarry_ulong(cf, z[toff + k], r[k], &z[toff + k]); \
                for (; k <= n; k++) \
                    cf = _addcarry_ulong(cf, overhang[k-(zn_stop-toff)], r[k], &overhang[k-(zn_stop-toff)]); \
            } \
        } \
    } \
    else \
    { \
        for (ulong i = stop_easy; i < zlen; i++) \
        { \
            ulong r[N + 1]; \
            ulong t[N + 1]; \
            ulong l = 0; \
            double xx = sd_fft_ctx_get_index(d + l*dstride, i); \
            ulong x = vec1d_reduce_to_0n(xx, Rffts[l].p, Rffts[l].pinv); \
 \
            CAT3(_big_mul, N, M)(r, t, crt_data_co_prime(Rcrts + np - 1, l), x); \
            for (l++; l < np; l++) \
            { \
                xx = sd_fft_ctx_get_index(d + l*dstride, i); \
                x = vec1d_reduce_to_0n(xx, Rffts[l].p, Rffts[l].pinv); \
                CAT3(_big_addmul, N, M)(r, t, crt_data_co_prime(Rcrts + np - 1, l), x); \
            } \
 \
            CAT(_reduce_big_sum, N)(r, t, crt_data_prod_primes(Rcrts + np - 1)); \
 \
            ulong toff = (i*bits)/FLINT_BITS; \
            ulong tshift = (i*bits)%FLINT_BITS; \
 \
            if (toff >= zn) \
                break; \
 \
            CAT(_add_to_answer_hard, N)(z, r, zn, toff, tshift); \
        } \
    } \
}

DEFINE_IT(4, 4, 3)
DEFINE_IT(5, 4, 4)
DEFINE_IT(6, 5, 4)
DEFINE_IT(7, 6, 5)
DEFINE_IT(8, 7, 6)
#undef DEFINE_IT

ulong next_fft_number(ulong p)
{
    ulong bits, l, q;
    bits = n_nbits(p);
    l = n_trailing_zeros(p - 1);
    q = p - (UWORD(2) << l);
    if (bits < 15)
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);
    if (n_nbits(q) == bits)
        return q;
    if (l < 5)
        return n_pow2(bits - 2) + 1;
    return n_pow2(bits) - n_pow2(l - 1) + 1;
}

/* fill x[i] = 2^i mod p for 0 <= i < len */
static void fill_slow_two_pow_tab(double* x, ulong len, double p, double pinv)
{
    double t = 1;
    x[0] = t;
    for (ulong i = 1; i < len; i++)
    {
        double q = vec1d_round(vec1d_mul(t, 2*pinv));
        t = vec1d_fnmadd(q, p, vec1d_add(t, t));
        x[i] = t;
    }
}

/*
    fill in  d[i*nvs + k/VEC_SZ][k%VEC_SZ] = 2^i mod Rffts[k].p
    for 0 <= k < VEC_SZ*nvs and 0 <= i < len.
*/
static void fill_vec_two_pow_tab(
    vec4d* x,
    sd_fft_ctx_struct* Rffts,
    ulong len,
    ulong nvs)
{
    ulong i, l;
    vec4d* ps;

    ps = (vec4d*) flint_aligned_alloc(32, 2*nvs*sizeof(vec4d));
    for (l = 0; l < nvs; l++)
    {
        /* just p */
        ps[2*l+0] = vec4d_set_d4(Rffts[4*l+0].p,
                                 Rffts[4*l+1].p,
                                 Rffts[4*l+2].p,
                                 Rffts[4*l+3].p);
        /* 2/p */
        ps[2*l+1] = vec4d_set_d4(Rffts[4*l+0].pinv,
                                 Rffts[4*l+1].pinv,
                                 Rffts[4*l+2].pinv,
                                 Rffts[4*l+3].pinv);
        ps[2*l+1] = vec4d_add(ps[2*l+1], ps[2*l+1]);
    }

    for (l = 0; l < nvs; l++)
        x[0*nvs + l] = vec4d_one();

    for (i = 1; i < len; i++)
    for (l = 0; l < nvs; l++)
    {
        vec4d t = x[(i-1)*nvs+l];
        vec4d p = ps[2*l+0];
        vec4d two_over_p = ps[2*l+1];
        vec4d q = vec4d_round(vec4d_mul(t, two_over_p));
        x[i*nvs+l] = vec4d_fnmadd(q, p, vec4d_add(t, t));
    }

    flint_aligned_free(ps);
}





void mpn_ctx_init(mpn_ctx_t R, ulong p)
{
    slong i;

    R->buffer = NULL;
    R->buffer_alloc = 0;

    for (i = 0; i < MPN_CTX_NCRTS; i++)
    {
        if (i > 0)
            p = next_fft_number(p);

        while (!n_is_prime(p))
            p = next_fft_number(p);

        /* ffts */
        sd_fft_ctx_init_prime(R->ffts + i, p);

        /* crts */
        if (i == 0)
        {
            crt_data_init(R->crts + 0, p, 1, 1);
            *crt_data_co_prime_red(R->crts + 0, 0) = 1;
            crt_data_co_prime(R->crts + 0, 0)[0] = 1;
            crt_data_prod_primes(R->crts + 0)[0] = p;
        }
        else
        {
            ulong pi;
            ulong len = R->crts[i - 1].coeff_len;
            ulong* t, * tt;

            t = FLINT_ARRAY_ALLOC(2*(len + 2), ulong);
            tt = t + (len + 2);

            t[len + 1] = 0;
            t[len] = mpn_mul_1(t, crt_data_prod_primes(R->crts + i - 1), len, p);

            /* leave enough room for (product of primes)*(number of primes) */
            len += 2;
            mpn_mul_1(tt, t, len, i + 1);
            while (tt[len - 1] == 0)
                len--;

            crt_data_init(R->crts + i, p, len, i + 1);

            /* set product of primes */
            flint_mpn_copyi(crt_data_prod_primes(R->crts + i), t, len);

            /* set cofactors */
            for (pi = 0; pi < i + 1; pi++)
            {
                ulong* cofac = crt_data_co_prime(R->crts + i, pi);
                mpn_divexact_1(cofac, t, len, R->crts[pi].prime);
                *crt_data_co_prime_red(R->crts + i, pi) =
                                      mpn_mod_1(cofac, len, R->crts[pi].prime);
            }

            flint_free(t);
        }
    }

    /* powers of two for slow mod */
    {
        ulong len = MPN_CTX_TWO_POWER_TAB_SIZE;
        double* x = FLINT_ARRAY_ALLOC(len*MPN_CTX_NCRTS, double);
        R->slow_two_pow_buffer = x;
        for (i = 0; i < MPN_CTX_NCRTS; i++)
        {
            R->slow_two_pow_tab[i] = x;
            fill_slow_two_pow_tab(x, len, R->ffts[i].p, R->ffts[i].pinv);
            x += len;
        }
    }

    /* powers of two for fast mod */
    {
        ulong len = MPN_CTX_TWO_POWER_TAB_SIZE;
        ulong max_nvs = n_cdiv(MPN_CTX_NCRTS, VEC_SZ);
        vec4d* x = (vec4d*) flint_aligned_alloc(32,
                                    max_nvs*(max_nvs + 1)/2*len*sizeof(vec4d));
        R->vec_two_pow_buffer = x;
        for (ulong nvs = 1; nvs <= max_nvs; nvs++)
        {
            R->vec_two_pow_tab[nvs - 1] = x;
            fill_vec_two_pow_tab(x, R->ffts, len, nvs);
            x += nvs*len;
        }
    }

    R->profiles_size = 0;

/*flint_printf("\n");*/
#define PUSH_PROFILE(np_, bits_, n, m) \
    i = R->profiles_size; \
    R->profiles[i].np        = np_; \
    R->profiles[i].bits      = bits_; \
    R->profiles[i].bn_bound  = crt_data_find_bn_bound(R->crts + np_ - 1, bits_); \
    R->profiles[i].to_ffts   = CAT3(mpn_to_ffts, np_, bits_); \
/*flint_printf("profile np = %wu, bits = %3wu, bn <= 0x%16wx\n", R->profiles[i].np, R->profiles[i].bits, R->profiles[i].bn_bound);*/ \
    R->profiles_size = i + 1;

    PUSH_PROFILE(4, 84, 4,3);
    PUSH_PROFILE(4, 88, 4,3);
    PUSH_PROFILE(4, 92, 4,3);
    PUSH_PROFILE(5,112, 4,4);
    PUSH_PROFILE(5,116, 4,4);
    PUSH_PROFILE(5,120, 4,4);
    PUSH_PROFILE(6,136, 5,4);
    PUSH_PROFILE(6,140, 5,4);
    PUSH_PROFILE(6,144, 5,4);
    PUSH_PROFILE(7,160, 6,5);
    PUSH_PROFILE(7,164, 6,5);
    PUSH_PROFILE(7,168, 6,5);
    PUSH_PROFILE(8,184, 7,6);
    PUSH_PROFILE(8,188, 7,6);
    PUSH_PROFILE(8,192, 7,6);

    FLINT_ASSERT(R->profiles_size <= MAX_NPROFILES);
}

#define VEC_SZ 4

void mpn_ctx_clear(mpn_ctx_t R)
{
    slong i;

    for (i = 0; i < MPN_CTX_NCRTS; i++)
    {
        sd_fft_ctx_clear(R->ffts + i);
        crt_data_clear(R->crts + i);
    }

    flint_free(R->slow_two_pow_buffer);
    flint_aligned_free(R->vec_two_pow_buffer);

    flint_aligned_free(R->buffer);
}


typedef struct {
    thread_pool_handle* handles;
    slong nhandles;
    ulong nthreads;
    ulong np;
    ulong bits;
    to_ffts_func to_ffts;
} profile_entry;

static void mpn_ctx_best_profile(
    const mpn_ctx_t R,
    profile_entry* P,
    ulong an, ulong bn)
{
    ulong i = 0;
    ulong best_i = 0;
    double best_score = 100000000.0*(an + bn);

    ulong thread_limit = 8;
    ulong zn = an + bn;

    if (zn < 2048)
        thread_limit = 1;
    else if (zn < 4096)
        thread_limit = 4;
    else if (zn < 8192)
        thread_limit = 5;
    else if (zn < 16384)
        thread_limit = 6;
    else if (zn < 32768)
        thread_limit = 7;

    P->nhandles = flint_request_threads(&P->handles, thread_limit);
    P->nthreads = 1 + P->nhandles;

    /*
        The first profile is supposed to have the biggest bn_bound. If the
        given bn is too large, we must fill in P->to_ffts = NULL because we
        don't have a fast mod function.
        We can also fill in P->to_ffts = NULL any time to not use the fast mod
        function and use the slow generic one instead.
    */

    if (FLINT_UNLIKELY(bn > R->profiles[i].bn_bound))
    {
        P->np = n_max(n_min(P->nthreads, 8), 4);
        P->bits = crt_data_find_bits(R->crts + P->np - 1, bn);
        P->to_ffts = NULL;
        return;
    }

got_one:

    /* maximize R->profiles[i].bits */

    FLINT_ASSERT(i < R->profiles_size);
    FLINT_ASSERT(bn <= R->profiles[i].bn_bound);

    while (i+1 < R->profiles_size &&
           bn <= R->profiles[i+1].bn_bound &&
           R->profiles[i+1].np == R->profiles[i].np)
    {
        i++;
    }

    ulong np = R->profiles[i].np;

    if (np % P->nthreads != 0)
        goto find_next;

    ulong bits = R->profiles[i].bits;
    ulong alen = n_cdiv(64*an, bits);
    ulong blen = n_cdiv(64*bn, bits);
    ulong zlen = alen + blen - 1;
    ulong ztrunc = n_round_up(zlen, BLK_SZ);
    ulong depth = n_max(LG_BLK_SZ, n_clog2(ztrunc));

    double ratio = (double)(ztrunc)/(double)(n_pow2(depth));
    double score = (1-0.25*ratio)*(1.0/1000000);
    score *= np*depth;
    score *= ztrunc;
    if (score < best_score)
    {
        best_i = i;
        best_score = score;
    }

find_next:

    do {
        i++;
        if (i >= R->profiles_size)
        {
            P->np = R->profiles[best_i].np;
            P->bits = R->profiles[best_i].bits;
            P->to_ffts = R->profiles[best_i].to_ffts;
            return;
        }
    } while (bn > R->profiles[i].bn_bound);

    goto got_one;
}

void* mpn_ctx_fit_buffer(mpn_ctx_t R, ulong n)
{
    if (n > R->buffer_alloc)
    {
        flint_aligned_free(R->buffer);
        n = n_round_up(n_max(n, R->buffer_alloc*17/16), 4096);
        R->buffer = flint_aligned_alloc(4096, n);
        R->buffer_alloc = n;
    }
    return R->buffer;
}

/* pointwise mul of a with b and m */
void sd_fft_lctx_point_mul(
    const sd_fft_lctx_t Q,
    double* a,
    const double* b,
    ulong m_,
    ulong depth)
{
    vec8d m = vec8d_set_d(vec1d_reduce_0n_to_pmhn((slong)m_, Q->p));
    vec8d n    = vec8d_set_d(Q->p);
    vec8d ninv = vec8d_set_d(Q->pinv);
    FLINT_ASSERT(depth >= LG_BLK_SZ);
    for (ulong I = 0; I < n_pow2(depth - LG_BLK_SZ); I++)
    {
        double* ax = a + sd_fft_ctx_blk_offset(I);
        const double* bx = b + sd_fft_ctx_blk_offset(I);
        ulong j = 0; do {
            vec8d x0, x1, b0, b1;
            x0 = vec8d_load(ax+j+0);
            x1 = vec8d_load(ax+j+8);
            b0 = vec8d_load(bx+j+0);
            b1 = vec8d_load(bx+j+8);
            x0 = vec8d_mulmod(x0, m, n, ninv);
            x1 = vec8d_mulmod(x1, m, n, ninv);
            x0 = vec8d_mulmod(x0, b0, n, ninv);
            x1 = vec8d_mulmod(x1, b1, n, ninv);
            vec8d_store(ax+j+0, x0);
            vec8d_store(ax+j+8, x1);
        } while (j += 16, j < BLK_SZ);
    }
}

void sd_fft_lctx_point_sqr(
    const sd_fft_lctx_t Q,
    double* a,
    ulong m_,
    ulong depth)
{
    vec8d m = vec8d_set_d(vec1d_reduce_0n_to_pmhn((slong)m_, Q->p));
    vec8d n    = vec8d_set_d(Q->p);
    vec8d ninv = vec8d_set_d(Q->pinv);
    FLINT_ASSERT(depth >= LG_BLK_SZ);

    for (ulong I = 0; I < n_pow2(depth - LG_BLK_SZ); I++)
    {
        double* ax = a + sd_fft_ctx_blk_offset(I);
        ulong j = 0; do {
            vec8d x0, x1;
            x0 = vec8d_load(ax+j+0);
            x1 = vec8d_load(ax+j+8);
            x0 = vec8d_mulmod(x0, x0, n, ninv);
            x1 = vec8d_mulmod(x1, x1, n, ninv);
            x0 = vec8d_mulmod(x0, m, n, ninv);
            x1 = vec8d_mulmod(x1, m, n, ninv);
            vec8d_store(ax+j+0, x0);
            vec8d_store(ax+j+8, x1);
        } while (j += 16, j < BLK_SZ);
    }
}

typedef struct {
    to_ffts_func to_ffts;
    sd_fft_ctx_struct* ffts;
    ulong stride;
    const vec4d* two_pow_tab;
    double* abuf;
    const ulong* a;
    ulong an;
    ulong atrunc;
    ulong a_start_easy;
    ulong a_stop_easy;
    ulong a_start_hard;
    ulong a_stop_hard;
    double* bbuf;
    const ulong* b;
    ulong bn;
    ulong btrunc;
    ulong b_start_easy;
    ulong b_stop_easy;
    ulong b_start_hard;
    ulong b_stop_hard;
    int squaring;
} mod_worker_struct;

void mod_worker_func(void* varg)
{
    mod_worker_struct* X = (mod_worker_struct*) varg;

    X->to_ffts(X->ffts, X->abuf, X->stride, X->a, X->an, X->atrunc, X->two_pow_tab,
             X->a_start_easy, X->a_stop_easy, X->a_start_hard, X->a_stop_hard);

    if (!X->squaring)
    {
        X->to_ffts(X->ffts, X->bbuf, X->stride, X->b, X->bn, X->btrunc, X->two_pow_tab,
                 X->b_start_easy, X->b_stop_easy, X->b_start_hard, X->b_stop_hard);
    }
}

typedef struct fft_worker_struct {
    sd_fft_ctx_struct* fctx;
    ulong cop;
    ulong depth;
    ulong ztrunc;
    double* abuf;
    ulong atrunc;
    double* bbuf;
    ulong btrunc;
    struct fft_worker_struct* next;
    int squaring;
} fft_worker_struct;

void fft_worker_func(void* varg)
{
    fft_worker_struct* X = (fft_worker_struct*) varg;
    sd_fft_lctx_t Q;
    ulong m;

    do {
        sd_fft_lctx_init(Q, X->fctx, X->depth);

        if (!X->squaring)
            sd_fft_lctx_fft_trunc(Q, X->bbuf, X->depth, X->btrunc, X->ztrunc);

        sd_fft_lctx_fft_trunc(Q, X->abuf, X->depth, X->atrunc, X->ztrunc);
        NMOD_RED2(m, X->cop >> (64 - X->depth), X->cop << X->depth, X->fctx->mod);
        m = nmod_inv(m, X->fctx->mod);

        if (X->squaring)
            sd_fft_lctx_point_sqr(Q, X->abuf, m, X->depth);
        else
            sd_fft_lctx_point_mul(Q, X->abuf, X->bbuf, m, X->depth);

        sd_fft_lctx_ifft_trunc(Q, X->abuf, X->depth, X->ztrunc);
        sd_fft_lctx_clear(Q, X->fctx);
    } while (X = X->next, X != NULL);
}

typedef struct mod_fft_worker_struct {
    ulong bits;
    sd_fft_ctx_struct* fctx;
    const double* two_pow_tab;
    ulong cop;
    ulong depth;
    ulong ztrunc;
    const ulong* a;
    ulong an;
    double* abuf;
    ulong atrunc;
    const ulong* b;
    ulong bn;
    double* bbuf;
    ulong btrunc;
    struct mod_fft_worker_struct* next;
    int squaring;
} mod_fft_worker_struct;

void mod_fft_worker_func(void* varg)
{
    mod_fft_worker_struct* X = (mod_fft_worker_struct*) varg;
    sd_fft_lctx_t Q;
    ulong m;

    do {
        sd_fft_lctx_init(Q, X->fctx, X->depth);

        if (!X->squaring)
        {
            slow_mpn_to_fft(Q, X->bbuf, X->btrunc, X->b, X->bn, X->bits, X->two_pow_tab);
            sd_fft_lctx_fft_trunc(Q, X->bbuf, X->depth, X->btrunc, X->ztrunc);
        }

        slow_mpn_to_fft(Q, X->abuf, X->atrunc, X->a, X->an, X->bits, X->two_pow_tab);
        sd_fft_lctx_fft_trunc(Q, X->abuf, X->depth, X->atrunc, X->ztrunc);

        NMOD_RED2(m, X->cop >> (64 - X->depth), X->cop << X->depth, X->fctx->mod);
        m = nmod_inv(m, X->fctx->mod);

        if (X->squaring)
            sd_fft_lctx_point_sqr(Q, X->abuf, m, X->depth);
        else
            sd_fft_lctx_point_mul(Q, X->abuf, X->bbuf, m, X->depth);

        sd_fft_lctx_ifft_trunc(Q, X->abuf, X->depth, X->ztrunc);
        sd_fft_lctx_clear(Q, X->fctx);
    } while (X = X->next, X != NULL);
}

typedef struct {
    from_ffts_func from_ffts;
    ulong* z;
    ulong zn;
    ulong zlen;
    sd_fft_ctx_struct* fctxs;
    double* abuf;
    ulong stride;
    crt_data_struct* crts;
    ulong bits;
    ulong start_easy;
    ulong stop_easy;
    ulong* overhang;
    ulong overhang_buffer[MPN_CTX_NCRTS];
} crt_worker_struct;

void crt_worker_func(void* varg)
{
    crt_worker_struct* X = (crt_worker_struct*) varg;

    X->from_ffts(X->z, X->zn, X->zlen, X->fctxs, X->abuf, X->stride, X->crts,
                X->bits, X->start_easy, X->stop_easy, X->overhang);
}


void mpn_ctx_mpn_mul(mpn_ctx_t R, ulong* z, const ulong* a, ulong an, const ulong* b, ulong bn)
{
    ulong zn, alen, blen, zlen, atrunc, btrunc, ztrunc, depth, stride;
    double* abuf;
    profile_entry P;
    ulong sz;
    void* worker_struct_buffer;
    int squaring;

    mpn_ctx_best_profile(R, &P, an, bn);

    sz =           sizeof(mod_worker_struct)*P.nthreads;
    sz = n_max(sz, sizeof(fft_worker_struct)*P.np);
    sz = n_max(sz, sizeof(mod_fft_worker_struct)*P.np);
    sz = n_max(sz, sizeof(crt_worker_struct)*P.nthreads);
    worker_struct_buffer = flint_malloc(sz);

    squaring = (a == b) && (an == bn);
    zn = an + bn;
    alen = n_cdiv(FLINT_BITS*an, P.bits);
    blen = n_cdiv(FLINT_BITS*bn, P.bits);
    zlen = alen + blen - 1;
    atrunc = n_round_up(alen, BLK_SZ);
    btrunc = n_round_up(blen, BLK_SZ);
    ztrunc = n_round_up(zlen, BLK_SZ);
    depth = n_max(LG_BLK_SZ, n_clog2(ztrunc));
    stride = n_round_up(sd_fft_ctx_data_size(depth), 128);

    FLINT_ASSERT(an > 0);
    FLINT_ASSERT(bn > 0);
    FLINT_ASSERT(0 <= flint_mpn_cmp_ui_2exp(
                                crt_data_prod_primes(R->crts + P.np - 1),
                                R->crts[P.np - 1].coeff_len, blen, 2*P.bits));

#define TIME_THIS 0

#if TIME_THIS
timeit_t timer, timer_overall;
flint_printf("------------ zn = %wu, nthreads = %wu np = %wu, bits = %wu, -------------\n", zn, nthreads, np, bits);
#endif

#if TIME_THIS
timeit_start(timer_overall);
#endif

    if (P.to_ffts != NULL)
    {
        ulong bits = P.bits;
        mod_worker_struct* wm;
        fft_worker_struct* wf;
        /* if i*bits + 32 < 64*an, then the index into a is always in bounds */
        ulong a_stop_easy = n_min(atrunc, (64*an - 33)/bits);
        /* if i*bits >= 64*an, then the index into a is always out of bounds */
        ulong a_stop_hard = n_min(atrunc, (64*an + bits - 1)/bits);
        /* ditto */
        ulong b_stop_easy = n_min(btrunc, (64*bn - 33)/bits);
        ulong b_stop_hard = n_min(btrunc, (64*bn + bits - 1)/bits);
        ulong rounding = (bits%8 == 0) ? 4 : (bits%4 == 0) ? 8 : 16;
        ulong nthreads = P.nthreads;
        double* bbuf;

        abuf = (double*) mpn_ctx_fit_buffer(R, 2*P.np*stride*sizeof(double));
        bbuf = abuf + P.np*stride;

#if TIME_THIS
timeit_start(timer);
#endif

        /* some fixups for loop unrollings: round down the easy stops */
        FLINT_ASSERT(bits%2 == 0);
        a_stop_easy &= -rounding;
        b_stop_easy &= -rounding;

        wm = (mod_worker_struct*) worker_struct_buffer;
        for (ulong i = 0; i < nthreads; i++)
        {
            mod_worker_struct* X = wm + i;
            X->to_ffts = P.to_ffts;
            X->ffts = R->ffts;
            X->stride = stride;
            X->two_pow_tab = R->vec_two_pow_tab[n_cdiv(P.np, VEC_SZ) - 1];
            X->abuf = abuf;
            X->a = a;
            X->an = an,
            X->atrunc = atrunc;
            X->bbuf = bbuf;
            X->b = b;
            X->bn = bn,
            X->btrunc = btrunc;
            X->a_start_easy = n_round_up((i+0)*a_stop_easy/nthreads, rounding);
            X->a_stop_easy  = n_round_up((i+1)*a_stop_easy/nthreads, rounding);
            X->b_start_easy = n_round_up((i+0)*b_stop_easy/nthreads, rounding);
            X->b_stop_easy  = n_round_up((i+1)*b_stop_easy/nthreads, rounding);
            /* only the last thread i = nthreads - 1 does the hard ends */
            X->a_start_hard = (i + 1 == nthreads) ? a_stop_easy : atrunc;
            X->a_stop_hard  = (i + 1 == nthreads) ? a_stop_hard : atrunc;
            X->b_start_hard = (i + 1 == nthreads) ? b_stop_easy : btrunc;
            X->b_stop_hard  = (i + 1 == nthreads) ? b_stop_hard : btrunc;
            X->squaring = squaring;
        }

        for (slong i = P.nhandles; i > 0; i--)
            thread_pool_wake(global_thread_pool, P.handles[i - 1], 0,
                                                       mod_worker_func, wm + i);
        mod_worker_func(wm + 0);

        for (slong i = P.nhandles; i > 0; i--)
            thread_pool_wait(global_thread_pool, P.handles[i - 1]);

#if TIME_THIS
timeit_stop(timer);
if (timer->wall > 50)
flint_printf("    mod: %wd\n", timer->wall);
#endif

#if TIME_THIS
timeit_start(timer);
#endif

        /*
            current scheduling:
                np = 5, nthreads = 3:
                thread0: p0, p3
                thread1: p1, p4
                thread2: p2

                np = 3, nthreads = 5:
                thread0: p0
                thread1: p1
                thread2: p2
                thread3: -
                thread4: -
        */

        wf = (fft_worker_struct*) worker_struct_buffer;

        for (ulong l = 0; l < P.np; l++)
        {
            fft_worker_struct* X = wf + l;
            X->fctx = R->ffts + l;
            X->cop = *crt_data_co_prime_red(R->crts + P.np - 1, l);
            X->depth = depth;
            X->ztrunc = ztrunc;
            X->abuf = abuf + l*stride;
            X->atrunc = atrunc;
            X->bbuf = bbuf + l*stride;
            X->btrunc = btrunc;
            X->next = (l + nthreads < P.np) ? X + nthreads : NULL;
            X->squaring = squaring;
        }

        for (ulong i = n_min(P.nhandles, P.np - 1); i > 0; i--)
            thread_pool_wake(global_thread_pool, P.handles[i - 1], 0,
                                                      fft_worker_func, wf + i);
        fft_worker_func(wf + 0);

        for (ulong i = n_min(P.nhandles, P.np - 1); i > 0; i--)
            thread_pool_wait(global_thread_pool, P.handles[i - 1]);

#if TIME_THIS
timeit_stop(timer);
if (timer->wall > 50)
flint_printf("    fft: %wd\n", timer->wall);
#endif
    }
    else
    {
        mod_fft_worker_struct* w = (mod_fft_worker_struct*) worker_struct_buffer;
        ulong bits = P.bits;
        ulong np = P.np;
        ulong nthreads = P.nthreads;
        double* bbuf;

        abuf = (double*) mpn_ctx_fit_buffer(R, (np+nthreads)*stride*sizeof(double));
        bbuf = abuf + np*stride;

#if TIME_THIS
timeit_start(timer);
#endif
        for (ulong l = 0; l < np; l++)
        {
            mod_fft_worker_struct* X = w + l;
            X->bits = bits;
            X->fctx = R->ffts + l;
            X->cop = *crt_data_co_prime_red(R->crts + np - 1, l);
            X->depth = depth;
            X->ztrunc = ztrunc;
            X->a = a;
            X->an = an;
            X->abuf = abuf + l*stride;
            X->atrunc = atrunc;
            X->b = b;
            X->bn = bn;
            X->bbuf = bbuf + (l%nthreads)*stride;
            X->btrunc = btrunc;
            X->two_pow_tab = R->slow_two_pow_tab[l];
            X->next = (l + nthreads < np) ? X + nthreads : NULL;
            X->squaring = squaring;
        }

        for (ulong i = n_min(P.nhandles, P.np - 1); i > 0; i--)
            thread_pool_wake(global_thread_pool, P.handles[i - 1], 0,
                                                   mod_fft_worker_func, w + i);
        mod_fft_worker_func(w + 0);

        for (ulong i = n_min(P.nhandles, P.np - 1); i > 0; i--)
            thread_pool_wait(global_thread_pool, P.handles[i - 1]);

#if TIME_THIS
timeit_stop(timer);
if (timer->wall > 50)
flint_printf("mod+fft: %wd\n", timer->wall);
#endif
    }

#if TIME_THIS
timeit_start(timer);
#endif

    {
        ulong n = R->crts[P.np-1].coeff_len;
        crt_worker_struct* w = (crt_worker_struct*) worker_struct_buffer;
        ulong nthreads = P.nthreads;
        ulong end_easy = (zn >= n+1 ? zn - (n+1) : UWORD(0))*64/P.bits;

        /* this is how must space was statically allocated in each struct */
        FLINT_ASSERT(n <= MPN_CTX_NCRTS);

        end_easy &= -BLK_SZ;

        FLINT_ASSERT(4 <= P.np && P.np <= 8);
        static from_ffts_func tab[8-4+1] = {_mpn_from_ffts_4,
                                            _mpn_from_ffts_5,
                                            _mpn_from_ffts_6,
                                            _mpn_from_ffts_7,
                                            _mpn_from_ffts_8};

        for (ulong l = 0; l < nthreads; l++)
        {
            crt_worker_struct* X = w + l;
            X->from_ffts = tab[P.np - 4];
            X->z = z;
            X->zn = zn;
            X->zlen = zlen;
            X->fctxs = R->ffts;
            X->abuf = abuf;
            X->stride = stride;
            X->crts = R->crts;
            X->bits = P.bits;
            X->start_easy = n_round_up((l+0)*end_easy/nthreads, BLK_SZ);
            X->stop_easy  = n_round_up((l+1)*end_easy/nthreads, BLK_SZ);
            X->overhang = (l + 1 == nthreads) ? NULL : X->overhang_buffer;
        }

        for (slong i = P.nhandles; i > 0; i--)
            thread_pool_wake(global_thread_pool, P.handles[i - 1], 0,
                                                   crt_worker_func, w + i);
        crt_worker_func(w + 0);

        for (slong i = P.nhandles; i > 0; i--)
            thread_pool_wait(global_thread_pool, P.handles[i - 1]);

        unsigned char cf = 0;
        for (ulong i = 1; i <= P.nhandles; i++)
        {
            ulong start = w[i].start_easy*P.bits/64;
            if (i == P.nhandles)
            {
                cf = flint_mpn_add_inplace_c(z + start, zn - start,
                                              w[i - 1].overhang_buffer, n, cf);
            }
            else
            {
                ulong stop = w[i].stop_easy*P.bits/64;
                if (stop > start)
                {
                    cf = flint_mpn_add_inplace_c(z + start, stop - start,
                                              w[i - 1].overhang_buffer, n, cf);
                }
                else
                {
                    for (ulong k = 0; k < n; k++)
                    {
                        FLINT_ASSERT(w[i].overhang_buffer[k] == 0);
                        w[i].overhang_buffer[k] = w[i - 1].overhang_buffer[k];
                    }
                }
            }
        }
    }

#if TIME_THIS
timeit_stop(timer);
if (timer->wall > 50)
flint_printf("    crt: %wd\n", timer->wall);
timeit_stop(timer_overall);
if (timer_overall->wall > 50)
flint_printf("      +: %wd\n", timer_overall->wall);
#endif

#undef TIME_THIS

    flint_free(worker_struct_buffer);
    flint_give_back_threads(P.handles, P.nhandles);
}

