/*
    Copyright (C) 2026 Fredrik Johansson

    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Negacyclic pointwise multiplication: z = x*y mod 2^N+1 on
    top-limb-0-or-1 (n+1)-limb operands, via CRT over sd_fft primes
    with b-bit digits, N = m*b, m a power of two. Pointwise stage
    for Schoenhage-Strassen convolutions (fmpz_poly mul_SS).
*/

#include <string.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod.h"
#include "fmpz.h"
#include "mpn_extras.h"
#include "fft_small.h"
#include "crt_helpers.h"
#include "longlong.h"

#define LG_BLK 7   /* LG_BLK_SZ */



/* Smallest np whose exact prime product covers the biased
   coefficient bound 2^(2b + depth + 1), among those with an
   instantiated input pass. Least np is fastest: every stage
   (input, transforms, pointwise, crt) scales with np. */
static slong _mulmod_2expp1_choose_m(mpn_ctx_struct * R,
                                     slong N);

slong
sd_fft_mpn_mulmod_2expp1_choose_np(mpn_ctx_struct * R, slong b, slong depth)
{
    fmpz_t P;
    slong np, need = 2 * b + depth + 1;

    fmpz_init_set_ui(P, R->ffts[0].p);
    for (np = 2; np <= 8; np++)
    {
        fmpz_mul_ui(P, P, R->ffts[np - 1].p);
        if ((slong) fmpz_bits(P) - 1 >= need
            && (b <= 50
                || _mpn_ctx_to_ffts_func((ulong) np, (ulong) b) != NULL))
        {
            fmpz_clear(P);
            return np;
        }
    }
    fmpz_clear(P);
    return -1;
}

void
sd_fft_mpn_mulmod_2expp1_ctx_init(sd_fft_mpn_mulmod_2expp1_ctx_struct * C, mpn_ctx_struct * R,
                    slong N, slong m)
{
    slong np;

    if (m == 0)
    {
        m = _mulmod_2expp1_choose_m(R, N);
        if (m == 0)
            flint_throw(FLINT_ERROR,
                "sd_fft_mpn_mulmod_2expp1_ctx_init: no supported "
                "digit splitting for N = %wd\n", N);
    }
    if (N % m != 0 || m < 16 || (m & (m - 1)) != 0)
        flint_throw(FLINT_ERROR, "sd_fft_mpn_mulmod_2expp1_ctx_init: need N = m*b "
                    "with m a power of two >= 16\n");
    np = sd_fft_mpn_mulmod_2expp1_choose_np(R, N / m, FLINT_BIT_COUNT(m) - 1);
    if (np < 0)
        flint_throw(FLINT_ERROR, "sd_fft_mpn_mulmod_2expp1_ctx_init: no supported "
                    "prime configuration for b = %wd, m = %wd\n",
                    N / m, m);
    C->R = R;
    C->N = N;
    C->m = m;
    C->b = N / m;
    C->np = np;
    C->depth = FLINT_BIT_COUNT(m) - 1;
    C->dsz = (slong) sd_fft_ctx_data_size(C->depth);
    C->biaspoly = NULL;
    for (slong k = 0; k < np; k++)
    {
        C->p[k] = C->R->ffts[k].p;
        nmod_init(&C->mod[k], C->p[k]);
        sd_fft_ctx_fit_depth(C->R->ffts + k, C->depth);
    }

    for (slong k = 0; k < np; k++)
    {
        /* 2m-th root: g^((p-1)/2m) */
        ulong g = n_primitive_root_prime(C->p[k]);
        ulong w = nmod_pow_ui(g, (C->p[k] - 1) / (2 * (ulong) m),
                              C->mod[k]);
        ulong minv = nmod_inv(m, C->mod[k]);
        /* the from_ffts contract (determined empirically with a
           planted-coefficient probe): stored residues carry the
           inverse of the reduced crt cofactor. The two-prime
           reconstruction is a plain Garner step instead, so its
           untwist stays cofactor-free and the export saves the
           corresponding multiplications. */
        if (C->np != 2)
            minv = nmod_mul(minv, nmod_inv(
                *crt_data_co_prime_red(C->R->crts + C->np - 1, k) % C->p[k],
                C->mod[k]), C->mod[k]);
        ulong winv = nmod_inv(w, C->mod[k]);
        C->tw[k] = flint_aligned_alloc(FLINT_FFT_SMALL_ALIGNMENT,
                        n_round_up(m * sizeof(double), FLINT_FFT_SMALL_ALIGNMENT));
        C->itw[k] = flint_aligned_alloc(FLINT_FFT_SMALL_ALIGNMENT,
                        n_round_up(m * sizeof(double), FLINT_FFT_SMALL_ALIGNMENT));
        ulong wp = 1, wip = minv;
        for (slong i = 0; i < m; i++)
        {
            C->tw[k][i] = (double)(slong)(wp <= C->p[k]/2 ? wp
                                                          : wp - C->p[k]);
            C->itw[k][i] = (double)(slong)(wip <= C->p[k]/2 ? wip
                                                            : wip - C->p[k]);
            wp = nmod_mul(wp, w, C->mod[k]);
            wip = nmod_mul(wip, winv, C->mod[k]);
        }
    }
    fmpz_init(C->Pbig);
    fmpz_init(C->Phalf);
    fmpz_one(C->Pbig);
    for (slong k = 0; k < np; k++)
        fmpz_mul_ui(C->Pbig, C->Pbig, C->p[k]);
    fmpz_fdiv_q_2exp(C->Phalf, C->Pbig, 1);
    {
        /* bias B = 2^(2b + depth + 2) exceeds any coefficient
           magnitude; biased coefficients stay below the four-prime
           product. biaspoly = B * (2^(b m) - 1)/(2^b - 1). */
        slong Bexp = 2 * C->b + C->depth;
        slong nc2 = (slong) C->R->crts[np - 1].coeff_len;
        C->whi = (C->b * (m - 1)) / FLINT_BITS + nc2 + 2;
        (void) Bexp;
        for (slong k = 0; k < np; k++)
        {
            ulong r = nmod_pow_ui(2, (ulong) Bexp, C->mod[k]);
            if (np != 2)
                r = nmod_mul(r, nmod_inv(
                *crt_data_co_prime_red(C->R->crts + np - 1, k) % C->p[k],
                C->mod[k]), C->mod[k]);
            C->biasres[k] = r;
        }
        C->biaspoly = flint_calloc(C->whi + 1, sizeof(ulong));
        for (slong i = 0; i < m; i++)
        {
            slong bit = C->b * i + Bexp, l = bit / FLINT_BITS, o = bit % FLINT_BITS;
            ulong cy = mpn_add_1(C->biaspoly + l, C->biaspoly + l,
                                 C->whi + 1 - l, UWORD(1) << o);
            (void) cy;
        }
    }
    C->crt2 = flint_malloc(m * sizeof(slong));
}


/* For a caller-fixed ring 2^N + 1, choose the digit count m over
   the instantiated digit sizes dividing N (m = N/b a power of
   two), minimizing np * m * (depth + 4 + b/16). For power-of-two N
   this switches from (3, 64) to (6, 128) at half the transform
   length once three-prime capacity runs out, flattening the large-N
   cost staircase. */
static slong
_mulmod_2expp1_choose_m(mpn_ctx_struct * R, slong N)
{
    /* Chunk sizes come from N's divisor lattice: every b = N / m for
       a power-of-two m is a candidate, admissible when the digit
       path covers it (b <= 50) or a vectorized conversion is
       instantiated. The score carries measured per-prime-count
       weights on top of the butterfly-and-conversion estimate:
       raced over the admissible geometries of six
       Schoenhage-Strassen ring sizes from 47104 to 2752512, per-unit
       costs sit at 0.80-0.93 for np <= 7 with a small two-prime
       credit, while np = 8 cliffs to 1.08-1.28 (its seven-word crt
       chains); the weights below reproduce the measured winner at
       every size raced, including the two-prime b = 42 geometry
       (16% over the previous choice) and the escape from the np = 8
       tops of the m-octaves (up to 37%). */
    static const slong npw[9] = { 0, 0, 13, 14, 14, 14, 14, 14, 19 };
    slong m, bestm = 0;
    double bestc = 0.0;

    for (m = 16; m <= N / 20; m *= 2)
    {
        slong b, depth, np;
        double cost;

        if (N % m != 0)
            continue;
        b = N / m;
        if (b < 20 || b > 200)
            continue;
        depth = FLINT_BIT_COUNT(m) - 1;
        np = sd_fft_mpn_mulmod_2expp1_choose_np(R, b, depth);
        if (np < 0)
            continue;
        cost = (double) npw[np] * (double) np * (double) m
               * ((double) depth + 4.0 + (double) b / 16.0);
        if (bestm == 0 || cost < bestc)
        {
            bestc = cost;
            bestm = m;
        }
    }
    return bestm;
}

void
sd_fft_mpn_mulmod_2expp1_ctx_clear(sd_fft_mpn_mulmod_2expp1_ctx_struct * C)
{
    for (slong k = 0; k < C->np; k++)
    {
        flint_aligned_free(C->tw[k]);
        flint_aligned_free(C->itw[k]);
    }
    flint_free(C->crt2);
    if (C->biaspoly) flint_free(C->biaspoly);
    fmpz_clear(C->Pbig); fmpz_clear(C->Phalf);
}

/* digits of the N-bit residue x (n limbs), b-bit fields, as doubles
   ready for the vectorized twist. Vectorized in the style of the
   fft_small input pass, across output digits: byte-granular unaligned
   loads put each field within one 64-bit word (b <= 56 keeps the field
   plus the sub-byte offset inside it), a per-lane variable shift
   aligns it, a mask isolates it, and the limited conversion (values
   below 2^52) produces the doubles. The scalar tail also covers the
   final fields, whose byte loads could run past the n+1 limbs. */
static void
digits_of_d(double * dig, nn_srcptr x, slong n, slong b,
                        slong m)
{
    const unsigned char * bytes = (const unsigned char *) x;
    ulong mask = (UWORD(1) << b) - 1;
    slong i = 0;
#if defined(__AVX2__)
    if (b <= 56)
    {
        /* last i whose 8-byte load ends within the n+1 limbs */
        slong safe = (FLINT_BITS * (slong)(n + 1) - 71) / b - 3;
        __m256i vmask = _mm256_set1_epi64x((slong) mask);
        for (; i + 4 <= m && i + 4 <= safe; i += 4)
        {
            ulong b0 = (ulong)(b*(i+0)), b1 = (ulong)(b*(i+1)),
                  b2 = (ulong)(b*(i+2)), b3 = (ulong)(b*(i+3));
            ulong w0, w1, w2, w3;
            memcpy(&w0, bytes + b0/8, 8);
            memcpy(&w1, bytes + b1/8, 8);
            memcpy(&w2, bytes + b2/8, 8);
            memcpy(&w3, bytes + b3/8, 8);
            __m256i w = _mm256_setr_epi64x((slong) w0, (slong) w1,
                                           (slong) w2, (slong) w3);
            __m256i sh = _mm256_setr_epi64x((slong)(b0 % 8),
                (slong)(b1 % 8), (slong)(b2 % 8), (slong)(b3 % 8));
            w = _mm256_and_si256(_mm256_srlv_epi64(w, sh), vmask);
            vec4d_store(dig + i, vec4n_convert_limited_vec4d(w));
        }
    }
#endif
    for (; i < m; i++)
    {
        slong bit = i * b, l = bit / FLINT_BITS, o = bit % FLINT_BITS;
        ulong lo = x[l] >> o;
        if (o + b > FLINT_BITS && l + 1 < n)
            lo |= x[l+1] << (64 - o);
        dig[i] = (double)(slong)(lo & mask);
    }
}


/* z = x*y mod 2^N+1; x, y in [0, 2^N] (top limb 0 or the value 2^N),
   all n = N/64 limbs + 1 top limb; z likewise */
/* transform one operand (n limbs plus a 0/1 top limb) into np data
   blocks of dsz doubles: input pass, twist by w^i, forward transform.
   The result depends only on the operand and the context and may be
   cached across sd_fft_mpn_mulmod_2expp1_mul_cached calls. */
static void
_mulmod_2expp1_transform_op(const sd_fft_mpn_mulmod_2expp1_ctx_struct * C,
            double * F, nn_srcptr x,
            sd_fft_mpn_mulmod_2expp1_scratch_struct * S)
{
    slong m = C->m, b = C->b, n = C->N / 64, np = C->np;
    slong dsz = C->dsz;
    ulong bits = (ulong) b;

    if (b > 50)
    {
        ulong stop_easy = n_min((ulong) m, (FLINT_BITS*(ulong)n - 33)/bits);
        stop_easy &= -(ulong) (bits % 8 == 0 ? 4 : bits % 4 == 0 ? 8 : 16);
        ulong stop_hard = n_min((ulong) m,
                                (FLINT_BITS*(ulong)n + bits - 1)/bits);
        const vec4d * tp = C->R->vec_two_pow_tab[(np + 3)/4 - 1];
        to_ffts_func tf = _mpn_ctx_to_ffts_func((ulong) np, bits);
        FLINT_ASSERT(tf != NULL);
        tf(C->R->ffts, F, dsz, x, n, (ulong) m, tp,
           0, stop_easy, stop_easy, stop_hard);
    }
    else
    {
        digits_of_d(S->digs, x, n, b, m);
        S->digs[0] -= (double)(x[n] != 0);
    }

    for (slong k = 0; k < np; k++)
    {
        sd_fft_ctx_struct * Q = C->R->ffts + k;
        double * X = F + k * dsz;
        vec1d pd = Q->p, pinv = Q->pinv;
        vec8d n8 = vec8d_set_d(pd), ninv8 = vec8d_set_d(pinv);

        if (b > 50)
        {
            if (x[n]) X[0] -= 1.0;
            for (slong i = 0; i + 8 <= m; i += 8)
            {
                vec8d tw8 = vec8d_load(C->tw[k] + i);
                vec8d_store(X + i,
                    vec8d_mulmod(vec8d_load(X + i), tw8, n8, ninv8));
            }
        }
        else
        {
            const double * dx = S->digs;
            for (slong i = 0; i + 8 <= m; i += 8)
            {
                vec8d tw8 = vec8d_load(C->tw[k] + i);
                vec8d_store(X + i,
                    vec8d_mulmod(vec8d_load(dx + i), tw8, n8, ninv8));
            }
        }
        sd_fft_trunc(Q, X, C->depth, m, m);
    }
}

slong
sd_fft_mpn_mulmod_2expp1_transformed_size(
            const sd_fft_mpn_mulmod_2expp1_ctx_struct * C)
{
    return C->np * C->dsz;
}

void
sd_fft_mpn_mulmod_2expp1_transform(sd_fft_mpn_mulmod_2expp1_ctx_struct * C,
            double * F, nn_srcptr y,
            sd_fft_mpn_mulmod_2expp1_scratch_struct * S)
{
    _mulmod_2expp1_transform_op(C, F, y, S);
}

void
sd_fft_mpn_mulmod_2expp1_mul_cached(sd_fft_mpn_mulmod_2expp1_ctx_struct * C,
            nn_ptr z, nn_srcptr x, const double * F,
            sd_fft_mpn_mulmod_2expp1_scratch_struct * S)
{
    slong m = C->m, n = C->N / 64, np = C->np;
    slong dsz = C->dsz;
    double * abuf = S->buf;

    _mulmod_2expp1_transform_op(C, abuf, x, S);

    for (slong k = 0; k < np; k++)
    {
        sd_fft_ctx_struct * Q = C->R->ffts + k;
        double * X = abuf + k * dsz;
        const double * Y = F + k * dsz;
        vec1d pd = Q->p, pinv = Q->pinv;
        {
            vec8d n8 = vec8d_set_d(pd), ninv8 = vec8d_set_d(pinv);
            for (slong i = 0; i + 8 <= m; i += 8)
            {
                vec8d t = vec8d_reduce_to_pm1n(vec8d_load(X + i),
                                               n8, ninv8);
                vec8d_store(X + i, vec8d_mulmod(t,
                    vec8d_load(Y + i), n8, ninv8));
            }
        }
        sd_ifft_trunc(Q, X, C->depth, m);
        {
            /* untwist, then bias into nonnegative residues for the
               unsigned from_ffts reconstruction */
            vec8d n8 = vec8d_set_d(pd), ninv8 = vec8d_set_d(pinv);
            vec1d br = (double)(slong)(C->biasres[k] <= C->p[k]/2
                    ? C->biasres[k]
                    : C->biasres[k] - C->p[k]);
            vec8d b8 = vec8d_set_d(br);
            for (slong i = 0; i + 8 <= m; i += 8)
            {
                vec8d t = vec8d_mulmod(vec8d_load(X + i),
                    vec8d_load(C->itw[k] + i), n8, ninv8);
                t = vec8d_add(t, b8);
                t = vec8d_reduce_to_pm1n(t, n8, ninv8);
                vec8d_store(X + i, t);
            }
        }
    }

    {
        nn_ptr w = S->w;
        slong ncoef = (slong) C->R->crts[np - 1].coeff_len;
        ulong E0 = 0;
        ulong E1 = (((ulong) C->whi >= (ulong) ncoef + 1)
                    ? (ulong) C->whi - ((ulong) ncoef + 1) : UWORD(0))
                   * FLINT_BITS / (ulong) C->b;
        E1 &= -(ulong) BLK_SZ;
        if (E1 > (ulong) m) E1 = (ulong) m & -(ulong) BLK_SZ;
        /* from_ffts zeroes its own window */
        if (np == 2)
        {
            /* no vectorized from_ffts instance exists for two primes;
               the shared Garner recomposition serves, on the plain
               residues of the cofactor-free untwist */
            _fft_small_np2_crt_recompose(C->R, w, 0, (ulong) C->whi,
                abuf, abuf + dsz, (ulong) m, (ulong) m, (ulong) C->b,
                UWORD_MAX);
            (void) E0; (void) E1;
        }
        else
        {
            /* CRT recomposition via mpn_mul.c instances */
            _mpn_ctx_from_ffts_func((ulong) np)(w, 0, (ulong) C->whi, 0, (ulong) m,
                C->R->ffts, abuf, dsz, C->R->crts, (ulong) C->b,
                E0, E1, NULL, NULL);
        }
        mpn_sub_n(w, w, C->biaspoly, C->whi);
        /* fold: same shape as the small path; w[n..] is the
           signed excess */
        {
            /* the signed excess spans b + depth + 2 bits above 2^N:
               take the full window, sign-extending past whi */
            slong elen = (C->b + C->depth + 2) / FLINT_BITS + 2;
            ulong ex[8];
            ulong sign;
            slong t;

            FLINT_ASSERT(elen <= 8);
            for (t = 0; t < elen; t++)
                ex[t] = (n + t < C->whi) ? w[n + t]
                        : (ulong)((slong) ex[t - 1] >> (FLINT_BITS - 1));
            sign = (slong) ex[elen - 1] < 0;
            if (sign)
            {
                for (t = 0; t < elen; t++)
                    ex[t] = ~ex[t];
                mpn_add_1(ex, ex, elen, 1);
            }
            z[n] = 0;
            if (!sign)
            {
                if (mpn_sub(z, w, n, ex, elen))
                {
                    if (mpn_add_1(z, z, n, 1))
                    {
                        flint_mpn_zero(z, n);
                        z[n] = 1;
                    }
                }
            }
            else
            {
                if (mpn_add(z, w, n, ex, elen))
                {
                    if (mpn_sub_1(z, z, n, 1))
                    {
                        flint_mpn_zero(z, n);
                        z[n] = 1;
                    }
                }
            }
            /* result is fully reduced: the unique representative in
               [0, 2^N], top limb set only for 2^N itself */
            FLINT_ASSERT(z[n] == 0
                         || (z[n] == 1 && flint_mpn_zero_p(z, n)));
        }
    }
}

void
sd_fft_mpn_mulmod_2expp1(sd_fft_mpn_mulmod_2expp1_ctx_struct * C, nn_ptr z, nn_srcptr x,
            nn_srcptr y, sd_fft_mpn_mulmod_2expp1_scratch_struct * S)
{
    double * bbuf = S->buf + C->np * C->dsz;

    _mulmod_2expp1_transform_op(C, bbuf, y, S);
    sd_fft_mpn_mulmod_2expp1_mul_cached(C, z, x, bbuf, S);
}

/* per-thread mutable state for sd_fft_mpn_mulmod_2expp1; the context itself is
   read-only during multiplication and may be shared */
void
sd_fft_mpn_mulmod_2expp1_scratch_init(sd_fft_mpn_mulmod_2expp1_scratch_struct * S,
                        const sd_fft_mpn_mulmod_2expp1_ctx_struct * C)
{
    S->buf = flint_aligned_alloc(FLINT_FFT_SMALL_ALIGNMENT,
                 n_round_up(2 * C->np * C->dsz * sizeof(double),
                            FLINT_FFT_SMALL_ALIGNMENT));
    S->w = flint_malloc((C->N/FLINT_BITS + 6) * sizeof(ulong));
    S->digs = (C->b <= 50)
        ? flint_aligned_alloc(64, 2 * C->m * sizeof(double)) : NULL;
}

void
sd_fft_mpn_mulmod_2expp1_scratch_clear(sd_fft_mpn_mulmod_2expp1_scratch_struct * S)
{
    flint_aligned_free(S->buf);
    flint_free(S->w);
    if (S->digs != NULL)
        flint_aligned_free(S->digs);
}
