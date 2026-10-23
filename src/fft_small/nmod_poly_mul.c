/*
    Copyright (C) 2022 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_pool.h"
#include "thread_support.h"
#include "mpn_extras.h"
#include "ulong_extras.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "crt_helpers.h"
#include "fft_small.h"

static void _mod_red(
    double* abuf, ulong atrunc,
    const ulong* a, ulong an,
    const sd_fft_ctx_struct* fft,
    nmod_t mod)
{
    double* aI;
    ulong j;

    FLINT_ASSERT(atrunc < an);
    /* flat over the transform: short lengths are simply fewer
       iterations, the addressing is unchanged */
    FLINT_ASSERT(atrunc % 8 == 0);

#if 1

    ulong tt = an%atrunc;

#define UNROLL 8

    {
        aI = abuf;

        vec8n nn = vec8n_set_n(mod.n);
        vec8d n = vec8d_set_d(fft->p);
        vec8d ninv = vec8d_set_d(fft->pinv);

        for (j = 0; j < atrunc; j += UNROLL)
        {
            if (j+UNROLL <= tt || j >= tt)
            {
                ulong k = j;
                FLINT_ASSERT(k+UNROLL-1 < an);
                vec8n t = vec8n_load_unaligned(a + k);

                if (mod.norm == 0)
                    for (k += atrunc; k < an; k += atrunc)
                        t = vec8n_addmod(t, vec8n_load_unaligned(a + k), nn);
                else
                    for (k += atrunc; k < an; k += atrunc)
                        t = vec8n_addmod_limited(t, vec8n_load_unaligned(a + k), nn);


                vec8d tlo = vec8n_convert_limited_vec8d(vec8n_bit_and(t, vec8n_set_n(n_pow2(32)-1)));
                vec8d thi = vec8n_convert_limited_vec8d(vec8n_bit_shift_right_32(t));
                vec8d_store_aligned(aI + j, vec8d_add(tlo, vec8d_mulmod(thi, vec8d_set_d(n_pow2(32)), n, ninv)));
            }
            else
            {
                for (ulong l = 0; l < UNROLL; l++)
                {
                    ulong k = j+l;
                    ulong c = a[k];
                    for (k += atrunc; k < an; k += atrunc)
                        c = nmod_add(c, a[k], mod);

                    aI[j+l] = (slong)(nmod_set_ui(c, fft->mod));
                }
            }
        }
    }

#else
    // wrong way!!
    if (modn <= fft->mod.n)
    {
        /* first pass fill in */

        for (i = 0; i < atrunc; i += BLK_SZ)
        {
            aI = sd_fft_ctx_blk_index(abuf, i/BLK_SZ);
            for (j = 0; j < BLK_SZ; j += 8)
            {
                vec8n t = vec8n_load_unaligned(a + i + j);
                vec8d_store_aligned(aI + j, vec8n_convert_limited_vec8d(t));
            }
        }

        vec8d n = vec8d_set_d(fft->p);

        /* second pass add to existing */

        for (i = atrunc; i + BLK_SZ <= an; i += BLK_SZ)
        {
            aI = sd_fft_ctx_blk_index(abuf, (i%atrunc)/BLK_SZ);
            for (j = 0; j < BLK_SZ; j += 8)
            {
                vec8n t = vec8n_load_unaligned(a + i + j);
                vec8d s = vec8n_convert_limited_vec8d(t);
                s = vec8d_add(s, vec8d_load_aligned(aI + j));
                s = vec8d_reduce_2n_to_n(s, n);
                vec8d_store_aligned(aI + j, s);
            }
        }

        aI = sd_fft_ctx_blk_index(abuf, (i%atrunc)/BLK_SZ);
        for (j = 0; j < an - i; j++)
            aI[j] = vec1d_reduce_2n_to_n(aI[j] + (slong)a[i + j], fft->p);
    }
    else
    {
        vec8d n = vec8d_set_d(fft->p);
        vec8d ninv = vec8d_set_d(fft->pinv);
        for (i = 0; i < atrunc; i += BLK_SZ)
        {
            aI = sd_fft_ctx_blk_index(abuf, i/BLK_SZ);
            for (j = 0; j < BLK_SZ; j += 8)
            {
                vec8n t = vec8n_load_unaligned(a + i + j);
                vec8d tlo = vec8n_convert_limited_vec8d(vec8n_bit_and(t, vec8n_set_n(n_pow2(32)-1)));
                vec8d thi = vec8n_convert_limited_vec8d(vec8n_bit_shift_right_32(t));
                vec8d_store_aligned(aI + j, vec8d_add(tlo, vec8d_mulmod(thi, vec8d_set_d(n_pow2(32)), n, ninv)));
            }
        }

        for (i = atrunc; i + BLK_SZ <= an; i += BLK_SZ)
        {
            aI = sd_fft_ctx_blk_index(abuf, (i%atrunc)/BLK_SZ);
            for (j = 0; j < BLK_SZ; j += 8)
            {
                vec8n t = vec8n_load_unaligned(a + i + j);
                vec8d tlo = vec8n_convert_limited_vec8d(vec8n_bit_and(t, vec8n_set_n(n_pow2(32)-1)));
                vec8d thi = vec8n_convert_limited_vec8d(vec8n_bit_shift_right_32(t));
                vec8d s = vec8d_add(tlo, vec8d_mulmod(thi, vec8d_set_d(n_pow2(32)), n, ninv));
                s = vec8d_add(s, vec8d_load_aligned(aI + j));
                s = vec8d_reduce_2n_to_n(s, n);
                vec8d_store_aligned(aI + j, s);
            }
        }

        aI = sd_fft_ctx_blk_index(abuf, (i%atrunc)/BLK_SZ);
        for (j = 0; j < an - i; j++)
            aI[j] = vec1d_reduce_2n_to_n(aI[j] + (slong)nmod_set_ui(a[i+j], fft->mod), fft->p);
    }
#endif
}

static void _mod(
    double* abuf, ulong atrunc,
    const ulong* a, ulong an,
    const sd_fft_ctx_struct* fft,
    nmod_t mod)
{
    double* aI;
    ulong i, j;

    if (atrunc < an)
    {
        _mod_red(abuf, atrunc, a, an, fft, mod);
        return;
    }

    if (mod.n <= fft->mod.n)
    {
        for (i = 0; i + BLK_SZ <= an; i += BLK_SZ)
        {
            aI = sd_fft_ctx_blk_index(abuf, i/BLK_SZ);
            for (j = 0; j < BLK_SZ; j += 8)
            {
                vec8n t = vec8n_load_unaligned(a + i + j);
FLINT_ASSERT(i+j < atrunc);
                vec8d_store_aligned(aI + j, vec8n_convert_limited_vec8d(t));
            }
        }

        aI = sd_fft_ctx_blk_index(abuf, i/BLK_SZ);
        for (j = 0; j < an - i; j++)
            aI[j] = (slong)a[i + j];
    }
    else
    {
        vec8d n = vec8d_set_d(fft->p);
        vec8d ninv = vec8d_set_d(fft->pinv);
        for (i = 0; i + BLK_SZ <= an; i += BLK_SZ)
        {
            aI = sd_fft_ctx_blk_index(abuf, i/BLK_SZ);
#if 1
            for (j = 0; j < BLK_SZ; j += 8)
            {
                vec8n t = vec8n_load_unaligned(a + i + j);
                vec8d tlo = vec8n_convert_limited_vec8d(vec8n_bit_and(t, vec8n_set_n(n_pow2(32)-1)));
                vec8d thi = vec8n_convert_limited_vec8d(vec8n_bit_shift_right_32(t));
                vec8d_store_aligned(aI + j, vec8d_add(tlo, vec8d_mulmod(thi, vec8d_set_d(n_pow2(32)), n, ninv)));
            }
#else
            for (j = 0; j < BLK_SZ; j += 1)
            {
                aI[j] = (slong)nmod_set_ui(a[i+j], fft->mod);
            }
#endif
        }

        aI = sd_fft_ctx_blk_index(abuf, i/BLK_SZ);
        for (j = 0; j < an - i; j++)
            aI[j] = (slong)nmod_set_ui(a[i+j], fft->mod);
    }

    for (i = an; i < atrunc; i++)
        abuf[i] = 0;
}



/* _crt_1, _crt_2, _crt_3 are special-cased below; the generic template
   handles four or more primes, with a fast path while the reconstructed
   values stay below 2^192 */
#define CRT_FN(NP) CAT3(_crt, NP, fn)
#define CRT_Z_TYPE ulong*
#define CRT_HEAD nmod_t mod = *(const nmod_t*) params;
#define CRT_EMIT(zi, r, NP, N, M) \
    if (N == 3 || (N == 4 && r[3] == 0)) \
    { \
        NMOD_RED3(z[zi], r[2], r[1], r[0], mod); \
    } \
    else \
    { \
        z[zi] = mpn_mod_1(r, N, mod.n); \
    }
#define CRT_TAIL
CRT_DEFINE(4, 4, 3)  /* 200 bits */
CRT_DEFINE(5, 4, 4)  /* 250 bits */
CRT_DEFINE(6, 5, 4)  /* 300 bits */
CRT_DEFINE(7, 6, 5)  /* 350 bits */
CRT_DEFINE(8, 7, 6)  /* 400 bits */
#undef CRT_FN
#undef CRT_Z_TYPE
#undef CRT_HEAD
#undef CRT_EMIT
#undef CRT_TAIL
#undef CRT_DEFINE


FLINT_FORCE_INLINE
void mullo_2x1(ulong * r1, ulong * r0, ulong a1, ulong a0, ulong b0)
{
    ulong t0, t1;
    umul_ppmm(t1, t0, a0, b0);
    t1 += a1 * b0;
    *r0 = t0;
    *r1 = t1;
}

static void _crt_1(
    ulong* z, ulong zl, ulong zi_start, ulong zi_stop,
    sd_fft_ctx_struct* Rffts, double* d, ulong dstride,
    crt_data_struct* FLINT_UNUSED(Rcrts),
    ulong FLINT_UNUSED(min_an_bn),
    nmod_t mod)
{
    ulong i, j, jstart, jstop;
    ulong Xs[BLK_SZ*1];
    int have_fft_prime;
    double n = 0.0, ninv = 0.0;

    /* Reconstructing from a single prime; the modulus is either an FFT
       prime itself or small enough to be a valid FFT prime. */
    FLINT_ASSERT(mod.n <= (UWORD(1) << 50));

    /* This applies in the `direct_fft` branch, or when `mod.n` is already
       in the context. In the latter case, `s2worker_func` passes
       `Rffts = ffts + offset`, making the matched prime `Rffts[0]`.  */
    have_fft_prime = (mod.n == Rffts[0].mod.n);
    if (!have_fft_prime)
    {
        n = mod.n;
        ninv = 1.0 / n;
    }

    for (i = n_round_down(zi_start, BLK_SZ); i < zi_stop; i += BLK_SZ)
    {
        jstart = (i < zi_start) ? zi_start - i : 0; \
        jstop = FLINT_MIN(BLK_SZ, zi_stop - i);

        /* Can write block directly to output */
        if (jstart == 0 && jstop == BLK_SZ)
        {
            if (have_fft_prime)
                _convert_block(z + i - zl, Rffts, d, dstride, 1, i/BLK_SZ);
            else
                _convert_block_1_mod(z + i - zl, Rffts, d, dstride, i/BLK_SZ, n, ninv);
        }
        else /* Write to temporary buffer, then copy truncated part of block */
        {
            if (have_fft_prime)
                _convert_block(Xs, Rffts, d, dstride, 1, i/BLK_SZ);
            else
                _convert_block_1_mod(Xs, Rffts, d, dstride, i/BLK_SZ, n, ninv);

            for (j = jstart; j < jstop; j += 1)
                z[i+j-zl] = Xs[j];
        }
    }
}

static void _crt_2(
    ulong* z, ulong zl, ulong zi_start, ulong zi_stop,
    sd_fft_ctx_struct* Rffts, double* d, ulong dstride,
    crt_data_struct* Rcrts, ulong min_an_bn,
    nmod_t mod)
{
    ulong np = 2;
    ulong n = 2;
    ulong m = 1;

    FLINT_ASSERT(n == Rcrts[np-1].coeff_len);

    if (n == m + 1)
    {
        for (ulong l = 0; l < np; l++) {
            FLINT_ASSERT(crt_data_co_prime(Rcrts + np - 1, l)[m] == 0);
        }
    }
    else
    {
        FLINT_ASSERT(n == m);
    }

    ulong Xs[BLK_SZ*2];

    /* If the output fits in 100 bits then certainly the modulus
       should have <= 63 bits. */
    FLINT_ASSERT(mod.n < (UWORD(1) << (FLINT_BITS - 1)));

    /* The output coefficients are bounded by min(an,bn) * (n-1)^2. If
       the high words are smaller than n (i.e. already reduced mod n), we
       can use NMOD_RED2 (or NMOD_RED2_NONFULLWORD), otherwise NMOD2_RED2.
       We could do something faster than NMOD2_RED2, but in practice the bound
       is rarely exceeded (we need products longer than 10^10 or so), so just
       keep it simple and fall back on NMOD2_RED2 when that occurs. */
    int high_reduced = 0;

    {
        ulong hi, lo;
        umul_ppmm(hi, lo, mod.n - 1, mod.n - 1);
        mullo_2x1(&hi, &lo, hi, lo, min_an_bn);
        if (hi < mod.n)
            high_reduced = 1;
    }

#if defined(__SIZEOF_INT128__)

    /* The two-limb CRT is slightly faster using __uint128_t than using
       the generic macros. */
    __uint128_t crt_M = ((__uint128_t) crt_data_prod_primes(Rcrts + np - 1)[0]) |
                        (((__uint128_t) crt_data_prod_primes(Rcrts + np - 1)[1]) << 64);
    ulong c0 = _crt_data_co_prime(Rcrts + np - 1, 0, 2)[0];
    ulong c1 = _crt_data_co_prime(Rcrts + np - 1, 1, 2)[0];

#define DO_CRT \
    ulong r[2]; \
    __uint128_t rr; \
    rr = ((__uint128_t) c0) * ((__uint128_t) (Xs[0*BLK_SZ + j])) + \
         ((__uint128_t) c1) * ((__uint128_t) (Xs[1*BLK_SZ + j])); \
    if (rr >= crt_M) \
        rr -= crt_M; \
    r[0] = rr; \
    r[1] = rr >> 64;

#else

#define DO_CRT \
    ulong r[2]; \
    ulong t[2]; \
    ulong l = 0; \
    CAT3(_big_mul, 2, 1)(r, t, _crt_data_co_prime(Rcrts + np - 1, l, n), Xs[l*BLK_SZ + j]); \
    CRT_UNROLL_PRIMES \
    for (l++; l < np; l++) \
        CAT3(_big_addmul, 2, 1)(r, t, _crt_data_co_prime(Rcrts + np - 1, l, n), Xs[l*BLK_SZ + j]); \
    CAT(_reduce_big_sum, 2)(r, t, crt_data_prod_primes(Rcrts + np - 1)); \

#endif

    for (ulong i = n_round_down(zi_start, BLK_SZ); i < zi_stop; i += BLK_SZ)
    {
        _convert_block(Xs, Rffts, d, dstride, np, i/BLK_SZ);

        ulong jstart = (i < zi_start) ? zi_start - i : 0;
        ulong jstop = FLINT_MIN(BLK_SZ, zi_stop - i);

        if (high_reduced)
        {
            for (ulong j = jstart; j < jstop; j += 1)
            {
                DO_CRT
                FLINT_ASSERT(r[1] < mod.n);
                NMOD_RED2_NONFULLWORD(z[i+j-zl], r[1], r[0], mod);
            }
        }
        else
        {
            for (ulong j = jstart; j < jstop; j += 1)
            {
                DO_CRT
                NMOD2_RED2(z[i+j-zl], r[1], r[0], mod);
            }
        }
    }
#undef DO_CRT
}

/* Ad hoc modular reduction to do better than NMOD2_RED2 and NMOD_RED3.
   This should eventually be moved out and used elsewhere. */

/* Precompute floor(2^(2*FLINT_BITS) / n) for Barrett-style modular reduction.
   This could be optimized (not necessary here, but relevant for other
   applications). */
static void
n_ll_rem_l_precomp(ulong * qhi, ulong * qlo, ulong n)
{
    ulong q[4];
    ulong a[4];
    a[0] = 0;
    a[1] = 0;
    a[2] = 1;
    mpn_divrem_1(q, 0, a, 3, n);
    *qlo = q[0];
    *qhi = q[1];
}

/* 2 -> 1 limb mod, n < 2^(FLINT_BITS-1), using linear combination + Barrett */
FLINT_FORCE_INLINE ulong
n_ll_rem_l_nonfullword(ulong xhi, ulong xlo, ulong n, ulong qhi, ulong qlo)
{
    ulong c2, c1, c0;
    FLINT_MPN_MUL_3P2X2(c2, c1, c0, qhi, qlo, xhi, xlo);
    (void) c1;
    (void) c0;
    xlo -= c2 * n;
    if (xlo >= n)
        xlo -= n;
    return xlo;
}

/* 3 -> 1 limb mod, n < 2^(FLINT_BITS-1), using linear combination + Barrett */
FLINT_FORCE_INLINE ulong
n_lll_rem_l_nonfullword(ulong y2, ulong y1, ulong y0, ulong n, ulong qhi, ulong qlo, ulong alpha2, ulong alpha1)
{
    ulong c2, c1, c0, t1, t0;
    ulong xhi, xlo;

    FLINT_ASSERT(n < (UWORD(1) << (FLINT_BITS - 1)));

    umul_ppmm(t1, t0, y2, alpha2);
    umul_ppmm(c1, c0, y1, alpha1);
    add_ssaaaa(xhi, xlo, t1, t0, c1, c0);
    add_ssaaaa(xhi, xlo, xhi, xlo, 0, y0);

    FLINT_MPN_MUL_3P2X2(c2, c1, c0, qhi, qlo, xhi, xlo);
    (void) c1;
    (void) c0;
    xlo -= c2 * n;
    if (xlo >= n)
        xlo -= n;

    return xlo;
}

/* 3 -> 1 limb mod, n >= 2^(FLINT_BITS-1), using linear combination + Granlund-Möller */
FLINT_FORCE_INLINE ulong
n_lll_rem_l_fullword(ulong y2, ulong y1, ulong y0, nmod_t mod, ulong alpha2, ulong alpha1)
{
    ulong c1, c0, t1, t0;
    ulong xhi, xlo;

    FLINT_ASSERT(mod.n >= (UWORD(1) << (FLINT_BITS - 1)));

    umul_ppmm(t1, t0, y2, alpha2);
    umul_ppmm(c1, c0, y1, alpha1);
    add_ssaaaa(xhi, xlo, t1, t0, c1, c0);
    add_ssaaaa(xhi, xlo, xhi, xlo, 0, y0);

    if (xhi >= mod.n) xhi -= mod.n;
    NMOD_RED2_FULLWORD(xlo, xhi, xlo, mod);

    return xlo;
}

/* For moduli just larger than 2^63, the conditional subtraction can be shown
   to be redundant. */
FLINT_FORCE_INLINE ulong
n_lll_rem_l_fullword_limited(ulong y2, ulong y1, ulong y0, nmod_t mod, ulong alpha2, ulong alpha1)
{
    ulong c1, c0, t1, t0;
    ulong xhi, xlo;

    FLINT_ASSERT(mod.n >= (UWORD(1) << (FLINT_BITS - 1)));
    FLINT_ASSERT(mod.n < (UWORD(1) << (FLINT_BITS - 1)) + (UWORD(1) << (FLINT_BITS / 2 - 2)));

    umul_ppmm(t1, t0, y2, alpha2);
    umul_ppmm(c1, c0, y1, alpha1);
    add_ssaaaa(xhi, xlo, t1, t0, c1, c0);
    add_ssaaaa(xhi, xlo, xhi, xlo, 0, y0);

    NMOD_RED2_FULLWORD(xlo, xhi, xlo, mod);

    return xlo;
}

static void _crt_3(
    ulong* z, ulong zl, ulong zi_start, ulong zi_stop,
    sd_fft_ctx_struct* Rffts, double* d, ulong dstride,
    crt_data_struct* Rcrts, ulong min_an_bn,
    nmod_t mod)
{
    ulong np = 3;
    ulong n = 3;

#if FLINT_WANT_ASSERT
    ulong m = 2;
    FLINT_ASSERT(n == Rcrts[np-1].coeff_len);
    for (ulong l = 0; l < np; l++)
        FLINT_ASSERT(crt_data_co_prime(Rcrts + np - 1, l)[m] == 0);
#endif

    ulong Xs[BLK_SZ*3];

    ulong qhi = 0, qlo = 0;
    ulong alpha1 = 0, alpha2 = 0;
    ulong hi, lo, u2, u1, u0;

    /* Coefficients before modular reduction can have 2 or 3 limbs. */
    int two_limbs = 0;
    /* High word of 2-limb input is reduced */
    int high_reduced = 0;
    /* Modulus has MSB set. */
    int mod_fullword = (mod.norm == 0);

    /* Bound coefficients before modular reduction. */
    umul_ppmm(hi, lo, mod.n - 1, mod.n - 1);
    FLINT_MPN_MUL_2X1(u2, u1, u0, hi, lo, min_an_bn);
    (void) u0;
    two_limbs = (u2 == 0);

    /* For moduli close to 2^63, we can avoid the high word reduction
       before NMOD_RED2. */
    int fullword_limited = 0;

    /* Branch 1 and 2 of mod code */
    if (two_limbs && !mod_fullword)
    {
        high_reduced = (u1 < mod.n);

        /* Extra precomputation for branch 2 */
        if (!high_reduced)
            n_ll_rem_l_precomp(&qhi, &qlo, mod.n);
    }
    else
    {
        /* Extra precomputation for branch 3 and 4 of mod code */
        alpha1 = nmod_set_ui(UWORD(1) << (FLINT_BITS / 2), mod);
        alpha1 = nmod_mul(alpha1, alpha1, mod);
        alpha2 = nmod_mul(alpha1, alpha1, mod);

        /* Extra precomputation for branch 3 */
        if (!mod_fullword)
            n_ll_rem_l_precomp(&qhi, &qlo, mod.n);

        if (mod_fullword && (mod.n < ((UWORD(1) << 63) + (UWORD(1) << 30))))
            fullword_limited = 1;
    }

#define DO_CRT \
    ulong r[3]; \
    ulong t[3]; \
    ulong l = 0; \
    CAT3(_big_mul, 3, 2)(r, t, _crt_data_co_prime(Rcrts + np - 1, l, n), Xs[l*BLK_SZ + j]); \
    CRT_UNROLL_PRIMES \
    for (l++; l < np; l++) \
        CAT3(_big_addmul, 3, 2)(r, t, _crt_data_co_prime(Rcrts + np - 1, l, n), Xs[l*BLK_SZ + j]); \
    CAT(_reduce_big_sum, 3)(r, t, crt_data_prod_primes(Rcrts + np - 1));

    for (ulong i = n_round_down(zi_start, BLK_SZ); i < zi_stop; i += BLK_SZ)
    {
        _convert_block(Xs, Rffts, d, dstride, np, i/BLK_SZ);

        ulong jstart = (i < zi_start) ? zi_start - i : 0;
        ulong jstop = FLINT_MIN(BLK_SZ, zi_stop - i);

        if (two_limbs && !mod_fullword)
        {
            if (high_reduced)
            {
                for (ulong j = jstart; j < jstop; j += 1)
                {
                    DO_CRT
                    FLINT_ASSERT(r[2] == 0);
                    NMOD_RED2_NONFULLWORD(z[i+j-zl], r[1], r[0], mod);
                }
            }
            else
            {
                for (ulong j = jstart; j < jstop; j += 1)
                {
                    DO_CRT
                    FLINT_ASSERT(r[2] == 0);
                    z[i+j-zl] = n_ll_rem_l_nonfullword(r[1], r[0], mod.n, qhi, qlo);
                }
            }
        }
        else
        {
            if (!mod_fullword)
            {
                for (ulong j = jstart; j < jstop; j += 1)
                {
                    DO_CRT
                    z[i+j-zl] = n_lll_rem_l_nonfullword(r[2], r[1], r[0], mod.n, qhi, qlo, alpha2, alpha1);
                }
            }
            else if (fullword_limited)
            {
                for (ulong j = jstart; j < jstop; j += 1)
                {
                    DO_CRT
                    z[i+j-zl] = n_lll_rem_l_fullword_limited(r[2], r[1], r[0], mod, alpha2, alpha1);
                }
            }
            else
            {
                for (ulong j = jstart; j < jstop; j += 1)
                {
                    DO_CRT
                    z[i+j-zl] = n_lll_rem_l_fullword(r[2], r[1], r[0], mod, alpha2, alpha1);
                }
            }
        }
    }
#undef DO_CRT
}



/* ------------------------------------------------------------------------ */
/* glue for the generic convolution engine                                  */
/* ------------------------------------------------------------------------ */

static void _mod_fn(double* xbuf, ulong xtrunc,
    const void* x, ulong xn, slong FLINT_UNUSED(xaux),
    const sd_fft_ctx_struct* fft, const void* params)
{
    _mod(xbuf, xtrunc, (const ulong*) x, xn, fft, *(const nmod_t*) params);
}

#define DEFINE_IT(NP) \
static void CAT3(_crt, NP, fn)(void* z, ulong zl, \
    ulong zi_start, ulong zi_stop, \
    sd_fft_ctx_struct* Rffts, double* d, ulong dstride, \
    crt_data_struct* Rcrts, ulong min_an_bn, ulong* FLINT_UNUSED(local), \
    const void* params) \
{ \
    CAT(_crt, NP)((ulong*) z, zl, zi_start, zi_stop, Rffts, d, dstride, \
                  Rcrts, min_an_bn, *(const nmod_t*) params); \
}
DEFINE_IT(1)
DEFINE_IT(2)
DEFINE_IT(3)
#undef DEFINE_IT

static fft_small_crt_func _nmod_crt_fn(ulong np)
{
    return np == 1 ? _crt_1_fn :
           np == 2 ? _crt_2_fn :
           np == 3 ? _crt_3_fn :
           np == 4 ? _crt_4_fn :
           np == 5 ? _crt_5_fn :
           np == 6 ? _crt_6_fn :
           np == 7 ? _crt_7_fn :
                     _crt_8_fn;
}

/* todo: even with np == 1, we may want threads for the conversion to fft */
static ulong _nmod_conv_s1_threads(ulong np, ulong tune_n)
{
    return tune_n < 1500 ? 1 : np;
}

static int _nmod_conv_s1_b_worker(ulong FLINT_UNUSED(np), ulong tune_n,
                                  int squaring)
{
    return tune_n > 5000 && !squaring;
}

static int _nmod_conv_s2_rethread(ulong np, ulong zn)
{
    return np*zn > 10000;
}

static const fft_small_conv_tuning _nmod_conv_tuning =
    { _nmod_conv_s1_threads, _nmod_conv_s1_b_worker, _nmod_conv_s2_rethread };

/* the cyclic and precomputed-operand drivers historically do not
   re-request threads for stage 2 */
static const fft_small_conv_tuning _nmod_conv_tuning_no_rethread =
    { _nmod_conv_s1_threads, _nmod_conv_s1_b_worker, NULL };

typedef struct {
    fft_small_op_struct* X;
    const ulong* a;
    ulong an;
    ulong itrunc;
    int transposed;
    nmod_t mod;
    const fft_small_plan_struct* P;
    ulong start_pi;
    ulong stop_pi;
} _op_fft_worker_struct;

static void _op_fft_worker_func(void* varg)
{
    _op_fft_worker_struct* W = (_op_fft_worker_struct*) varg;
    const fft_small_plan_struct* P = W->P;

    for (ulong i = W->start_pi; i < W->stop_pi; i++)
    {
        sd_fft_ctx_struct* Q = P->ffts + P->offset + i;
        double* xbuf = W->X->data + P->stride*i;

        if (W->transposed)
        {
            /* the conversion zero-fills up to ztrunc, which
               sd_ifft_trunc_tr reads in full */
            _mod(xbuf, P->ztrunc, W->a, W->an, Q, W->mod);
            sd_ifft_trunc_tr(Q, xbuf, P->depth, P->ztrunc);
        }
        else
        {
            _mod(xbuf, W->itrunc, W->a, W->an, Q, W->mod);
            sd_fft_trunc(Q, xbuf, P->depth, W->itrunc, P->ztrunc);
        }
    }
}

/* convert (a, an) and apply the forward transform truncated to itrunc, or,
   transposed, the transpose of the inverse transform at full length */
static void _op_fft_nmod(fft_small_op_t X, const ulong* a, ulong an,
                    ulong itrunc, int transposed, nmod_t mod,
                    const fft_small_plan_t P)
{
    ulong i;
    thread_pool_handle* handles = NULL;
    slong nworkers = 0;
    ulong nthreads;
    _op_fft_worker_struct args[MPN_CTX_NCRTS];

    FLINT_ASSERT(X->np == P->np && X->offset == P->offset &&
                 X->depth == P->depth && X->stride == P->stride);
    FLINT_ASSERT(itrunc % 8 == 0);
    FLINT_ASSERT(itrunc <= P->ztrunc);
    FLINT_ASSERT(!transposed || an <= P->ztrunc);

    /* thread over the primes; the threshold mirrors the fused drivers'
       stage 1 gate */
    if (P->np >= 2 && an >= 1500)
        nworkers = flint_request_threads(&handles, P->np);
    nthreads = FLINT_MIN((ulong)(nworkers + 1), P->np);

    for (i = 0; i < nthreads; i++)
    {
        _op_fft_worker_struct* W = args + i;
        W->X = X;
        W->a = a;
        W->an = an;
        W->itrunc = itrunc;
        W->transposed = transposed;
        W->mod = mod;
        W->P = P;
        W->start_pi = (i+0)*P->np/nthreads;
        W->stop_pi  = (i+1)*P->np/nthreads;
    }

    for (i = nthreads - 1; i > 0; i--)
        thread_pool_wake(global_thread_pool, handles[i - 1], 0,
                         _op_fft_worker_func, args + i);
    _op_fft_worker_func(args + 0);
    for (i = nthreads - 1; i > 0; i--)
        thread_pool_wait(global_thread_pool, handles[i - 1]);

    flint_give_back_threads(handles, nworkers);

    X->itrunc = itrunc;
    X->domain = FFT_SMALL_OP_PRIMAL;
}

void fft_small_fft_nmod(fft_small_op_t X, const ulong* a, ulong an,
                    ulong itrunc, nmod_t mod, const fft_small_plan_t P)
{
    _op_fft_nmod(X, a, an, itrunc, 0, mod, P);
}


/*
    Transposes of fft_small_fft_nmod and fft_small_ifft, for the transposed
    convolutions (see sd_fft_tr.c and sd_ifft_tr.c): fft_small_ifft_tr_nmod
    (here) converts (a, an) and applies the transpose of the inverse
    transform, leaving a PRIMAL operand ready for the pointwise stage;
    fft_small_fft_tr (in op.c) applies the transpose of the forward transform
    of an operand of truncation itrunc, consuming the normalization the
    pointwise stage applied. A transposed product is therefore

        fft_small_fft_nmod(B, b, bn, btrunc, mod, P);
        fft_small_ifft_tr_nmod(X, v, vn, mod, P);
        fft_small_op_mul(X, X, B, P);
        fft_small_fft_tr(X, atrunc, P);
        fft_small_export_nmod_range(z, X, 0, an, mod, P);

    which is the transpose, with respect to a, of the multiplication by b
    performed by fft_small_fft_nmod / fft_small_op_mul / fft_small_ifft /
    fft_small_export_nmod.
*/
void fft_small_ifft_tr_nmod(fft_small_op_t X, const ulong* a, ulong an,
                    nmod_t mod, const fft_small_plan_t P)
{
    _op_fft_nmod(X, a, an, P->ztrunc, 1, mod, P);
}

typedef struct {
    ulong* z;
    const fft_small_op_struct* X;
    nmod_t mod;
    const fft_small_plan_struct* P;
    fft_small_crt_func crt_fn;
    ulong start_zi;
    ulong stop_zi;
} _op_export_worker_struct;

static void _op_export_worker_func(void* varg)
{
    _op_export_worker_struct* W = (_op_export_worker_struct*) varg;
    const fft_small_plan_struct* P = W->P;

    W->crt_fn((void*) W->z, P->zl, W->start_zi, W->stop_zi,
              P->ffts + P->offset, W->X->data, P->stride,
              P->crts + P->offset, P->bound_c, NULL, &W->mod);
}

void fft_small_export_nmod_range(ulong* z, const fft_small_op_t X,
                    ulong zl, ulong zh, nmod_t mod, const fft_small_plan_t P)
{
    ulong i, o;
    thread_pool_handle* handles = NULL;
    slong nworkers = 0;
    ulong nthreads;
    _op_export_worker_struct args[8];

    FLINT_ASSERT(P->zl <= zl && zh <= P->zh);

    if (zl >= zh)
        return;

    if (P->np*(zh - zl) > 10000)
        nworkers = flint_request_threads(&handles, 8);
    nthreads = nworkers + 1;

    o = zl;
    for (i = 0; i < nthreads; i++)
    {
        _op_export_worker_struct* W = args + i;
        W->z = z - (zl - P->zl);   /* worker offsets by P->zl */
        W->X = X;
        W->mod = mod;
        W->P = P;
        W->crt_fn = _nmod_crt_fn(P->np);
        W->start_zi = o;
        ulong newo = n_round_down(zl + (i+1)*(zh - zl)/nthreads, BLK_SZ);
        o = i+1 < nthreads ? FLINT_MAX(o, newo) : zh;
        W->stop_zi = o;
    }

    for (i = nworkers; i > 0; i--)
        thread_pool_wake(global_thread_pool, handles[i - 1], 0,
                         _op_export_worker_func, args + i);
    _op_export_worker_func(args + 0);
    for (i = nworkers; i > 0; i--)
        thread_pool_wait(global_thread_pool, handles[i - 1]);

    flint_give_back_threads(handles, nworkers);
}

void fft_small_export_nmod(ulong* z, const fft_small_op_t X, nmod_t mod,
                    const fft_small_plan_t P)
{
    ulong i, o;
    thread_pool_handle* handles = NULL;
    slong nworkers = 0;
    ulong nthreads;
    _op_export_worker_struct args[8];

    /* the number of products accumulated onto one output coefficient is
       bounded by the len_bound the plan was created with; thread over
       BLK_SZ-aligned output ranges as the engine's stage 2 does, with
       the same rethread heuristic */
    if (P->np*(P->zh - P->zl) > 10000)
        nworkers = flint_request_threads(&handles, 8);
    nthreads = nworkers + 1;

    o = P->zl;
    for (i = 0; i < nthreads; i++)
    {
        _op_export_worker_struct* W = args + i;
        W->z = z;
        W->X = X;
        W->mod = mod;
        W->P = P;
        W->crt_fn = _nmod_crt_fn(P->np);
        W->start_zi = o;
        ulong newo = n_round_down(P->zl + (i+1)*(P->zh - P->zl)/nthreads, BLK_SZ);
        o = i+1 < nthreads ? FLINT_MAX(o, newo) : P->zh;
        W->stop_zi = o;
    }

    for (i = nworkers; i > 0; i--)
        thread_pool_wake(global_thread_pool, handles[i - 1], 0,
                         _op_export_worker_func, args + i);
    _op_export_worker_func(args + 0);
    for (i = nworkers; i > 0; i--)
        thread_pool_wait(global_thread_pool, handles[i - 1]);

    flint_give_back_threads(handles, nworkers);
}

/*
    The coefficients [zl, zh) of the product of (a, an) and (b, bn), from a
    convolution of length L = 2^depth chosen by _fft_small_plan_set_window:
    usually L >= an + bn - 1, but L may be shorter when the wrapped
    coefficients of the product modulo x^L - 1 all fall below zl.
*/
void _nmod_poly_mul_mid_mpn_ctx(
    ulong* z, ulong zl, ulong zh,
    const ulong* a, ulong an,
    const ulong* b, ulong bn,
    nmod_t mod,
    mpn_ctx_t R)
{
    ulong modbits = FLINT_BITS - mod.norm;
    ulong zn = an + bn - 1;
    ulong mn = n_min(an, bn);
    ulong mx = n_max(an, bn);
    ulong atrunc, btrunc;
    int squaring, success;
    fft_small_plan_t P;
    fft_small_conv_arg_struct A;

    FLINT_ASSERT(an > 0);
    FLINT_ASSERT(bn > 0);

    if (zl >= zh)
        return;

    if (zh > zn)
    {
        if (zl >= zn)
        {
            flint_mpn_zero(z, zh - zl);
            return;
        }

        flint_mpn_zero(z + zn - zl, zh - zn);
        zh = zn;
    }

    FLINT_ASSERT(zl < zh);
    FLINT_ASSERT(zh <= zn);

    squaring = (a == b) && (an == bn);

    atrunc = _fft_small_len_trunc(an);
    btrunc = _fft_small_len_trunc(bn);

    /* at most min(an, bn) products, each < 2^(2*modbits), accumulate onto
       one output coefficient */
    /* the direct-fft branch is gated on the longer operand, as in
       _nmod_poly_mul_tr_mpn_ctx; nothing here depends on the operand order */
    success = fft_small_plan_init_nmod(P, R, zl, zh, zn,
                        n_max(atrunc, btrunc), mn, 2*modbits, mod, mx);
    FLINT_ASSERT(success);
    (void) success;

    /* operand truncations are capped by the plan's transform: the
       wraparound importer folds any excess length */
    atrunc = n_min(atrunc, P->ztrunc);
    btrunc = n_min(btrunc, P->ztrunc);

    A.a = a; A.an = an; A.aaux = 0; A.atrunc = atrunc;
    A.b = b; A.bn = bn; A.baux = 0; A.btrunc = btrunc;
    A.bfft = NULL;
    A.bfft_stride = 0;
    A.squaring = squaring;
    A.transposed = 0;
    A.tune_n = mn;
    A.min_an_bn = mn;
    A.mod_fn = _mod_fn;
    A.s2_finish = NULL;
    A.crt_fn = _nmod_crt_fn(P->np);
    A.params = &mod;
    A.tuning = &_nmod_conv_tuning;

    _fft_small_conv(z, P, &A);

    fft_small_plan_clear(P);
}


/*
    The transposed product.

    Summary: writes to z the zn coefficients
       z[j] = sum_i a[i] * v[i + j], 0 <= j < zn,
    with v read as padded with zeros to length an + zn - 1 (so vn <= an + zn - 1
    is required). This is the transpose of the multiplication by (a, an) of
    polynomials of length zn, and costs the same; with a = rev(f) it is the
    range [fn - 1, fn - 1 + zn) of f*v, at a lower cost than the fn by vn
    product _nmod_poly_mul_mid_mpn_ctx would transform for the same range.

    ********************************************************************************

    Let `M` be the matrix of the multiplication by `(a, an)`,

        M : b  |-->  a * b,     R^zn -> R^(an + zn - 1).

    Then `M[r][j]` is the coefficient of `x^r` in `a*x^j`, that is
    `a[r - j]`, so its transpose maps `v` of length `an + zn - 1` to

        (M^T v)[j] = sum_r a[r - j] v[r] = sum_i a[i] v[i + j],
                                                       0 <= j < zn.

    This is what the function below computes, with `v` read as padded with
    zeros to length `an + zn - 1`. With `a` the reversal of `(f, fn)`,
    `a[i] = f[fn - 1 - i]`, it is

        sum_i f[fn - 1 - i] v[i + j] = (f*v)[fn - 1 + j],

    i.e. the range `[fn - 1, fn - 1 + zn)` of `f*v`: the middle product of
    `f` and `v` is the transpose of a multiplication by `f` reversed. This
    is Tellegen's principle applied to polynomial multiplication; see
    [BLS2003]_. The window `[zl, zh)` of any product is of this form once
    the operand it slides over is shifted and padded (see
    _nmod_poly_mulmid_fft_small in nmod_poly/mullow_fft_small.c).

    ********************************************************************************

    Transposing the multiplication pipeline

        b  |-->  export o ifft o (fft(a) *) o fft

    step by step gives

        v  |-->  export o fft_tr o (fft(a) *) o ifft_tr

    where the pointwise stage is unchanged (it is its own transpose, and it
    still applies the normalization exactly once) and only the two
    transforms are replaced by the transposed ones of sd_fft_tr.c and
    sd_ifft_tr.c. In particular `a` is transformed exactly as in a
    multiplication by `a`. Both ends of the pipeline are in the natural
    coefficient order: `ifft_tr` consumes `v` indexed by increasing powers
    of `x`, as `ifft` produces its output, and `fft_tr` produces the result
    indexed the same way as `fft` consumes its input. The
    bit-reversed-and-shuffled order only appears on the evaluation side,
    where all three of `fft(a)`, the output of `ifft_tr` and the input of
    `fft_tr` agree.

    This is what the convolution engine of conv_engine.c does in its
    transposed mode, so the transposed product is threaded, buffered and
    tuned exactly as the multiplication it transposes; the fused driver
    below follows _nmod_poly_mul_mid_mpn_ctx.

    ********************************************************************************

    The point for middle products is the transform length. The transposed
    pipeline transforms `a` of length `an` and produces `zn` coefficients
    from a convolution of length `an + zn - 1`, so it costs exactly one
    product of lengths `an` and `zn`. Computing the same range as a window
    of the full product, which is what

        _nmod_poly_mul_mid_mpn_ctx(z, fn - 1, fn - 1 + zn, v, vn, f, fn, mod, R)

    does, transforms `v` of length `vn` and needs a convolution of length
    `fn + vn - 1`. In the balanced case `fn = zn`, the usual one for Newton
    iteration and power series division, that is a transform of length
    `3*fn` against `2*fn`, so up to a factor of about 1.5 in the transform
    work and 2 in the pointwise and conversion work; the wraparound of
    _fft_small_plan_set_window closes part of that gap when a power of two
    sits in the right place. See nmod_poly/profile/p-mul_mulmid_tune.c.

    No wraparound is used here: the plan is built for the window `[0, zn)`
    of a convolution of length `an + zn - 1`, so _fft_small_plan_set_window
    chooses `2^depth` at least that length. The cyclic trick would fold `v`
    into fewer slots, and the transpose of a cyclic convolution is a cyclic
    correlation, which is not the transposed product.
*/
void _nmod_poly_mul_tr_mpn_ctx(
    ulong* z, ulong zn,
    const ulong* a, ulong an,
    const ulong* v, ulong vn,
    nmod_t mod,
    mpn_ctx_t R)
{
    ulong modbits = FLINT_BITS - mod.norm;
    ulong pn = an + zn - 1;          /* length of the transposed product */
    ulong atrunc, ztrunc;
    int success;
    fft_small_plan_t P;
    fft_small_conv_arg_struct A;

    FLINT_ASSERT(an > 0);
    FLINT_ASSERT(vn > 0);
    FLINT_ASSERT(zn > 0);
    FLINT_ASSERT(vn <= pn);

    atrunc = _fft_small_len_trunc(an);
    ztrunc = _fft_small_len_trunc(zn);

    /* the plan of the multiplication of lengths an and zn, exported on
       [0, zn): every output coefficient is a sum of at most an products,
       each < 2^(2*modbits). The direct-fft branch is gated on the longer
       length: its setup cost is paid against the whole transform, and
       one-sided shapes (zn or an < 100) are where it pays most */
    success = fft_small_plan_init_nmod(P, R, 0, zn, pn,
                        n_max(atrunc, ztrunc), an, 2*modbits, mod,
                        n_max(an, zn));
    FLINT_ASSERT(success);
    (void) success;

    /* the plan's depth fixes the granularity of the operand truncations */
    atrunc = n_min(_fft_small_op_trunc(an, P), P->ztrunc);
    ztrunc = n_min(_fft_small_op_trunc(zn, P), P->ztrunc);

    /* transposed: the engine's a = v goes through ifft_tr at full length
       and comes back through fft_tr truncated to ztrunc; its b = a is
       transformed as in a multiplication */
    A.a = v; A.an = vn; A.aaux = 0; A.atrunc = ztrunc;
    A.b = a; A.bn = an; A.baux = 0; A.btrunc = atrunc;
    A.bfft = NULL;
    A.bfft_stride = 0;
    A.squaring = 0;
    A.transposed = 1;
    A.tune_n = n_min(an, zn);
    A.min_an_bn = an;       /* every output is a sum of at most an products */
    A.mod_fn = _mod_fn;
    A.s2_finish = NULL;
    A.crt_fn = _nmod_crt_fn(P->np);
    A.params = &mod;
    A.tuning = &_nmod_conv_tuning;

    _fft_small_conv(z, P, &A);

    fft_small_plan_clear(P);
}

#if 0
static void _nmod_poly_mul_mod_xpnm1_naive(
    ulong* z, ulong zn,
    const ulong* a, ulong an,
    const ulong* b, ulong bn,
    ulong lgN,
    nmod_t mod,
    mpn_ctx_t R)
{
    ulong N = n_pow2(lgN);
    FLINT_ASSERT(zn <= N);

    ulong* t = FLINT_ARRAY_ALLOC(an + bn - 1, ulong);

    if (an >= bn)
        _nmod_poly_mul(t, a, an, b, bn, mod);
    else
        _nmod_poly_mul(t, b, bn, a, an, mod);

    for (ulong i = 0; i < zn; i++)
    {
        ulong c = 0;
        for (ulong j = i; j < an + bn - 1; j += N)
            c = nmod_add(c, t[j], mod);
        z[i] = c;
    }

    flint_free(t);
}
#endif

/*
Set `z` to the cyclic convolution of `a` and `b` modulo `mod`
with length `N = 2^depth`.
*/
static void _nmod_poly_mul_mod_xpnm1(
    ulong* z, ulong ztrunc,
    const ulong* a, ulong an,
    const ulong* b, ulong bn,
    ulong depth,
    nmod_t mod,
    mpn_ctx_t R)
{
    ulong N = n_pow2(depth);
    ulong modbits = FLINT_BITS - mod.norm;
    fft_small_plan_t P;
    fft_small_conv_arg_struct A;
    int success;

    FLINT_ASSERT(an > 0);
    FLINT_ASSERT(bn > 0);
    FLINT_ASSERT(ztrunc <= N);

    /* a cyclic convolution accumulates up to N products onto one output
       coefficient */
    success = fft_small_plan_init_nmod_cyclic(P, R, depth, ztrunc, N,
                                              2*modbits, mod);
    FLINT_ASSERT(success);
    (void) success;

    A.a = a; A.an = an; A.aaux = 0; A.atrunc = N;
    A.b = b; A.bn = bn; A.baux = 0; A.btrunc = N;
    A.bfft = NULL;
    A.bfft_stride = 0;
    A.squaring = 0;
    A.transposed = 0;
    A.tune_n = bn;
    A.min_an_bn = FLINT_MIN(an, bn);
    A.mod_fn = _mod_fn;
    A.s2_finish = NULL;
    A.crt_fn = _nmod_crt_fn(P->np);
    A.params = &mod;
    A.tuning = &_nmod_conv_tuning_no_rethread;

    _fft_small_conv(z, P, &A);

    fft_small_plan_clear(P);
}


/*
    Stores the transform of (b, bn), reduced modulo x^N - 1, for cyclic
    convolutions of length N = 2^depth, where depth >= 4, in which at most
    len_bound products accumulate onto one output coefficient: N covers any
    use, and a smaller bound, such as min(an, bn) when bn <= N and the other
    operand has length at most an <= N, may save a prime. btrunc is the
    truncation of the transforms the precomputation will be used with (rounded
    up as the operand truncations). direct_len gates the direct transform
    modulo a prime as in fft_small_plan_init_nmod (0 disables it). Returns 0 if
    no plan is available.
*/
int _mul_precomp_init(
    mul_precomp_struct* M,
    const ulong * b, ulong bn, ulong btrunc,
    ulong depth, ulong len_bound, ulong direct_len,
    nmod_t mod,
    mpn_ctx_t R)
{
    ulong modbits = FLINT_BITS - mod.norm;
    ulong N = n_pow2(depth);

    FLINT_ASSERT(depth >= 4);
    FLINT_ASSERT(bn > 0);

    /* the window [0, N) of a length N convolution: a cyclic plan */
    if (!fft_small_plan_init_nmod(M->P, R, 0, N, N, N, len_bound,
                                  2*modbits, mod, direct_len))
        return 0;

    M->bn = bn;
    M->btrunc = _fft_small_len_trunc(btrunc);
    M->transposed = 0;

    /* the forward transform reads b zero-padded (or folded) to this
       truncation, which does not change the transform */
    fft_small_op_init(M->bfft, M->P);
    fft_small_fft_nmod(M->bfft, b, bn, n_min(_fft_small_len_trunc(bn), N), mod, M->P);

    return 1;
}


/*
    Coefficients [zl, zh) of a * b, with b given by M; returns 0 if the
    window needs a transform longer than the one of M.

    M is only read: the window is fitted on a copy of the plan, which takes
    its scratch and its primes from R, so that concurrent calls with the
    same M are fine if each uses its own R, for instance the default context
    of its thread. R must have the primes of the context M was made with,
    as all default contexts do; a plan transforming directly modulo mod.n
    keeps its own tables.
*/
int _nmod_poly_mul_mid_precomp(
    ulong* z, ulong zl, ulong zh,
    const ulong* a, ulong an,
    const mul_precomp_struct* M,
    nmod_t mod,
    mpn_ctx_t R)
{
    ulong bn = M->bn;
    ulong zn = an + bn - 1;
    ulong atrunc;
    fft_small_plan_t P;
    fft_small_conv_arg_struct A;

    FLINT_ASSERT(an > 0);
    FLINT_ASSERT(bn > 0);
    FLINT_ASSERT(!M->transposed);
    FLINT_ASSERT(M->btrunc <= n_pow2(M->P->depth));

    if (zl >= zh)
        return 1;

    if (zh > zn)
    {
        if (zl >= zn)
        {
            flint_mpn_zero(z, zh - zl);
            return 1;
        }

        flint_mpn_zero(z + zn - zl, zh - zn);
        zh = zn;
    }

    FLINT_ASSERT(zl < zh);
    FLINT_ASSERT(zh <= zn);

    atrunc = n_min(_fft_small_len_trunc(an), n_pow2(M->P->depth));

    *P = *M->P;
    P->R = R;
    P->crts = R->crts;
    if (!P->use_direct_fft)
        P->ffts = R->ffts;

    if (!_fft_small_plan_fit_window(P, zl, zh, zn, atrunc, (ulong) 4))
        return 0;

    /* cap AFTER the window is fitted: fitting recomputes ztrunc, which
       can shrink on a reused plan, and coefficients past it cannot
       reach a linear window */
    atrunc = n_min(atrunc, P->ztrunc);

    A.a = a; A.an = an; A.aaux = 0; A.atrunc = atrunc;
    A.b = NULL; A.bn = bn; A.baux = 0; A.btrunc = 0;
    A.bfft = M->bfft->data;
    A.bfft_stride = M->bfft->stride;
    A.squaring = 0;
    A.transposed = 0;
    A.tune_n = bn;
    A.min_an_bn = FLINT_MIN(an, bn);
    A.mod_fn = _mod_fn;
    A.s2_finish = NULL;
    A.crt_fn = _nmod_crt_fn(P->np);
    A.params = &mod;
    A.tuning = &_nmod_conv_tuning_no_rethread;

    _fft_small_conv(z, P, &A);

    return 1;
}

/*
    The precomputation of the transposed products with a fixed (v, vn):
    the transposed products z[j] = sum_i r[i] v[i + j], 0 <= j < zn, for
    any (r, rn) with rn + zn - 1 <= vn (see _nmod_poly_mul_tr_mpn_ctx),
    are the transposes of products of length at most vn, which one plan
    for the window [0, vn) of a convolution of length vn serves, and they
    share the images of v through the transpose of the inverse transform,
    which are stored. With r the reversal of a, zero padded in front, this
    is a window of a * v: the middle products of evaluations and
    extrapolations at geometric progressions with a fixed polynomial v.
    Against the window of the full product with a stored transform of v
    (_mul_precomp_init), the transform has length about vn instead of a
    power of two, or of the length of the full product.

    len_bound bounds the number of nonzero coefficients of r (the number of
    products accumulating onto one output coefficient), and direct_len
    gates the direct transform modulo a prime as in
    fft_small_plan_init_nmod. Returns 0 if no plan is available.
*/
int _mul_precomp_init_tr(
    mul_precomp_struct* M,
    const ulong * v, ulong vn,
    ulong len_bound, ulong direct_len,
    nmod_t mod,
    mpn_ctx_t R)
{
    ulong modbits = FLINT_BITS - mod.norm;

    FLINT_ASSERT(vn > 0);

    if (!fft_small_plan_init_nmod(M->P, R, 0, vn, vn, _fft_small_len_trunc(vn),
                                  len_bound, 2*modbits, mod, direct_len))
        return 0;

    M->bn = vn;
    M->btrunc = 0;
    M->transposed = 1;

    fft_small_op_init(M->bfft, M->P);
    fft_small_ifft_tr_nmod(M->bfft, v, vn, mod, M->P);

    return 1;
}

/*
    The transposed product z[j] = sum_i r[i] v[i + j], 0 <= j < zn, with v
    given by M (_mul_precomp_init_tr); returns 0 if rn + zn - 1 exceeds
    the length of v. r must have at most len_bound nonzero coefficients
    (the len_bound of M), which leaves room for zeros in front of a
    reversed operand. The precomputed images of v are multiplied by the
    transform of r and the transpose of the forward transform gives the zn
    coefficients: two transforms, of the length of M. M is only read, as
    for _nmod_poly_mul_mid_precomp.
*/
int _nmod_poly_mul_tr_precomp(
    ulong* z, ulong zn,
    const ulong* r, ulong rn,
    const mul_precomp_struct* M,
    nmod_t mod,
    mpn_ctx_t R)
{
    fft_small_plan_t P;
    fft_small_conv_arg_struct A;

    FLINT_ASSERT(M->transposed);
    FLINT_ASSERT(rn > 0 && zn > 0);

    if (rn + zn - 1 > M->bn)
        return 0;

    *P = *M->P;
    P->R = R;
    P->crts = R->crts;
    if (!P->use_direct_fft)
        P->ffts = R->ffts;

    P->zl = 0;
    P->zh = zn;

    A.a = NULL; A.an = 0; A.aaux = 0;
    A.atrunc = n_min(_fft_small_op_trunc(zn, P), P->ztrunc);
    A.b = r; A.bn = rn; A.baux = 0;
    A.btrunc = n_min(_fft_small_op_trunc(rn, P), P->ztrunc);
    A.bfft = M->bfft->data;
    A.bfft_stride = M->bfft->stride;
    A.squaring = 0;
    A.transposed = 2;
    A.tune_n = n_min(rn, zn);
    A.min_an_bn = n_min(rn, M->P->bound_c);
    A.mod_fn = _mod_fn;
    A.s2_finish = NULL;
    A.crt_fn = _nmod_crt_fn(P->np);
    A.params = &mod;
    A.tuning = &_nmod_conv_tuning_no_rethread;

    _fft_small_conv(z, P, &A);

    return 1;
}


static void _nmod_poly_mul_mod_xpnm1_precomp(
    ulong* z, ulong ztrunc,
    const ulong* a, ulong an,
    mul_precomp_struct* M,
    nmod_t mod,
    mpn_ctx_t FLINT_UNUSED(R))
{
    fft_small_conv_arg_struct A;

    FLINT_ASSERT(an > 0);
    FLINT_ASSERT(R == M->P->R);
    FLINT_ASSERT(ztrunc <= M->btrunc);

    /* output window [0, ztrunc) of a cyclic convolution, with the ffts
       truncated at M->btrunc as the precomputed transform was */
    M->P->zl = 0;
    M->P->zh = ztrunc;
    M->P->zn = n_pow2(M->P->depth);
    M->P->ztrunc = M->btrunc;

    A.a = a; A.an = an; A.aaux = 0; A.atrunc = M->btrunc;
    A.b = NULL; A.bn = 0; A.baux = 0; A.btrunc = 0;
    A.bfft = M->bfft->data;
    A.bfft_stride = M->bfft->stride;
    A.squaring = 0;
    A.transposed = 0;
    A.tune_n = an;
    A.min_an_bn = an;  /* Todo: we might want to store bn and min by this here,
                          in case this is much smaller than an. */
    A.mod_fn = _mod_fn;
    A.s2_finish = NULL;
    A.crt_fn = _nmod_crt_fn(M->P->np);
    A.params = &mod;
    A.tuning = &_nmod_conv_tuning_no_rethread;

    _fft_small_conv(z, M->P, &A);
}



/* z -= a mod x^N-1, write coeffs [0,ztrunc) */
static void _nmod_poly_sub_mod_xpNm1(
    ulong* z, ulong ztrunc,
    const ulong* a, ulong an,
    ulong N, nmod_t mod)
{
    FLINT_ASSERT(ztrunc <= an);
    FLINT_ASSERT(ztrunc <= N);

    for (ulong i = 0; i < ztrunc; i++)
    {
        ulong k = i;
        ulong c = nmod_sub(a[k], z[i], mod);
        for (k += N; k < an; k += N)
            c = nmod_add(c, a[k], mod);
        z[i] = c;
    }
}

/*
    for divrem(a, b)

    an = length(a)
    bn = length(b)
    qn = length(q) = an - bn + 1

    choose a precision Bn of B(x) = B[0] + ... + B[Bn-1]*x^(Bn-1) with

        rev(B) = rev(b)^-1 mod x^Bn = B[Bn-1] + B[Bn-2]*x + ... + B[0]*x^(Bn-1)

    then
        (a[an-1] + a[an-2]*x + ... + a[an-qn]*x^(qn-1))
       *
        (B[Bn-1] + B[Bn-2]*x + ... + B[0]*x^(Bn-1))
       =
        q[qn-1] + q[qn-2]*x + ... + q[0]*x^(qn-1)  mod x^qn

    therefore need Bn >= qn, or, the same thing, Bn >= an - bn + 1

    in terms of non-reversed polys,

        _mul_mid(q, an+Bn-1-qn, an+Bn-1, a, an, B, Bn)

    or, the same thing,

        _mul_mid(q, Bn-1, Bn-1+qn, a+an-qn, qn, B, Bn)

    will calculate q. Then, find r via

        r = a - b*q mod x^N-1 where N >= bn - 1
*/
void _nmod_poly_divrem_mpn_ctx(
    ulong* q,
    ulong* r,
    const ulong* a, ulong an,
    const ulong* b, ulong bn,
    nmod_t mod,
    mpn_ctx_t R)
{
    ulong qn = an - bn + 1;

    FLINT_ASSERT(an >= bn);
    FLINT_ASSERT(bn > 1);
    FLINT_ASSERT(qn > 0);

    /* choose precision */
    ulong Bn = qn;

    ulong lgN = n_max(LG_BLK_SZ, n_clog2(bn-1));
    ulong N = n_pow2(lgN);

    ulong* B = FLINT_ARRAY_ALLOC(Bn, ulong);
    ulong* t = FLINT_ARRAY_ALLOC(FLINT_MAX(N, bn), ulong);

    _nmod_poly_reverse(t, b, bn, bn);
    _nmod_poly_inv_series(B, t, bn, Bn, mod);
    _nmod_poly_reverse(B, B, Bn, Bn);

    _nmod_poly_mul_mid_mpn_ctx(q, Bn-1, Bn-1+qn, a+an-qn, qn, B, Bn, mod, R);
    _nmod_poly_mul_mod_xpnm1(r, bn-1, q, qn, b, bn, lgN, mod, R);
    _nmod_poly_sub_mod_xpNm1(r, bn-1, a, an, N, mod);

    flint_free(B);
    flint_free(t);
}

void _nmod_poly_divrem_precomp_init(
    nmod_poly_divrem_precomp_struct* M,
    const ulong* b, ulong bn,
    ulong Bn,
    nmod_t mod,
    mpn_ctx_t R)
{
    ulong* B = FLINT_ARRAY_ALLOC(Bn, ulong);
    ulong* t = FLINT_ARRAY_ALLOC(bn, ulong);

    _nmod_poly_reverse(t, b, bn, bn);
    _nmod_poly_inv_series(B, t, bn, Bn, mod);
    _nmod_poly_reverse(B, B, Bn, Bn);


    /* TODO: both precomputations keep the bound 2^depth on the number of
       accumulated products and do not allow the direct transform modulo a
       prime (direct_len = 0), as before _mul_precomp_init had these
       parameters. Passing a nonzero direct_len (and a tighter bound where
       the operand lengths allow) would let a prime modulus that fft_small
       transforms directly use a single prime instead of two or three;
       to be measured. */
    ulong lgQ = n_max(LG_BLK_SZ, n_clog2(2*Bn-1));
    int success = _mul_precomp_init(M->quo_maker, B, Bn, Bn, lgQ,
                                     n_pow2(lgQ), 0, mod, R);

    ulong lgN = n_max(LG_BLK_SZ, n_clog2(bn-1));
    ulong N = n_pow2(lgN);
    success = success && _mul_precomp_init(M->rem_maker, b, bn, N, lgN,
                                           N, 0, mod, R);
    FLINT_ASSERT(success);
    (void) success;

    flint_free(B);
    flint_free(t);
}


int _nmod_poly_divrem_precomp(
    ulong* q,
    ulong* r,
    const ulong* a, ulong an,
    nmod_poly_divrem_precomp_struct* M,
    nmod_t mod,
    mpn_ctx_t R)
{
    ulong Bn = M->quo_maker->bn;
    ulong bn = M->rem_maker->bn;
    ulong qn = an-bn+1;

    FLINT_ASSERT(an >= bn);
    FLINT_ASSERT(bn > 1);
    FLINT_ASSERT(qn > 0);

    /* the stored inverse has precision Bn */
    if (qn > Bn)
        return 0;

    if (!_nmod_poly_mul_mid_precomp(q, Bn-1, Bn-1+qn, a+an-qn, qn, M->quo_maker, mod, R))
        return 0;
    _nmod_poly_mul_mod_xpnm1_precomp(r, bn-1, q, qn, M->rem_maker, mod, R);
    _nmod_poly_sub_mod_xpNm1(r, bn-1, a, an, n_pow2(M->rem_maker->P->depth), mod);
    return 1;
}



/*
definition of _mul_mid(z, zl, zh, a, an, b, bn)

              h
[sum z[i]*x^i]  :=  sum  z[i]*x^(i-l)
  i           l    l<=i<h

i.e. the coeffs in [zl, zh)
*/
void _nmod_poly_mul_mid_classical(
    ulong* z, slong zl, slong zh,
    const ulong* a, slong an,
    const ulong* b, slong bn,
    nmod_t mod)
{
    for (slong i = zl; i < zh; i++)
    {
        slong jstart = z_max(0, i - (bn - 1));
        slong jstop = z_min(i + 1, an);
        ulong zi = 0;
        for (slong j = jstart; j < jstop; j++)
            zi = nmod_addmul(zi, a[j], b[i - j], mod);
        z[i - zl] = zi;
    }
}

/*
**** Karasuba ****

with deg(Ai), deg(Bi) < k, consider

P := (A0 + A1*x^k)*(B0 + B1*x^k)

define

P2 = A1*B1
P1 = (A0 + A1)*(B0 + B1)
P0 = A0*B0

then

# P = P0 + (P1 - P2 - P0)*x^k + P2*x^2k
*/




/*
**** Karasuba for middle product ****

with deg(Ai), deg(Bi) < k, consider
P := (A0 + A1*x^k + A2*x^2k + A3*x^3k)*(B0 + B1*x^k)

                            h
we would like to compute [P]
                            l
define

P0 := ((A0 + A1) + (A1 + A2)*x^k)*(B1)
P1 := (A1 + A2*x^k)*(B0 - B1)
P2 := ((A1 + A2) + (A2 + A3)*x^k)*(B0)

so that

P0 = (A0*B1 + A1*B1) + (A1*B1 + A2*B1)*x^k
P1 = (A1*B0 - A1*B1) + (A2*B0 - A2*B1)*x^k
P2 = (A1*B0 + A2*B0) + (A2*B0 + A3*B0)*x^k

and

   h            2k             h-2k
[P]  = [P0 + P1]    + [P2 - P1]    * x^(3k-l)
   l            l-k            k

In order to calculate the rhs, we need

    2k         h-2k        max(2k,h-2k)
[P0]     , [P2]     ,  [P1]
    l-k        k           min(l-k,k)


requires k <= l < 3k < h <= 4k

*/


static void _nmod_poly_mul_mid_unbalanced(
    ulong* z, slong zl, slong zh,
    const ulong* a, slong an,
    const ulong* b, slong bn,
    nmod_t mod)
{
    FLINT_ASSERT(zl < zh);
    FLINT_ASSERT(bn < an);
    flint_mpn_zero(z, zh - zl);

    ulong* t = FLINT_ARRAY_ALLOC(2*bn, ulong);

    slong i;
    for (i = 0; i*bn < an; i++)
    {
        slong zl_new, zh_new, an_new;

        // produces coefficient for powers x^[zl_new + i*k, zh_new + i*k)

        zl_new = z_max(zl - i*bn, 0);
        an_new = z_min(bn, an - i*bn);
        zh_new = z_min(zh - i*bn, an_new + bn - 1);

        _nmod_poly_mul_mid(t, zl_new, zh_new, a + i*bn, an_new, b, bn, mod);
        ulong* Z = z + zl_new + i*bn - zl;
        _nmod_vec_add(Z, Z, t, zh_new - zl_new, mod);
    }

    flint_free(t);
}


void _nmod_poly_mul_mid(
    ulong* z, slong zl, slong zh,
    const ulong* a, slong an,
    const ulong* b, slong bn,
    nmod_t mod)
{
    if (zl >= zh)
        return;

    if (an < bn)
    {
        FLINT_SWAP(const ulong *, a, b);
        FLINT_SWAP(ulong, an, bn);
    }

    if (zl > bn - 1)
    {
        if (an > zl - (bn - 1))
        {
            an -= zl - (bn - 1);
            a  += zl - (bn - 1);
            zh -= zl - (bn - 1);
            zl -= zl - (bn - 1);
            _nmod_poly_mul_mid(z, zl, zh, a, an, b, bn, mod);
        }
        else
        {
            flint_mpn_zero(z, zh - zl);
        }
        return;
    }

    if (zh < an)
    {
        an = zh;
        _nmod_poly_mul_mid(z, zl, zh, a, an, b, bn, mod);
        return;
    }

    FLINT_ASSERT(zl < bn && bn <= an && an <= zh);

    if (an >= 2*bn)
    {
        _nmod_poly_mul_mid_unbalanced(z, zl, zh, a, an, b, bn, mod);
        return;
    }

    if (zl < an + bn + 1)
    {
        if (zh > 0)
        {
            /*
                middle product or three pieces
                +----------------+
                |             |\ |
                |      1      | \|
                |             |3 |
                |-------------+--|
                |\     2         |
                +----------------+
            */
        }
        else if (0)
        {
            /*
                two pieces
                +----------------+
                |             |\ |
                |             | \|
                |      1      |2 |
                |             |  |
                |             |  |
                +-------------+--+
            */
        }
        else
        {
            /*
                two pieces
                +----------------+
                |           \    |
                |            \   |
                |------------ \  |
                |              \ |
                |               \|
                +----------------+
            */
        }
    }
    else
    {
        if (zh > 0)
        {
            /*
                two pieces
                +----------------+
                | |              |
                | |              |
                |1|      2       |
                | |              |
                |\|              |
                +----------------+
            */

        }
    }

    _nmod_poly_mul_mid_classical(z, zl, zh, a, an, b, bn, mod);
    return;
}
