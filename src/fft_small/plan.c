/*
    Copyright (C) 2022 Daniel Schultz
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "ulong_extras.h"
#include "nmod.h"
#include "fft_small.h"

int _fft_small_plan_set_bound(fft_small_plan_t P, ulong c, ulong e, ulong np_max)
{
    mpn_ctx_struct * R = P->R;
    ulong np;

    P->bound_c = c;
    P->bound_e = e;

    for (np = 1; ; np++)
    {
        if (np > np_max)
            return 0;

        if (flint_mpn_cmp_ui_2exp(crt_data_prod_primes(R->crts + np - 1),
              R->crts[np - 1].coeff_len, c, e) >= 0)
        {
            break;
        }
    }

    P->np = np;
    P->offset = 0;
    P->ffts = R->ffts;
    P->crts = R->crts;
    P->use_direct_fft = 0;

    return 1;
}

void _fft_small_plan_set_window(fft_small_plan_t P,
                    ulong zl, ulong zh, ulong zn, ulong xtrunc_max,
                    ulong min_depth)
{
    ulong ztrunc, depth, i;

    FLINT_ASSERT(zl < zh);
    FLINT_ASSERT(zh <= zn);

    /*
        if there is a power of two 2^d between zh and zn with good wrap around
            i.e. max(operand lengths, zh) <= 2^d <= zn with zn - 2^d <= zl
        then use d as the depth, otherwise the usual with no wrap around.
        Truncation granularity is one block at block depths and the
        transform length itself below one block, where sd_fft_trunc
        dispatches to the register basecases. */
    depth = n_flog2(zn);
    i = n_pow2(depth);
    if (xtrunc_max <= i && zh <= i && i <= zn && zn <= zl + i)
    {
        depth = n_max(depth, min_depth);
        /* if the floor raised the depth, the transform length grows
           with it; the wraparound conditions are monotone in i, so
           they continue to hold */
        ztrunc = n_pow2(depth);
    }
    else
    {
        /* the small headroom above zn covers the signed
           representation's bias digits, which span a few chunks past
           the product and were absorbed by block rounding before */
        depth = n_max(min_depth, n_clog2(zn + 4));
        ztrunc = (depth < LG_BLK_SZ) ? n_pow2(depth)
                                     : n_round_up(zn, BLK_SZ);
    }

    P->zl = zl;
    P->zh = zh;
    P->zn = zn;
    P->ztrunc = ztrunc;
    P->depth = depth;
    P->bits = 0;
}

int _fft_small_plan_fit_window(fft_small_plan_t P,
                    ulong zl, ulong zh, ulong zn, ulong xtrunc_max,
                    ulong min_depth)
{
    ulong N = n_pow2(P->depth);
    ulong ztrunc;

    FLINT_ASSERT(zl < zh);
    FLINT_ASSERT(zh <= zn);

    /* a fixed depth below the conversions' floor cannot serve */
    if (P->depth < min_depth)
        return 0;

    ztrunc = (P->depth < LG_BLK_SZ) ? N : n_round_up(zn, BLK_SZ);

    if (xtrunc_max <= N && zh <= N && N <= zn && zn <= zl + N)
    {
        ztrunc = N;
    }
    else if (ztrunc <= N)
    {
    }
    else
    {
        return 0;
    }

    P->zl = zl;
    P->zh = zh;
    P->zn = zn;
    P->ztrunc = ztrunc;

    return 1;
}

void _fft_small_plan_set_normalizers(fft_small_plan_t P)
{
    ulong i, m, cop;
    ulong np = P->np;
    ulong depth = P->depth;

    /* slabs padded to whole blocks: block-granular readers read a
       full block and discard beyond the emitted range; op init zeroes
       the pad. No cost at block depths. */
    P->stride = n_round_up(sd_fft_ctx_data_size(depth), BLK_SZ);

    for (i = 0; i < np; i++)
    {
        ulong ioff = i + P->offset;
        sd_fft_ctx_struct* Q = P->ffts + ioff;

        cop = np == 1 ? 1 : *crt_data_co_prime_red(P->crts + np - 1, ioff);
        NMOD_RED2(m, cop >> (FLINT_BITS - depth), cop << depth, Q->mod);
        m = nmod_inv(m, Q->mod);

        P->m[i] = m;
    }
}

/*
    Plan for mpn (integer) operands: choose a prime count np in [4, 8]
    (the range the mpn chinese remaindering implements) and an even chunk
    size `bits` with

        len_bound * ceil(FLINT_BITS*min(an_max, bn_max)/bits) * 2^(2*bits)
            <= prod_of_primes(np),

    i.e. enough room for len_bound products accumulated in transform
    space, maximizing bits for each np and choosing np by the same
    transform-cost heuristic np*depth*ztrunc the multiplication driver
    scores profiles with. The geometry is the full product of the
    maximal operand lengths; zl/zh/zn are in coefficient units and are
    not consulted by fft_small_export_mpn, which windows in limbs.
*/
static int _plan_init_mpn_np2(fft_small_plan_t P, mpn_ctx_t R,
                    ulong an_max, ulong bn_max, ulong len_bound);

static int
_plan_init_mpn_wide(fft_small_plan_t P, mpn_ctx_t R,
                    ulong an_max, ulong bn_max, ulong len_bound)
{
    ulong np, best_np = 0, best_bits = 0;
    ulong min_bn = FLINT_MIN(an_max, bn_max);
    double best_score = 0.0;

    FLINT_ASSERT(an_max > 0);
    FLINT_ASSERT(bn_max > 0);
    FLINT_ASSERT(len_bound > 0);

    /* prefer (np, bits) pairs with a vectorized packing profile: the
       operand conversions otherwise fall back to the scalar packing
       path, which costs more than the transforms it feeds */
    {
        ulong i;
        for (i = 0; i < R->profiles_size; i++)
        {
            ulong bits = R->profiles[i].bits;
            ulong pnp = R->profiles[i].np;
            const crt_data_struct* C = R->crts + pnp - 1;
            ulong t, thi, alen, blen, zlen, ztrunc, depth;
            double score;

            if (pnp < 4)
                continue;

            umul_ppmm(thi, t, len_bound, n_cdiv(FLINT_BITS*min_bn, bits));
            if (thi != 0 ||
                flint_mpn_cmp_ui_2exp(crt_data_prod_primes(C),
                                      C->coeff_len, t, 2*bits) < 0)
                continue;

            alen = n_cdiv(FLINT_BITS*an_max, bits);
            blen = n_cdiv(FLINT_BITS*bn_max, bits);
            zlen = alen + blen - 1;
            /* the (np, bits) selection deliberately scores at block
               granularity even though set_window will shorten the
               chosen plan: scoring short-aware makes the model trade
               prime count for depth (e.g. np 7 at depth 5 over np 4
               at depth 6), which measures slower -- the crt and
               conversion costs it underweights dominate at these
               sizes. Block scoring keeps the minimal-np narrow-chunk
               choices that measure best. */
            ztrunc = n_round_up(zlen, BLK_SZ);
            depth = n_max(LG_BLK_SZ, n_clog2(ztrunc));

            {
                double ratio = (double) ztrunc / (double) n_pow2(depth);
                score = (1 - 0.25*ratio) * (double) pnp * (double) depth
                                         * (double) ztrunc;
            }
            if (best_np == 0 || score < best_score)
            {
                best_np = pnp;
                best_bits = bits;
                best_score = score;
            }
        }

        if (best_np != 0)
            goto have_choice;
    }

    /* np = 3 has no vectorized packing profile: its imports take the
       slow per-coefficient path, which the narrow chunks (65..69 bits
       against the 3-prime product) keep affordable at small sizes */
    for (np = 3; np <= MPN_CTX_NCRTS; np++)
    {
        const crt_data_struct* C = R->crts + np - 1;
        ulong bits, alen, blen, zlen, ztrunc, depth;
        double score;

        /* the slow packing reads two_pow up to index bits - 32, and the
           reconstruction requires whole limbs to fit in a chunk */
        for (bits = (MPN_CTX_TWO_POWER_TAB_SIZE + 31) & ~UWORD(1);
             bits > 64;
             bits -= 2)
        {
            ulong t, thi;

            umul_ppmm(thi, t, len_bound, n_cdiv(FLINT_BITS*min_bn, bits));
            if (thi != 0)
                continue;

            if (flint_mpn_cmp_ui_2exp(crt_data_prod_primes(C),
                                      C->coeff_len, t, 2*bits) >= 0)
            {
                break;
            }
        }

        if (bits <= 64)
            continue;

        alen = n_cdiv(FLINT_BITS*an_max, bits);
        blen = n_cdiv(FLINT_BITS*bn_max, bits);
        zlen = alen + blen - 1;
        ztrunc = n_round_up(zlen, BLK_SZ);   /* scored at block, as above */
        depth = n_max(LG_BLK_SZ, n_clog2(ztrunc));

        {
            double ratio = (double) ztrunc / (double) n_pow2(depth);
            score = (1 - 0.25*ratio) * (double) np * (double) depth
                                     * (double) ztrunc;
        }
        if (best_np == 0 || score < best_score)
        {
            best_np = np;
            best_bits = bits;
            best_score = score;
        }
    }

    if (best_np == 0)
        return 0;

have_choice:
    {
        ulong bits = best_bits;
        ulong alen = n_cdiv(FLINT_BITS*an_max, bits);
        ulong blen = n_cdiv(FLINT_BITS*bn_max, bits);
        ulong zlen = alen + blen - 1;

        P->R = R;
        P->sign = 0;

        /* the block-rounded operand bound keeps the wraparound branch
           off for mpn plans, whose chunked conversions assume linear
           transforms; short depths come from the non-wrap branch */
        _fft_small_plan_set_window(P, 0, zlen, zlen,
                    n_max(n_round_up(alen, BLK_SZ), n_round_up(blen, BLK_SZ)),
                    (ulong) 4);

        P->bits = bits;
        P->bound_c = len_bound;
        P->bound_e = 2*bits;
        P->np = best_np;
        P->offset = 0;
        P->ffts = R->ffts;
        P->crts = R->crts;
        P->use_direct_fft = 0;

        _fft_small_plan_set_normalizers(P);
    }

    return 1;
}

/* Return whether to compute the FFT directly modulo `mod`, rather than
   modulo precomputed FFT primes followed by CRT reconstruction.
   Direct computation requires `mod.n` to be prime and `mod.n - 1`
   to have sufficiently high 2-valuation.
   This is usually faster when CRT would need multiple primes, but it pays
   the setup cost of initializing and clearing a temporary `sd_fft_ctx_t`.  */
static int
_nmod_poly_should_directly_fft(ulong bn, ulong depth, nmod_t mod)
{
    if (bn < 1500)
        return 0;

    if (mod.n <= 2 || mod.n > (UWORD(1) << 50))
        return 0;

    if (NMOD_BITS(mod) < 20)
        return 0;

    if (!fft_small_mulmod_satisfies_bounds(mod.n))
        /* should be implied by the 50-bit bound, but just in case */
        return 0;

    if (depth > SD_FFT_CTX_W2TAB_SIZE)
        /* unlikely, the convolution length would have to be massive */
        return 0;

    if (n_trailing_zeros(mod.n - 1) < depth)
        return 0;

    return n_is_prime(mod.n); /* check the most expensive condition last */
}

/* np/offset selection shared by the nmod plans. len_bound * 2^prod_bits
   is guaranteed to fit four primes for prod_bits <= 2*modbits and any
   realistic len_bound; larger accumulation bounds passed in by bilinear
   applications can use up to MPN_CTX_NCRTS primes, beyond which we
   return 0. */
static int _fft_small_plan_set_bound_nmod(fft_small_plan_t P,
                    ulong len_bound, ulong prod_bits, nmod_t mod,
                    ulong direct_len)
{
    mpn_ctx_struct * R = P->R;
    ulong modbits = FLINT_BITS - mod.norm;
    ulong i;

    P->bound_c = len_bound;
    P->bound_e = prod_bits;

    /* first see if mod.n is one of R->ffts[i].mod.n */
    if (modbits == 50)
    {
        for (i = 0; i < MPN_CTX_NCRTS; i++)
        {
            if (mod.n == R->ffts[i].mod.n)
            {
                P->np = 1;
                P->offset = i;
                P->ffts = R->ffts;
                P->crts = R->crts;
                P->use_direct_fft = 0;
                return 1;
            }
        }
    }

    if (direct_len != 0 &&
        _nmod_poly_should_directly_fft(direct_len, P->depth, mod))
    {
        sd_fft_ctx_init_prime(P->direct_fft, mod.n);
        P->np = 1;
        P->offset = 0;
        P->ffts = P->direct_fft;
        P->crts = R->crts;
        P->use_direct_fft = 1;
        return 1;
    }

    /* need prod_of_primes >= len_bound * 2^prod_bits */
    return _fft_small_plan_set_bound(P, len_bound, prod_bits, MPN_CTX_NCRTS);
}

int fft_small_plan_init_nmod(fft_small_plan_t P, mpn_ctx_t R,
                    ulong zl, ulong zh, ulong zn, ulong xtrunc_max,
                    ulong len_bound, ulong prod_bits, nmod_t mod,
                    ulong direct_len)
{
    P->R = R;
    P->sign = 0;

    /* the nmod conversions are per-coefficient and handle short
       transforms; the register basecases serve depths from 4 */
    _fft_small_plan_set_window(P, zl, zh, zn, xtrunc_max, (ulong) 4);

    if (!_fft_small_plan_set_bound_nmod(P, len_bound, prod_bits, mod, direct_len))
        return 0;

    _fft_small_plan_set_normalizers(P);
    return 1;
}

int fft_small_plan_init_nmod_cyclic(fft_small_plan_t P, mpn_ctx_t R,
                    ulong depth, ulong zh,
                    ulong len_bound, ulong prod_bits, nmod_t mod)
{
    ulong N = n_pow2(depth);

    FLINT_ASSERT(depth >= 4);
    FLINT_ASSERT(zh <= N);

    P->R = R;
    P->sign = 0;
    P->depth = depth;
    P->zl = 0;
    P->zh = zh;
    P->zn = N;
    P->ztrunc = N;

    if (!_fft_small_plan_set_bound_nmod(P, len_bound, prod_bits, mod, 0))
        return 0;

    _fft_small_plan_set_normalizers(P);
    return 1;
}

int fft_small_plan_init_mpn(fft_small_plan_t P, mpn_ctx_t R,
                    ulong an_max, ulong bn_max, ulong len_bound,
                    int is_signed)
{
    ulong len_bound_with_sign = is_signed ? 2 * len_bound : len_bound;

    if (len_bound <= 2
        && FLINT_MAX(an_max, bn_max) <= FLINT_MPN_MUL_NP2_MAX_AN
        && _plan_init_mpn_np2(P, R, an_max, bn_max, len_bound_with_sign))
        return 1;

    return _plan_init_mpn_wide(P, R, an_max, bn_max, len_bound_with_sign);
}

static int _plan_init_mpn_np2(fft_small_plan_t P, mpn_ctx_t R,
                    ulong an_max, ulong bn_max, ulong len_bound)
{
    ulong bits, depth, alen, blen, zlen;

    if (!_fft_small_np2_choose(R, an_max, bn_max, len_bound,
                               &bits, &depth))
        return 0;

    if (!_plan_init_mpn_wide(P, R, an_max, bn_max, len_bound))
        return 0;
    alen = n_cdiv(FLINT_BITS * an_max, bits);
    blen = n_cdiv(FLINT_BITS * bn_max, bits);
    zlen = alen + blen - 1;

    P->np = 2;
    P->offset = 0;
    P->bits = bits;
    P->depth = depth;
    P->stride = sd_fft_ctx_data_size(depth);
    P->zl = 0;
    P->zh = zlen;
    P->zn = zlen;
    P->ztrunc = n_round_up(zlen, BLK_SZ);
    P->m[0] = 1;
    P->m[1] = 1;
    return 1;
}

int fft_small_plan_init_mpn_np(fft_small_plan_t P, mpn_ctx_t R,
                    ulong an_max, ulong bn_max, ulong len_bound,
                    int is_signed, ulong np)
{
    ulong len_bound_with_sign = is_signed ? 2 * len_bound : len_bound;

    if (np == 0)
        return fft_small_plan_init_mpn(P, R, an_max, bn_max,
                                       len_bound, is_signed);
    if (np == 2)
        return _plan_init_mpn_np2(P, R, an_max, bn_max, len_bound_with_sign);

    return _plan_init_mpn_wide(P, R, an_max, bn_max, len_bound_with_sign);
}

/*
    Plan for cyclic integer convolutions: products of nonnegative operands
    of at most *nn limbs, reduced mod 2^(FLINT_BITS * *nn) - 1. Chooses a
    packing profile (np, bits) and a depth with

        bits * 2^depth == FLINT_BITS * nn,

    rounding *nn up to the nearest admissible length (written back): the
    operands then split into exactly 2^depth chunks, and the coefficient
    convolution's wrap mod x^(2^depth) - 1 is the integer wrap mod
    B^nn - 1. Operands longer than nn limbs are not folded by the packing
    and must be reduced by the caller.

    A full cyclic convolution accumulates up to 2^depth chunk products on
    one coefficient, so capacity is checked for len_bound * 2^depth
    products of 2*bits-bit values; as in the linear plans, len_bound > 1
    accounts for bilinear accumulation done in transform space.

    fft_small_export_mpn is destructive on the transform data, as
    output conversion is throughout fft_small; it returns a
    nonnegative representative
    sum(c_i 2^(bits i)) < len_bound * 2^(FLINT_BITS*nn + 2*bits + depth),
    exported exactly by any zn with

        FLINT_BITS*zn >= FLINT_BITS*nn + 2*bits + depth
                         + FLINT_BIT_COUNT(len_bound);

    the caller folds the zn - nn top limbs back around to reduce.

    Among the admissible profiles the plan takes the least padded length,
    then the fewest primes: the narrow-chunk minimal-np choices are the
    ones the linear driver's measurements favor at equal transform cost
    (see fft_small_plan_init_mpn).
*/
int fft_small_plan_init_mpn_cyclic(fft_small_plan_t P, mpn_ctx_t R,
                    ulong * nn, ulong len_bound)
{
    ulong best_np = 0, best_bits = 0, best_depth = 0, best_nn = 0;
    ulong i;

    FLINT_ASSERT(*nn > 0);
    FLINT_ASSERT(len_bound > 0);

    for (i = 0; i < R->profiles_size; i++)
    {
        ulong bits = R->profiles[i].bits;
        ulong pnp = R->profiles[i].np;
        const crt_data_struct* C = R->crts + pnp - 1;
        ulong depth, nn1, t, thi;

        if (pnp < 4)
            continue;

        /* smallest depth >= LG_BLK_SZ with bits * 2^depth >=
           FLINT_BITS * nn; bits divisible by 4 and 2^depth >= 16 keep
           bits * 2^depth divisible by FLINT_BITS, so nn1 is exact */
        depth = n_max(LG_BLK_SZ,
                      n_clog2(n_cdiv(FLINT_BITS * (*nn), bits)));
        if (bits >> (FLINT_BITS - depth) != 0)
            continue;
        nn1 = (bits << depth) / FLINT_BITS;

        /* room for len_bound * 2^depth products of 2*bits-bit chunks */
        umul_ppmm(thi, t, len_bound, n_pow2(depth));
        if (thi != 0 ||
            flint_mpn_cmp_ui_2exp(crt_data_prod_primes(C),
                                  C->coeff_len, t, 2*bits) < 0)
            continue;

        /* least padded length first, but np <= 4 stays in the
           single-vector packing/crt regime and measures ~10% faster
           per limb, so it is allowed ~6% extra padding (fold + linear
           pipeline timings across 83k-1M bit moduli) */
        {
            ulong score = nn1 * (pnp <= 4 ? 16 : 17);
            ulong best_score = best_nn * (best_np <= 4 ? 16 : 17);

            if (best_np == 0 || score < best_score ||
                (score == best_score && pnp < best_np))
            {
                best_np = pnp;
                best_bits = bits;
                best_depth = depth;
                best_nn = nn1;
            }
        }
    }

    /* two-prime candidate: nn1 = bits 2^depth / 64 minimal subject to
       len_bound 2^depth 2^(2 bits) <= p1 p2; cheaper conversions and
       Garner-only crt outweigh the longer transforms (per-shape op
       profiling), which the score constant 13 against the wide 16/17
       reflects */
    {
        ulong d, b;
        int stop = 0;

        for (d = n_max(LG_BLK_SZ, 6); d <= 26 && !stop; d++)
        {
            for (b = 45; b >= 40; b--)
            {
                ulong nn1 = (b << d) >> 6;
                if (nn1 < *nn)
                    continue;
                if (2*b + d + FLINT_BIT_COUNT(len_bound) > R->np2_prodbits)
                    continue;
                if (best_np == 0 ||
                    13 * nn1 < (best_np <= 4 ? 16 : 17) * best_nn)
                {
                    best_np = 2;
                    best_bits = b;
                    best_depth = d;
                    best_nn = nn1;
                }
                stop = 1;   /* larger d only pads more */
                break;
            }
        }
    }

    if (best_np == 0)
        return 0;

    *nn = best_nn;

    P->R = R;
    P->sign = 0;
    P->depth = best_depth;
    P->zl = 0;
    P->zh = n_pow2(best_depth);
    P->zn = n_pow2(best_depth);
    P->ztrunc = n_pow2(best_depth);
    P->bits = best_bits;
    /* stored as in fft_small_plan_init_mpn: the convolution-length
       factor entered the capacity check above but is not folded into
       bound_c */
    P->bound_c = len_bound;
    P->bound_e = 2*best_bits;
    P->np = best_np;
    P->offset = 0;
    P->ffts = R->ffts;
    P->crts = R->crts;
    P->use_direct_fft = 0;

    _fft_small_plan_set_normalizers(P);

    /* the two-prime export applies its own normalization and Garner
       constants; the pointwise stays unscaled */
    if (P->np == 2)
    {
        P->m[0] = 1;
        P->m[1] = 1;
    }

    return 1;
}

void fft_small_plan_clear(fft_small_plan_t P)
{
    if (P->use_direct_fft)
        sd_fft_ctx_clear(P->direct_fft);
}
