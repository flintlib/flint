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
                    ulong zl, ulong zh, ulong zn, ulong xtrunc_max)
{
    ulong ztrunc, depth, i;

    FLINT_ASSERT(zl < zh);
    FLINT_ASSERT(zh <= zn);

    ztrunc = n_round_up(zn, BLK_SZ);
    /*
        if there is a power of two 2^d between zh and zn with good wrap around
            i.e. max(operand lengths, zh) <= 2^d <= zn with zn - 2^d <= zl
        then use d as the depth, otherwise the usual with no wrap around
    */
    depth = n_flog2(zn);
    i = n_pow2(depth);
    if (xtrunc_max <= i && zh <= i && i <= zn && zn <= zl + i)
    {
        ztrunc = i;
    }
    else
    {
        depth = n_max(LG_BLK_SZ, n_clog2(ztrunc));
    }

    P->zl = zl;
    P->zh = zh;
    P->zn = zn;
    P->ztrunc = ztrunc;
    P->depth = depth;
    P->bits = 0;
}

int _fft_small_plan_fit_window(fft_small_plan_t P,
                    ulong zl, ulong zh, ulong zn, ulong xtrunc_max)
{
    ulong N = n_pow2(P->depth);
    ulong ztrunc;

    FLINT_ASSERT(zl < zh);
    FLINT_ASSERT(zh <= zn);

    ztrunc = n_round_up(zn, BLK_SZ);

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

    P->stride = n_round_up(sd_fft_ctx_data_size(depth), 128);

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
int fft_small_plan_init_mpn(fft_small_plan_t P, mpn_ctx_t R,
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

    for (np = 4; np <= MPN_CTX_NCRTS; np++)
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
        ztrunc = n_round_up(zlen, BLK_SZ);
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

        _fft_small_plan_set_window(P, 0, zlen, zlen,
                    n_max(n_round_up(alen, BLK_SZ), n_round_up(blen, BLK_SZ)));

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

    _fft_small_plan_set_window(P, zl, zh, zn, xtrunc_max);

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

    FLINT_ASSERT(depth >= LG_BLK_SZ);
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

void fft_small_plan_clear(fft_small_plan_t P)
{
    if (P->use_direct_fft)
        sd_fft_ctx_clear(P->direct_fft);
}
