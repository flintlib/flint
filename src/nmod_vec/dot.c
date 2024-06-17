/*
    Copyright (C) 2011, 2021 Fredrik Johansson
    Copyright (C) 2024 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "machine_vectors.h"


dot_params_t _nmod_vec_dot_params(ulong len, nmod_t mod)
{
    ulong t2, t1, t0, u1, u0;

    dot_params_t params = {_DOT0, UWORD(0)};

    if ((mod.n & (mod.n - 1)) == 0)
        params.method = _DOT_POW2;
    else
    {
        umul_ppmm(t1, t0, mod.n - 1, mod.n - 1);
        umul_ppmm(t2, t1, t1, len);
        umul_ppmm(u1, u0, t0, len);
        add_ssaaaa(t2, t1, t2, t1, UWORD(0), u1);

        if (t2 != 0) // three limbs
        {
            /* we can accumulate 8 terms if n == mod.n is such that        */
            /*      8 * (n-1)**2 < 2**(2*FLINT_BITS), this is equivalent to    */
            /*      n <= ceil(sqrt(2**(2*FLINT_BITS-3)))                     */
#if (FLINT_BITS == 64)
            if (mod.n <= UWORD(6521908912666391107))
#else
            if (mod.n <= UWORD(1518500250))
#endif
                params.method = _DOT3_ACC;
            else
                params.method = _DOT3;
        }

        else if (t1 != 0) // two limbs
        {
#if (FLINT_BITS == 64)
            if (mod.n <= UWORD(1515531528) && len <= WORD(380368697))
            {
                // see end of file for these constraints; they imply 2 limbs
                params.method = _DOT2_32_SPLIT;
                NMOD_RED(params.pow2_precomp, (1L << DOT_SPLIT_BITS), mod);
            }
            else
#endif
            if (mod.n <= (UWORD(1) << (FLINT_BITS / 2)))
                params.method = _DOT2_HALF;
            else
                params.method = _DOT2;
        }

        // single limb
        else if (u0 != 0)
            params.method = _DOT1;

        // remaining case: u0 == 0
        // <=> mod.n == 1 or len == 0
        // => dot product is zero, _DOT0 already set
    }

    return params;
}

ulong
_nmod_vec_dot(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod, dot_params_t params)
{
    ulong res = UWORD(0);   /* covers _DOT0 */
    slong i;

#if defined(__AVX2__)
    if (params.method == _DOT1
        || (params.method == _DOT_POW2 && mod.n < (UWORD(1) << (FLINT_BITS / 2))))
    {
        vec4n dp = vec4n_zero();

        slong i = 0;
        for ( ; i+31 < len; i += 32)
        {
            dp = vec4n_add(dp, vec4n_mul(vec4n_load_unaligned(vec1+i+ 0), vec4n_load_unaligned(vec2+i+ 0)));
            dp = vec4n_add(dp, vec4n_mul(vec4n_load_unaligned(vec1+i+ 4), vec4n_load_unaligned(vec2+i+ 4)));
            dp = vec4n_add(dp, vec4n_mul(vec4n_load_unaligned(vec1+i+ 8), vec4n_load_unaligned(vec2+i+ 8)));
            dp = vec4n_add(dp, vec4n_mul(vec4n_load_unaligned(vec1+i+12), vec4n_load_unaligned(vec2+i+12)));
            dp = vec4n_add(dp, vec4n_mul(vec4n_load_unaligned(vec1+i+16), vec4n_load_unaligned(vec2+i+16)));
            dp = vec4n_add(dp, vec4n_mul(vec4n_load_unaligned(vec1+i+20), vec4n_load_unaligned(vec2+i+20)));
            dp = vec4n_add(dp, vec4n_mul(vec4n_load_unaligned(vec1+i+24), vec4n_load_unaligned(vec2+i+24)));
            dp = vec4n_add(dp, vec4n_mul(vec4n_load_unaligned(vec1+i+28), vec4n_load_unaligned(vec2+i+28)));
        }

        for ( ; i + 3 < len; i += 4)
            dp = vec4n_add(dp, vec4n_mul(vec4n_load_unaligned(vec1+i), vec4n_load_unaligned(vec2+i)));

        res = dp[0] + dp[1] + dp[2] + dp[3];

        for (; i < len; i++)
            res += vec1[i] * vec2[i];

        NMOD_RED(res, res, mod);
    }
    else if (params.method == _DOT_POW2)  // cannot use avx 32-bit mul
        _NMOD_VEC_DOT1(res, i, len, vec1[i], vec2[i], mod)
#else // if defined(__AVX2__)
    if (params.method == _DOT1 || params.method == _DOT_POW2)
        _NMOD_VEC_DOT1(res, i, len, vec1[i], vec2[i], mod)
#endif // if defined(__AVX2__)

    else if (params.method == _DOT2_32_SPLIT)
#if defined(__AVX2__)
    {
        const vec4n low_bits = vec4n_set_n(DOT_SPLIT_MASK);
        vec4n dp_lo = vec4n_zero();
        vec4n dp_hi = vec4n_zero();

        slong i = 0;
        for ( ; i+31 < len; i += 32)
        {
            dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(vec1+i+ 0), vec4n_load_unaligned(vec2+i+ 0)));
            dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(vec1+i+ 4), vec4n_load_unaligned(vec2+i+ 4)));
            dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(vec1+i+ 8), vec4n_load_unaligned(vec2+i+ 8)));
            dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(vec1+i+12), vec4n_load_unaligned(vec2+i+12)));
            dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(vec1+i+16), vec4n_load_unaligned(vec2+i+16)));
            dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(vec1+i+20), vec4n_load_unaligned(vec2+i+20)));
            dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(vec1+i+24), vec4n_load_unaligned(vec2+i+24)));
            dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(vec1+i+28), vec4n_load_unaligned(vec2+i+28)));

            dp_hi = vec4n_add(dp_hi, vec4n_bit_shift_right(dp_lo, DOT_SPLIT_BITS));
            dp_lo = vec4n_bit_and(dp_lo, low_bits);
        }

        for ( ; i + 3 < len; i += 4)
            dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(vec1+i), vec4n_load_unaligned(vec2+i)));

        dp_hi = vec4n_add(dp_hi, vec4n_bit_shift_right(dp_lo, DOT_SPLIT_BITS));
        dp_lo = vec4n_bit_and(dp_lo, low_bits);

        ulong hsum_lo = dp_lo[0] + dp_lo[1] + dp_lo[2] + dp_lo[3];
        const ulong hsum_hi = dp_hi[0] + dp_hi[1] + dp_hi[2] + dp_hi[3] + (hsum_lo >> DOT_SPLIT_BITS);
        hsum_lo &= DOT_SPLIT_MASK;

        for (; i < len; i++)
            hsum_lo += vec1[i] * vec2[i];

        NMOD_RED(res, params.pow2_precomp * hsum_hi + hsum_lo, mod);
    }
#else // if defined(__AVX2__)
        _NMOD_VEC_DOT2_32_SPLIT(res, i, len, vec1[i], vec2[i], mod, params.pow2_precomp)
#endif // if defined(__AVX2__)

    else if (params.method == _DOT2_HALF)
        _NMOD_VEC_DOT2_HALF(res, i, len, vec1[i], vec2[i], mod)
    else if (params.method == _DOT2)
        _NMOD_VEC_DOT2(res, i, len, vec1[i], vec2[i], mod)
    else if (params.method == _DOT3_ACC)
        _NMOD_VEC_DOT3_ACC(res, i, len, vec1[i], vec2[i], mod)
    else if (params.method == _DOT3)
        _NMOD_VEC_DOT3(res, i, len, vec1[i], vec2[i], mod)

    return res;
}

ulong
_nmod_vec_dot_ptr(nn_srcptr vec1, const nn_ptr * vec2, slong offset,
                            slong len, nmod_t mod, dot_params_t params)
{
    ulong res;
    slong i;
    NMOD_VEC_DOT(res, i, len, vec1[i], vec2[i][offset], mod, params);
    return res;
}

static ulong
nmod_fmma(ulong a, ulong b, ulong c, ulong d, nmod_t mod)
{
    a = nmod_mul(a, b, mod);
    NMOD_ADDMUL(a, c, d, mod);
    return a;
}

ulong
_nmod_vec_dot_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod, dot_params_t params)
{
    // FIXME advantage of this? (if yes, no similar advantage in non-rev version?)
    if (len <= 2 && params.method > _DOT1)
    {
        if (len == 2)
            return nmod_fmma(vec1[0], vec2[1], vec1[1], vec2[0], mod);
        if (len == 1)
            return nmod_mul(vec1[0], vec2[0], mod);
        return 0;
    }

    ulong res = UWORD(0);   /* covers _DOT0 */
    slong i;

#if defined(__AVX2__)
    if (params.method == _DOT1
        || (params.method == _DOT_POW2 && mod.n < (UWORD(1) << (FLINT_BITS / 2))))
    {
        vec4n dp = vec4n_zero();

        slong i = 0;
        for ( ; i+31 < len; i += 32)
        {
            const ulong ii = len - 32 - i; // >= 0
            dp = vec4n_add(dp, vec4n_mul(vec4n_load_unaligned(vec1+i+ 0), vec4n_permute_3_2_1_0(vec4n_load_unaligned(vec2+ii+28))));
            dp = vec4n_add(dp, vec4n_mul(vec4n_load_unaligned(vec1+i+ 4), vec4n_permute_3_2_1_0(vec4n_load_unaligned(vec2+ii+24))));
            dp = vec4n_add(dp, vec4n_mul(vec4n_load_unaligned(vec1+i+ 8), vec4n_permute_3_2_1_0(vec4n_load_unaligned(vec2+ii+20))));
            dp = vec4n_add(dp, vec4n_mul(vec4n_load_unaligned(vec1+i+12), vec4n_permute_3_2_1_0(vec4n_load_unaligned(vec2+ii+16))));
            dp = vec4n_add(dp, vec4n_mul(vec4n_load_unaligned(vec1+i+16), vec4n_permute_3_2_1_0(vec4n_load_unaligned(vec2+ii+12))));
            dp = vec4n_add(dp, vec4n_mul(vec4n_load_unaligned(vec1+i+20), vec4n_permute_3_2_1_0(vec4n_load_unaligned(vec2+ii+ 8))));
            dp = vec4n_add(dp, vec4n_mul(vec4n_load_unaligned(vec1+i+24), vec4n_permute_3_2_1_0(vec4n_load_unaligned(vec2+ii+ 4))));
            dp = vec4n_add(dp, vec4n_mul(vec4n_load_unaligned(vec1+i+28), vec4n_permute_3_2_1_0(vec4n_load_unaligned(vec2+ii+ 0))));
        }

        for ( ; i + 3 < len; i += 4)
            dp = vec4n_add(dp, vec4n_mul(vec4n_load_unaligned(vec1+i), vec4n_permute_3_2_1_0(vec4n_load_unaligned(vec2+len-4-i))));

        res = dp[0] + dp[1] + dp[2] + dp[3];

        for (; i < len; i++)
            res += vec1[i] * vec2[len-1-i];

        NMOD_RED(res, res, mod);
    }
    else if (params.method == _DOT_POW2)  // cannot use avx 32-bit mul
        _NMOD_VEC_DOT1(res, i, len, vec1[i], vec2[len-1-i], mod)
#else // if defined(__AVX2__)
    if (params.method == _DOT1 || params.method == _DOT_POW2)
        _NMOD_VEC_DOT1(res, i, len, vec1[i], vec2[len-1-i], mod)
#endif // if defined(__AVX2__)

    else if (params.method == _DOT2_32_SPLIT)
#if defined(__AVX2__)
    {
        const vec4n low_bits = vec4n_set_n(DOT_SPLIT_MASK);
        vec4n dp_lo = vec4n_zero();
        vec4n dp_hi = vec4n_zero();

        slong i = 0;
        for ( ; i+31 < len; i += 32)
        {
            const ulong ii = len - 32 - i; // >= 0
            dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(vec1+i+ 0), vec4n_permute_3_2_1_0(vec4n_load_unaligned(vec2+ii+28))));
            dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(vec1+i+ 4), vec4n_permute_3_2_1_0(vec4n_load_unaligned(vec2+ii+24))));
            dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(vec1+i+ 8), vec4n_permute_3_2_1_0(vec4n_load_unaligned(vec2+ii+20))));
            dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(vec1+i+12), vec4n_permute_3_2_1_0(vec4n_load_unaligned(vec2+ii+16))));
            dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(vec1+i+16), vec4n_permute_3_2_1_0(vec4n_load_unaligned(vec2+ii+12))));
            dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(vec1+i+20), vec4n_permute_3_2_1_0(vec4n_load_unaligned(vec2+ii+ 8))));
            dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(vec1+i+24), vec4n_permute_3_2_1_0(vec4n_load_unaligned(vec2+ii+ 4))));
            dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(vec1+i+28), vec4n_permute_3_2_1_0(vec4n_load_unaligned(vec2+ii+ 0))));

            dp_hi = vec4n_add(dp_hi, vec4n_bit_shift_right(dp_lo, DOT_SPLIT_BITS));
            dp_lo = vec4n_bit_and(dp_lo, low_bits);
        }

        for ( ; i + 3 < len; i += 4)
            dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(vec1+i), vec4n_permute_3_2_1_0(vec4n_load_unaligned(vec2+len-4-i))));

        dp_hi = vec4n_add(dp_hi, vec4n_bit_shift_right(dp_lo, DOT_SPLIT_BITS));
        dp_lo = vec4n_bit_and(dp_lo, low_bits);

        ulong hsum_lo = dp_lo[0] + dp_lo[1] + dp_lo[2] + dp_lo[3];
        const ulong hsum_hi = dp_hi[0] + dp_hi[1] + dp_hi[2] + dp_hi[3] + (hsum_lo >> DOT_SPLIT_BITS);
        hsum_lo &= DOT_SPLIT_MASK;

        for (; i < len; i++)
            hsum_lo += vec1[i] * vec2[len-1-i];

        NMOD_RED(res, params.pow2_precomp * hsum_hi + hsum_lo, mod);
    }
#else // if defined(__AVX2__)
        _NMOD_VEC_DOT2_32_SPLIT(res, i, len, vec1[i], vec2[len-1-i], mod, params.pow2_precomp)
#endif // if defined(__AVX2__)

    else if (params.method == _DOT2_HALF)
        _NMOD_VEC_DOT2_HALF(res, i, len, vec1[i], vec2[len-1-i], mod)
    else if (params.method == _DOT2)
        _NMOD_VEC_DOT2(res, i, len, vec1[i], vec2[len-1-i], mod)
    else if (params.method == _DOT3_ACC)
        _NMOD_VEC_DOT3_ACC(res, i, len, vec1[i], vec2[len-1-i], mod)
    else if (params.method == _DOT3)
        _NMOD_VEC_DOT3(res, i, len, vec1[i], vec2[len-1-i], mod)

    return res;
}




/*----------------------------------------*/
/* notes concerning the different methods */
/*----------------------------------------*/

// Why no vectorization in the general NMOD_VEC_DOT macro?
// attempts at vectorized versions (2024-06-16, for methods _DOT1,
// _DOT2_32_SPLIT) did not show an advantage except in "regular" cases where
// memory accesses are fast (typically, expr = v[i] or expr = v[len - 1 -i]).
// For these, there is dedicated code anyway.

// 2024-06-16 _DOT2_HALF is slightly faster than _DOT2
// 2024-06-16 _DOT3_ACC is slightly faster than _DOT3

/*---------------------------------------------*/
/* dot product for small modulus via splitting */
/*---------------------------------------------*/

// in short: with current DOT_SPLIT_BITS value 56,
// -> modulus n up to about 2**30.5
//       (more precisely, n <= 1515531528)
// -> length of dot product up to at least 380368697
//       (more precisely, len * (n-1)**3 < 2**120 + 2**56 - 2**112)

// APPROACH:
//
// Let n = mod.n, s = DOT_SPLIT_BITS
// As input, take pow2_precomp == 2**s % n
//
// -> avoiding modular reductions altogether, compute dp_lo and dp_hi such that
// the dot product without modular reduction is dp  =  dp_lo + 2**s * dp_hi
// -> finally, compute (dp_lo + pow2_precomp * dp_hi)  %  n
// -> done through repeating this: accumulate a few terms,
// move higher bits to dp_hi and keep lower ones in dp_lo

// PARAMETER CONSTRAINTS:
//
// 2024-06-16: currently, the code accumulates 8 terms as this showed slightly better performance
//
// -> constraint (C0-8):
// if we accumulate 8 terms (each a product of two integers reduced modulo n)
// on top of an s-bit integer, we require
//     2**s - 1 + 8 * (n-1)**2  <  2**64
// so one can take any modulus with
//     n <= 1 + floor(sqrt(2**61 - 2**(s-3)))
// in particular, n-1 < 2**30.5, (n-1)**2 < 2**61, (n-1)**3 < 2**91.5
//
// -> constraint (C0-4):
// similarly, if we accumulate 4 terms on top of an s-bit integer, we require
//     2**s - 1 + 4 * (n-1)**2  <  2**64
// so one can take any modulus with
//     n <= 1 + floor(sqrt(2**62 - 2**(s-2)))
// in particular, n-1 < 2**30.5, (n-1)**2 < 2**61, (n-1)**3 < 2**91.5
//
// -> constraint (C1):
// in the above representation of dp we will use a ulong for dp_hi,
// so we require      len * (n-1)**2 <= 2**s * (2**64 - 1)
// which is less restrictive than the below (C2)
//
// -> constraint (C2):
// for dp_lo + pow2_precomp * dp_hi to fit in a single word, we require
//      2**s - 1 + (n-1) dp_hi < 2**64.
// Since dp_hi <= len * (n-1)**2 / 2**s, it suffices to ensure
//     len * (n-1)**3 < 2**s * (2**64 + 1 - 2**s)
//
// sage: for s in range(40,64):
// ....:     nmax8 = 1 + floor(sqrt(2**61 - 2**(s-3)))             # (C0-8)
// ....:     nmax4 = 1 + floor(sqrt(2**62 - 2**(s-2)))             # (C0-4)
// ....:     lenmax4 = floor(2**s * (2**64 - 1) / (nmax4-1)**2)     # (C1)
// ....:     lenmax4_bis = ceil(2**s * (2**64 + 1 - 2**s) / (nmax4-1)**3) - 1       # (C2)
// ....:     lenmax8 = floor(2**s * (2**64 - 1) / (nmax8-1)**2)     # (C1)
// ....:     lenmax8_bis = ceil(2**s * (2**64 + 1 - 2**s) / (nmax8-1)**3) - 1       # (C2)
// ....:     print(f"{s}\t{nmax.nbits()}\t{nmax8}\t{lenmax8_bis}\t{nmax4}\t{lenmax4_bis}")
// ....:
// s       nbits   nmax8          (C2) for nmax8    nmax4           (C2) for nmax4
// 40      31      1518500205      5792             2147483584      2048
// 41      31      1518500160      11585            2147483520      4096
// 42      31      1518500069      23170            2147483392      8192
// 43      31      1518499888      46340            2147483136      16384
// 44      31      1518499526      92681            2147482624      32768
// 45      31      1518498802      185363           2147481600      65536
// 46      31      1518497354      370728           2147479552      131072
// 47      31      1518494458      741458           2147475456      262145
// 48      31      1518488665      1482921          2147467264      524292
// 49      31      1518477080      2965866          2147450880      1048592
// 50      31      1518453909      5931822          2147418111      2097216
// 51      31      1518407566      11864007         2147352572      4194560
// 52      31      1518314875      23729463         2147221488      8389632
// 53      31      1518129478      47464722         2146959296      16781313
// 54      31      1517758614      94952640         2146434816      33570828
// 55      31      1517016615      189998167        2145385471      67174496
// 56      31      1515531528      380368697        2143285240      134480642
// 57      31      1512556978      762233438        2139078592      269490216
// 58      31      1506590261      1530504392       2130640379      541115017
// 59      31      1494585366      3085595597       2113662895      1090922784
// 60      31      1470281545      6273201268       2079292102      2217911575
// 61      31      1420426920      12986760413      2008787014      4591513178
// 62      31      1315059793      28054608908      1859775394      9918802104
// 63      31      1073741825      68719476736      1518500250      24296004047

