/*
    Copyright (C) 2024, 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Selection of the dot product method. This is not an inline function in
    nmod_vec.h: the choice depends on the kernels the library was compiled
    with (AVX2, AVX-512, IFMA), which a calling translation unit compiled
    with other flags cannot see. The split methods _DOT2_SPLIT and
    _DOT2_HALF have no length bound (their final step is a two-word
    reduction, see _NMOD_VEC_DOT_SPLIT_COMBINE in nmod_vec.h).
*/

#include "nmod.h"
#include "nmod_vec.h"

#define _FIXED_LEN_MOD_BOUNDS(fixedlen, onelimb_bnd, twolimb_bnd) \
        if (len == fixedlen)                                      \
        {                                                         \
            if (mod.n <= UWORD(onelimb_bnd))                      \
            {                                                     \
                dot_params_t params = {_DOT1, UWORD(0)};          \
                return params;                                    \
            }                                                     \
            if (mod.n <= UWORD(twolimb_bnd))                      \
            {                                                     \
                dot_params_t params = {_DOT2, UWORD(0)};          \
                return params;                                    \
            }                                                     \
            dot_params_t params = {_DOT3, UWORD(0)};              \
            return params;                                        \
        }

dot_params_t _nmod_vec_dot_params(ulong len, nmod_t mod)
{
    if (len == 0 || mod.n == 1)
    {
        dot_params_t params = {_DOT0, UWORD(0)};
        return params;
    }
    if ((mod.n & (mod.n - 1)) == 0)
    {
        dot_params_t params = {_DOT_POW2, UWORD(0)};
        return params;
    }
    // from here on len >= 1, n > 1 not power of 2

    // short dot products: we use only _DOT1, _DOT2, _DOT3 in that case
    if (len <= 11)
    {
#if FLINT_BITS == 64
        // 64 bits:  k limbs  <=>  n <= ceil(2**(32*k) / sqrt(len))
        _FIXED_LEN_MOD_BOUNDS(11, 1294981365,  5561902608746059656);
        _FIXED_LEN_MOD_BOUNDS(10, 1358187914,  5833372668713515885);
        _FIXED_LEN_MOD_BOUNDS( 9, 1431655766,  6148914691236517206);
        _FIXED_LEN_MOD_BOUNDS( 8, 1518500250,  6521908912666391107);
        _FIXED_LEN_MOD_BOUNDS( 7, 1623345051,  6972213902555716131);
        _FIXED_LEN_MOD_BOUNDS( 6, 1753413057,  7530851732716320753);
        _FIXED_LEN_MOD_BOUNDS( 5, 1920767767,  8249634742471189718);
        _FIXED_LEN_MOD_BOUNDS( 4, 2147483648,  9223372036854775808);
        _FIXED_LEN_MOD_BOUNDS( 3, 2479700525, 10650232656628343402);
        _FIXED_LEN_MOD_BOUNDS( 2, 3037000500, 13043817825332782213);
#else  // FLINT_BITS == 64
        // 32 bits: k limbs  <=>  n <= ceil(2**(16*k) / sqrt(len))
        _FIXED_LEN_MOD_BOUNDS(11, 19760, 1294981365);
        _FIXED_LEN_MOD_BOUNDS(10, 20725, 1358187914);
        _FIXED_LEN_MOD_BOUNDS( 9, 21846, 1431655766);
        _FIXED_LEN_MOD_BOUNDS( 8, 23171, 1518500250);
        _FIXED_LEN_MOD_BOUNDS( 7, 24771, 1623345051);
        _FIXED_LEN_MOD_BOUNDS( 6, 26755, 1753413057);
        _FIXED_LEN_MOD_BOUNDS( 5, 29309, 1920767767);
        _FIXED_LEN_MOD_BOUNDS( 4, 32768, 2147483648);
        _FIXED_LEN_MOD_BOUNDS( 3, 37838, 2479700525);
        _FIXED_LEN_MOD_BOUNDS( 2, 46341, 3037000500);
#endif  // FLINT_BITS == 64
        // remains len == 1
        if (mod.n <= (UWORD(1) << FLINT_BITS / 2))
        {
            dot_params_t params = {_DOT1, UWORD(0)};
            return params;
        }
        dot_params_t params = {_DOT2, UWORD(0)};
        return params;
    }

    if (mod.n <= UWORD(1) << (FLINT_BITS / 2)) // implies <= 2 limbs
    {
        const ulong t0 = (mod.n - 1) * (mod.n - 1);
        ulong u1 = n_mulhi(t0, len);
        if (u1 == 0)  // 1 limb
        {
            dot_params_t params = {_DOT1, UWORD(0)};
            return params;
        }

        // u1 != 0 <=> 2 limbs
#if defined(NMOD_VEC_DOT2_SPLIT_MAX_LEN)
        if (len >= NMOD_VEC_DOT2_SPLIT_MAX_LEN)
        {
            dot_params_t params = {_DOT2, UWORD(0)};
            return params;
        }
#endif
#if (FLINT_BITS == 64) // _SPLIT: see the end of dot.c for this constraint
        if (mod.n <= UWORD(1515531528))
        {
            ulong pow2_precomp;
            NMOD_RED(pow2_precomp, (UWORD(1) << DOT_SPLIT_BITS), mod);
            dot_params_t params = {_DOT2_SPLIT, pow2_precomp};
            return params;
        }
#endif // FLINT_BITS == 64
#if (FLINT_BITS == 64) && defined(__AVX2__)
        ulong pow2_precomp;
        NMOD_RED(pow2_precomp, (UWORD(1) << DOT_SPLIT_BITS), mod);
        dot_params_t params = {_DOT2_HALF, pow2_precomp};
        return params;
#elif defined(NMOD_VEC_DOT2_SPLIT_MAX_LEN)
        dot_params_t params = {_DOT2, UWORD(0)};
        return params;
#else // (FLINT_BITS == 64) && defined(__AVX2__)
        dot_params_t params = {_DOT2_HALF, UWORD(0)};
        return params;
#endif // (FLINT_BITS == 64) && defined(__AVX2__)
    }
    // from here on, mod.n > 2**(FLINT_BITS / 2)
    // --> unreduced dot cannot fit in 1 limb

    ulong t2, t1, t0, u1;
    umul_ppmm(t1, t0, mod.n - 1, mod.n - 1);
    umul_ppmm(t2, t1, t1, len);
    u1 = n_mulhi(t0, len);
    add_ssaaaa(t2, t1, t2, t1, UWORD(0), u1);

    if (t2 == 0) // 2 limbs
    {
#if NMOD_VEC_HAVE_DOT_U52
        if (mod.n <= (UWORD(1) << 52) && len >= NMOD_VEC_DOT_U52_MIN_LEN)
        {
            dot_params_t params = {_DOT_U52, UWORD(0)};
            return params;
        }
#endif
#if NMOD_VEC_DOT_SPLIT_LIMBS_ENABLED
        // up to NMOD_VEC_DOT_SPLIT_LIMBS_MAX_BITS bits (with IFMA, _DOT_U64
        // is faster beyond)
        if (len >= NMOD_VEC_DOT_SPLIT_LIMBS_MIN_LEN
                && mod.n <= NMOD_VEC_DOT_SPLIT_LIMBS_MAX_N)
        {
            dot_params_t params = {_DOT_SPLIT_LIMBS, UWORD(0)};
            return params;
        }
#endif
#if NMOD_VEC_HAVE_DOT_U64
        if (mod.n > (UWORD(1) << 52) && len >= NMOD_VEC_DOT_U64_MIN_LEN)
        {
            dot_params_t params = {_DOT_U64, UWORD(0)};
            return params;
        }
#endif
        dot_params_t params = {_DOT2, UWORD(0)};
        return params;
    }

    // 3 limbs:
#if NMOD_VEC_DOT_SPLIT_LIMBS_ENABLED
    if (len >= NMOD_VEC_DOT_SPLIT_LIMBS_MIN_LEN
            && mod.n <= NMOD_VEC_DOT_SPLIT_LIMBS_MAX_N)
    {
        dot_params_t params = {_DOT3_SPLIT_LIMBS, UWORD(0)};
        return params;
    }
#endif
#if NMOD_VEC_HAVE_DOT_U64
    if (len >= NMOD_VEC_DOT3_U64_MIN_LEN)
    {
        dot_params_t params = {_DOT3_U64, UWORD(0)};
        return params;
    }
#endif
#if (FLINT_BITS == 64)
    if (mod.n <= UWORD(6521908912666391107))  // room for accumulating 8 terms
#else
    if (mod.n <= UWORD(1518500250))           // room for accumulating 8 terms
#endif
    {
        dot_params_t params = {_DOT3_ACC, UWORD(0)};
        return params;
    }

    dot_params_t params = {_DOT3, UWORD(0)};
    return params;
}

#undef _FIXED_LEN_MOD_BOUNDS
