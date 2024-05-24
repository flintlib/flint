/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdint.h>
#include "n_mod.h"
#include "n_mod_vec.h"

#if FLINT64
typedef const uint32_t * UH_srcptr;
#else
typedef const uint16_t * UH_srcptr;
#endif

/* ctx->nu <= 2^{FLINT_BITS / 4}, assumes that result fits inside a limb */
ulong
_n_mod_vec_dot_0(nn_srcptr ap, nn_srcptr bp, slong len, n_mod_ctx_srcptr ctx)
{
    ulong res = 0;
    slong ix;
    UH_srcptr _ap = (UH_srcptr) ap, _bp = (UH_srcptr) bp;

    if (len <= WORD(0))
        FLINT_UNREACHABLE;

#if defined(__GNUC__)
# pragma GCC unroll 4
#endif
    for (ix = 0; ix < 2 * len; ix++)
        res += _ap[ix] * _bp[ix];

    return n_mod_set_ui(res, ctx);
}

/* ctx->nu <= 2^{FLINT_BITS / 2} */
ulong
_n_mod_vec_dot_1(nn_srcptr ap, nn_srcptr bp, slong len, n_mod_ctx_srcptr ctx)
{
    ulong r0 = 0, r1 = 0;
    slong ix;

    if (len <= WORD(0))
        FLINT_UNREACHABLE;

#if defined(__GNUC__)
# pragma GCC unroll 4
#endif
    for (ix = 0; ix < len; ix++)
    {
#if defined(__GNUC__)
        r1 += __builtin_add_overflow(r0, ap[ix] * bp[ix], &r0);
#else
        ulong prod = ap[ix] * bp[ix];
        add_ssaaaa(r1, r0, r1, r0, 0, prod);
#endif
    }

    return n_mod_set_uiui(r1, r0, ctx);
}

/* Result fits in 2 limbs */
ulong
_n_mod_vec_dot_2(nn_srcptr ap, nn_srcptr bp, slong len, n_mod_ctx_srcptr ctx)
{
#if defined(__GNUC__) && FLINT64
    __uint128_t res = 0;
#else
    ulong r0 = 0, r1 = 0;
#endif
    slong ix;

    if (len <= WORD(0))
        FLINT_UNREACHABLE;

#if defined(__GNUC__)
# pragma GCC unroll 4
#endif
    for (ix = 0; ix < len; ix++)
    {
#if defined(__GNUC__) && FLINT64
        __uint128_t a = ap[ix], b = bp[ix];
        res += a * b;
#else
        ulong t0, t1;
        umul_ppmm(t1, t0, ap[ix], bp[ix]);
        add_ssaaaa(r1, r0, r1, r0, t1, t0);
#endif
    }

#if defined(__GNUC__) && FLINT64
    return n_mod_set_uiui((ulong) (res >> FLINT_BITS), (ulong) res, ctx);
#else
    return n_mod_set_uiui(r1, r0, ctx);
#endif
}

/* Result fits in 3 limbs */
ulong
_n_mod_vec_dot_3(nn_srcptr ap, nn_srcptr bp, slong len, n_mod_ctx_srcptr ctx)
{
#if defined(__GNUC__) && FLINT64
    __uint128_t res = 0;
#else
    ulong r0 = 0, r1 = 0;
#endif
    slong ix;

    if (len <= WORD(0))
        FLINT_UNREACHABLE;

#if defined(__GNUC__)
# pragma GCC unroll 4
#endif
    for (ix = 0; ix < len; ix++)
    {
#if defined(__GNUC__) && FLINT64
        __uint128_t a = ap[ix], b = bp[ix];
        res += a * b;
#else
        ulong t0, t1;
        umul_ppmm(t1, t0, ap[ix], bp[ix]);
        add_ssaaaa(r1, r0, r1, r0, t1, t0);
#endif
    }

#if defined(__GNUC__) && FLINT64
    return n_mod_set_uiuiui((ulong) (res >> FLINT_BITS), (ulong) res, ctx);
#else
    return n_mod_set_uiuiui(r1, r0, ctx);
#endif
}
