/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdint.h>
#include "nmod.h"
#include "nmod/pow_impl.h"

ulong
_nmod_redc_pow_ui(ulong a, ulong exp, const nmod_redc_ctx_t ctx)
{
    return _redc_pow_ui(a, exp, ctx);
}

ulong
_nmod_redc_fast_pow_ui(ulong a, ulong exp, const nmod_redc_ctx_t ctx)
{
    return _redc_fast_pow_ui(a, exp, ctx);
}

ulong
_nmod_pow_ui_binexp(ulong a, ulong exp, nmod_t mod)
{
    return _pow_ui_binexp(a, exp, mod.n, mod.ninv, mod.norm);
}

ulong
_nmod_pow_ui_binexp_halflimb(ulong a, ulong exp, nmod_t mod)
{
    return _pow_ui_binexp_halflimb(a, exp, mod.n, mod.ninv, mod.norm);
}

ulong
_nmod_pow_ui_redc(ulong a, ulong exp, nmod_t mod)
{
    return _pow_ui_redc(a, exp, mod.n, mod.ninv, mod.norm);
}

ulong
_nmod_pow_ui_redc_half(ulong a, ulong exp, nmod_t mod)
{
    return _pow_ui_redc_half(a, exp, mod.n, mod.ninv, mod.norm);
}

ulong
nmod_pow_ui(ulong a, ulong exp, nmod_t mod)
{
    FLINT_ASSERT(a < mod.n);

    return _pow_ui(a, exp, mod);
}

/* Powering with base = 2 */

#if FLINT_BITS == 64
#define LG_FLINT_BITS 6
#else
#define LG_FLINT_BITS 5
#endif

ulong
_nmod_2_pow_ui_binexp(ulong exp, nmod_t mod)
{
    ulong x, bit;
    unsigned int ebits;

    if (exp < FLINT_BITS)
        return nmod_set_ui(UWORD(1) << exp, mod);

    ebits = FLINT_BITS - flint_clz(exp);
    bit = UWORD(1) << (ebits - LG_FLINT_BITS);
    x = nmod_set_ui(UWORD(1) << (exp >> (ebits - LG_FLINT_BITS)), mod);

    if (mod.norm == 0)
    {
        while (bit >>= 1)
        {
            x = _nmod_mul_fullword(x, x, mod);
            if (bit & exp)
                x = nmod_add(x, x, mod);
        }
    }
    else
    {
        while (bit >>= 1)
        {
            x = nmod_mul(x, x, mod);
            if (bit & exp)
                x = nmod_add(x, x, mod);
        }
    }

    return x;
}

ulong
_nmod_redc_2_pow_ui(ulong exp, const nmod_redc_ctx_t ctx)
{
    ulong x;
    ulong bit;
    unsigned int ebits;

    if (exp < FLINT_BITS)
        return nmod_redc_set_ui(UWORD(1) << exp, ctx);

    ebits = FLINT_BITS - flint_clz(exp);
    bit = UWORD(1) << (ebits - LG_FLINT_BITS);
    x = UWORD(1) << (exp >> (ebits - LG_FLINT_BITS));
    x = nmod_redc_set_ui(x, ctx);

    while (bit >>= 1)
    {
        x = nmod_redc_mul(x, x, ctx);
        if (bit & exp)
            x = nmod_redc_add(x, x, ctx);
    }

    return x;
}

ulong
_nmod_redc_fast_2_pow_ui(ulong exp, const nmod_redc_ctx_t ctx)
{
    ulong x;
    ulong bit;
    unsigned int ebits;

    if (exp < FLINT_BITS)
        return nmod_redc_set_ui(UWORD(1) << exp, ctx);

    ebits = FLINT_BITS - flint_clz(exp);
    bit = UWORD(1) << (ebits - LG_FLINT_BITS);
    x = UWORD(1) << (exp >> (ebits - LG_FLINT_BITS));
    x = nmod_redc_set_ui(x, ctx);

    while (bit >>= 1)
    {
        x = nmod_redc_fast_mul(x, x, ctx);
        if (bit & exp)
            x = nmod_redc_fast_mul_two(x, ctx);
    }

    return x;
}

ulong nmod_2_pow_ui(ulong exp, nmod_t mod)
{
    ulong x;
    nmod_redc_ctx_t ctx;

    /* Montgomery arithmetic requires an odd modulus */
    if (exp < (UWORD(1) << 20) || mod.n % 2 == 0)
    {
        return _nmod_2_pow_ui_binexp(exp, mod);
    }
    else
    {
        nmod_redc_ctx_init_nmod(ctx, mod);

        if (nmod_redc_can_use_fast(ctx))
            x = _nmod_redc_fast_2_pow_ui(exp, ctx);
        else
            x = _nmod_redc_2_pow_ui(exp, ctx);

        return nmod_redc_get_nmod(x, ctx);
    }
}

