/*
    Copyright (C) 2025 Fredrik Johansson
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Exponentiation modulo a word: shared by nmod/pow.c (nmod_pow_ui) and
   ulong_extras/powmod2_ui_preinv.c (n_powmod2_ui_preinv). */

#ifndef NMOD_POW_IMPL_H
#define NMOD_POW_IMPL_H

#include "nmod.h"

/*
    The algorithms are static functions, called both by the exported
    functions of nmod/pow.c and by nmod_pow_ui and n_powmod2_ui_preinv
    through _pow_ui: with a shared library, a call to an exported function
    would go through the PLT, which costs up to a third of the time for
    small exponents. They are not inlined into _pow_ui, which keeps its
    fast paths short, and they take the fields of the modulus as separate
    arguments, in registers (an nmod_t argument is passed in memory and
    would be copied for the call). All require exp >= 1.
*/

FLINT_STATIC_NOINLINE ulong
_redc_pow_ui(ulong a, ulong exp, const nmod_redc_ctx_t ctx)
{
    ulong x;

    while ((exp & 1) == 0)
    {
        a = nmod_redc_mul(a, a, ctx);
        exp >>= 1;
    }

    x = a;
    while (exp >>= 1)
    {
        a = nmod_redc_mul(a, a, ctx);
        if (exp & 1)
            x = nmod_redc_mul(x, a, ctx);
    }

    return x;
}

FLINT_STATIC_NOINLINE ulong
_redc_fast_pow_ui(ulong a, ulong exp, const nmod_redc_ctx_t ctx)
{
    ulong x;

    while ((exp & 1) == 0)
    {
        a = nmod_redc_fast_mul(a, a, ctx);
        exp >>= 1;
    }

    x = a;
    while (exp >>= 1)
    {
        a = nmod_redc_fast_mul(a, a, ctx);
        if (exp & 1)
            x = nmod_redc_fast_mul(x, a, ctx);
    }

    return x;
}

FLINT_STATIC_NOINLINE ulong
_pow_ui_binexp(ulong a, ulong exp, ulong n, ulong ninv, flint_bitcnt_t norm)
{
    ulong x;

    a <<= norm;
    n <<= norm;

    while ((exp & 1) == 0)
    {
        a = n_mulmod_preinv(a, a, n, ninv, norm);
        exp >>= 1;
    }

    x = a;

    while (exp >>= 1)
    {
        a = n_mulmod_preinv(a, a, n, ninv, norm);

        if (exp & 1)
            x = n_mulmod_preinv(x, a, n, ninv, norm);
    }

    return x >> norm;
}

/* n <= 2^(FLINT_BITS/2) */
FLINT_STATIC_NOINLINE ulong
_pow_ui_binexp_halflimb(ulong a, ulong exp, ulong n, ulong ninv, flint_bitcnt_t norm)
{
    nmod_t mod;
    ulong x;

    mod.n = n;
    mod.ninv = ninv;
    mod.norm = norm;

    FLINT_ASSERT(NMOD_IS_HALFLIMB(mod));

    while ((exp & 1) == 0)
    {
        a = _nmod_mul_halflimb(a, a, mod);
        exp >>= 1;
    }

    x = a;

    while (exp >>= 1)
    {
        a = _nmod_mul_halflimb(a, a, mod);

        if (exp & 1)
            x = _nmod_mul_halflimb(x, a, mod);
    }

    return x;
}

/* odd n */
FLINT_STATIC_NOINLINE ulong
_pow_ui_redc(ulong a, ulong exp, ulong n, ulong ninv, flint_bitcnt_t norm)
{
    nmod_t mod;
    ulong x;
    nmod_redc_ctx_t ctx;

    mod.n = n;
    mod.ninv = ninv;
    mod.norm = norm;

    nmod_redc_ctx_init_nmod(ctx, mod);
    a = nmod_redc_set_nmod(a, ctx);

    if (nmod_redc_can_use_fast(ctx))
        x = _redc_fast_pow_ui(a, exp, ctx);
    else
        x = _redc_pow_ui(a, exp, ctx);

    return nmod_redc_get_nmod(x, ctx);
}

/* odd n < 2^(FLINT_BITS/2 - 2): half-limb Montgomery form, values in [0, 2n) */
FLINT_STATIC_NOINLINE ulong
_pow_ui_redc_half(ulong a, ulong exp, ulong n, ulong ninv, flint_bitcnt_t norm)
{
    nmod_t mod;
    ulong x;
    nmod_redc_ctx_t ctx;

    mod.n = n;
    mod.ninv = ninv;
    mod.norm = norm;

    nmod_redc_half_ctx_init_nmod(ctx, mod);
    FLINT_ASSERT(nmod_redc_half_can_use_fast(ctx));
    a = nmod_redc_half_set_nmod(a, ctx);

    while ((exp & 1) == 0)
    {
        a = nmod_redc_half_fast_mul(a, a, ctx);
        exp >>= 1;
    }

    x = a;
    while (exp >>= 1)
    {
        a = nmod_redc_half_fast_mul(a, a, ctx);
        if (exp & 1)
            x = nmod_redc_half_fast_mul(x, a, ctx);
    }

    return nmod_redc_half_get_nmod(x, ctx);
}

/*
    Exponents from which Montgomery form, including the conversions, is
    faster than binary exponentiation in the standard representation, for
    odd n (exponents up to 4 use the inline code below in any case):
    - 2^(FLINT_BITS/2) < n < 2^(FLINT_BITS - 2): from NMOD_POW_UI_REDC_CUTOFF;
    - n >= 2^(FLINT_BITS - 2), where Montgomery multiplication needs a final
      correction: from NMOD_POW_UI_REDC_FULLWORD_CUTOFF;
    - n <= 2^(FLINT_BITS/2): binary exponentiation with the one-limb
      product _nmod_mul_halflimb stays faster up to larger exponents.
      Half-limb Montgomery form (n < 2^(FLINT_BITS/2 - 2)) pays off from
      NMOD_POW_UI_REDC_HALF_CUTOFF, full Montgomery form (other n) from
      NMOD_POW_UI_REDC_HALFLIMB_CUTOFF.
    The crossovers differ between machines (Cascade Lake, Ice Lake, Zen 4,
    Apple M4); the cutoffs are chosen so that no machine loses against
    binary exponentiation, at the price of smaller gains on some of them.
*/
#define NMOD_POW_UI_REDC_CUTOFF 8
#define NMOD_POW_UI_REDC_FULLWORD_CUTOFF 16
#define NMOD_POW_UI_REDC_HALF_CUTOFF 1024
#define NMOD_POW_UI_REDC_HALFLIMB_CUTOFF 2048

FLINT_FORCE_INLINE ulong
_pow_ui(ulong a, ulong exp, nmod_t mod)
{
    if (exp <= 1)
        return (exp == 1) ? a : (mod.n != 1);

    if (NMOD_IS_HALFLIMB(mod))
    {
        if (exp <= 4)
        {
            ulong x = _nmod_mul_halflimb(a, a, mod);

            if (exp == 2)
                return x;
            else if (exp == 3)
                return _nmod_mul_halflimb(a, x, mod);
            else
                return _nmod_mul_halflimb(x, x, mod);
        }

        if (mod.n % 2 == 1 && mod.n > 1)
        {
            if (exp >= NMOD_POW_UI_REDC_HALF_CUTOFF
                    && mod.norm >= FLINT_BITS / 2 + 2)
                return _pow_ui_redc_half(a, exp, mod.n, mod.ninv, mod.norm);

            if (exp >= NMOD_POW_UI_REDC_HALFLIMB_CUTOFF)
                return _pow_ui_redc(a, exp, mod.n, mod.ninv, mod.norm);
        }

        return _pow_ui_binexp_halflimb(a, exp, mod.n, mod.ninv, mod.norm);
    }

    if (exp <= 4)
    {
        ulong x = nmod_mul(a, a, mod);

        if (exp == 2)
            return x;
        else if (exp == 3)
            return nmod_mul(a, x, mod);
        else
            return nmod_mul(x, x, mod);
    }

    if (mod.n % 2 == 1 && exp >= ((mod.norm >= 2) ? NMOD_POW_UI_REDC_CUTOFF
                                                   : NMOD_POW_UI_REDC_FULLWORD_CUTOFF))
        return _pow_ui_redc(a, exp, mod.n, mod.ninv, mod.norm);

    return _pow_ui_binexp(a, exp, mod.n, mod.ninv, mod.norm);
}

#endif
