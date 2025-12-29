/*
    Copyright (C) 2008 Peter Shrimpton
    Copyright (C) 2009 William Hart
    Copyright (C) 2014, 2015 Dana Jacobsen
    Copyright (C) 2015 Kushagra Singh
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdint.h>
#include "ulong_extras.h"
#include "nmod.h"

#define SMALL_ODDPRIME_LIMIT 32768
#define WITNESS_BASE_HASH_SIZE 98304
#define NUM_OVERSIZE_BASES 4903

/* To keep this file readable, the lookup tables have been
   placed in a seprate header file. */
#include "is_prime_tables.h"

/* Branch-free hash table lookup */
static int u32_is_base2_pseudoprime(uint32_t x)
{
    uint32_t c = base2_32bit_pseudoprimes_hash_keys[x >> 25];
    uint32_t h = (x ^ c) % 2560;

    /* Linear probing */
    return (base2_32bit_pseudoprimes[h] == x) | (base2_32bit_pseudoprimes[h + 1] == x);
}

/* Faster arithmetic for 32-bit probable prime test on 64-bit machines. */
#if FLINT_BITS == 64

static uint32_t
u32_addmod(uint32_t x, uint32_t y, uint32_t n)
{
    return (n - y > x ? x + y : x + y - n);
}

static uint32_t
u32_rem_u64_shoup(uint64_t x, uint32_t n, uint64_t ninv)
{
    return n_mulmod_shoup(1, x, ninv, n);
}

static uint32_t
u32_2_powmod(uint32_t exp, uint32_t n, uint64_t ninv)
{
    uint32_t x, bit;
    unsigned int ebits;

    if (exp < FLINT_BITS)
        return u32_rem_u64_shoup(UWORD(1) << exp, n, ninv);

    ebits = FLINT_BITS - flint_clz(exp);
    bit = UWORD(1) << (ebits - 6);

    x = u32_rem_u64_shoup(((uint64_t) 1) << (exp >> (ebits - 6)), n, ninv);

    while (bit >>= 1)
    {
        x = u32_rem_u64_shoup((uint64_t) x * x, n, ninv);

        if (bit & exp)
            x = u32_addmod(x, x, n);
    }

    return x;
}

static int u32_is_base2_probabprime(uint32_t n)
{
    uint32_t d = n-1, s = 0;
    uint64_t ninv = n_mulmod_precomp_shoup(UWORD(1), n);

    s = flint_ctz(d);
    d >>= s;

    uint64_t t = d;
    uint64_t y;

    y = u32_2_powmod(t, n, ninv);
    if (y == 1)
        return 1;

    t <<= 1;
    while ((t != n - 1) && (y != n - 1))
    {
        y = u32_rem_u64_shoup(y * y, n, ninv);
        t <<= 1;
    }

    return y == n - 1;
}

#else

static int u32_is_base2_probabprime(uint32_t n)
{
    return n_is_strong_probabprime2_preinv(n, n_preinvert_limb(n), 2, (n - 1) >> flint_ctz(n - 1));
}

#endif

static uint32_t get_witness_base(uint64_t n)
{
    uint32_t hash, b;

    /* The specific hash function used to generate the table. */
    hash = ((uint32_t) ((n * 314159265) >> 15) % WITNESS_BASE_HASH_SIZE);

    /* Read 3 bytes = 24 bits. */
    b = witness_base_tab[3 * hash] | (witness_base_tab[3 * hash + 1] << 8) | (witness_base_tab[3 * hash + 2] << 16);

    /* A small b value encodes an index into the oversize table for bases > 24 bits. */
    if (b < NUM_OVERSIZE_BASES)
        b = oversize[b];

    return b;
}

static int
redc_fast_equal_redc(ulong a, ulong b, nmod_redc_ctx_t ctx)
{
    return (a == b) || (a == b + ctx->mod.n);
}

/* Assumes that a < n. In particular, this is fine with 32-bit bases
   when used for n > 2^32. */
static int
_n_is_strong_probabprime_redc(ulong a, ulong d, ulong one_red, nmod_redc_ctx_t ctx)
{
    ulong t = d, n = ctx->mod.n, y, minus_one_red = n - one_red;

    FLINT_ASSERT(a < n);
    FLINT_ASSERT(a >= 2);
    FLINT_ASSERT(a != n - 1);

    if (nmod_redc_can_use_fast(ctx))
    {
        if (a == 2)
            y = _nmod_redc_fast_2_pow_ui(t, ctx);
        else
            y = _nmod_redc_fast_pow_ui(nmod_redc_set_nmod(a, ctx), t, ctx);

        if (redc_fast_equal_redc(y, one_red, ctx))
            return 1;
        t <<= 1;

        while ((t != n - 1) && !redc_fast_equal_redc(y, minus_one_red, ctx))
        {
            y = nmod_redc_fast_mul(y, y, ctx);
            t <<= 1;
        }

        return redc_fast_equal_redc(y, minus_one_red, ctx);
    }
    else
    {
        if (a == 2)
            y = _nmod_redc_2_pow_ui(t, ctx);
        else
            y = _nmod_redc_pow_ui(nmod_redc_set_nmod(a, ctx), t, ctx);

        if (y == one_red)
            return 1;
        t <<= 1;

        while ((t != n - 1) && (y != minus_one_red))
        {
            y = nmod_redc_mul(y, y, ctx);
            t <<= 1;
        }

        return y == minus_one_red;
    }
}

static int n_is_oddprime_small2(ulong n)
{
    ulong q = n / 2;
    ulong x = (q & 63);
    return (flint_odd_prime_lookup[q / 64] & (((uint64_t) 1) << x)) >> x;
}

int
n_is_prime_odd_no_trial(ulong n)
{
    FLINT_ASSERT(n % 2 != 0);

    if (n <= UINT32_MAX)
    {
        if (n < SMALL_ODDPRIME_LIMIT)
            return n_is_oddprime_small2(n);

        return u32_is_base2_probabprime(n) && !u32_is_base2_pseudoprime(n);
    }
    else
    {
        ulong d, norm, one_red;
        nmod_redc_ctx_t ctx;

        nmod_redc_ctx_init_ui(ctx, n);
        one_red = nmod_redc_set_ui(1, ctx);

        d = n - 1;
        norm = flint_ctz(d);
        d >>= norm;

        if (!_n_is_strong_probabprime_redc(2, d, one_red, ctx))
            return 0;

        return _n_is_strong_probabprime_redc(get_witness_base(n), d, one_red, ctx);
    }
}

int
n_is_prime(ulong n)
{
    if (n < SMALL_ODDPRIME_LIMIT)
        return (n % 2 ? n_is_oddprime_small2(n) : n == 2);

    if (!(n % 2) || !(n % 3) || !(n % 5) || !(n % 7))       return 0;

    if (!(n % 11) || !(n % 13) || !(n % 17) || !(n % 19) ||
        !(n % 23) || !(n % 29) || !(n % 31) || !(n % 37) ||
        !(n % 41) || !(n % 43) || !(n % 47) || !(n % 53))   return 0;

    if (n > UINT32_MAX &&
        (!(n % 59) || !(n % 61) || !(n % 67) || !(n % 71) || !(n % 73) ||
         !(n % 79) || !(n % 83) || !(n % 89) || !(n % 97) || !(n %101) ||
         !(n %103) || !(n % 107) || !(n %109) || !(n %113) || !(n % 127) ||
         !(n %131) || !(n % 137) || !(n %139) || !(n %149)))
        return 0;

    return n_is_prime_odd_no_trial(n);
}

