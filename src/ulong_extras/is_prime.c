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
#define NUM_32BIT_PSEUDOPRIMES 2314
#define NUM_OVERSIZE_BASES 221
#define WITNESS_BASE_HASH_SIZE 131072
#define WITNESS_BASE_COMPRESSED_LEN 43691

/* To keep this file readable, the lookup tables have been
   placed in a seprate header file. */
#include "is_prime_tables.h"

static int u32_is_base2_pseudoprime(uint32_t x)
{
    const uint32_t * arr = base2_32bit_pseudoprimes;

    int a, b;
    unsigned int bc = FLINT_BITS - flint_clz(x);

    a = base2_32bit_pseudoprimes_start[bc - 1];
    b = base2_32bit_pseudoprimes_start[bc];

    while (a < b)
    {
        int m = a + (b - a) / 2;
        uint32_t mid_val = arr[m];

        if (mid_val == x)
            return 1;
        else if (mid_val < x)
            a = m + 1;
        else
            b = m;
    }

    return 0;
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
    uint32_t hash, base;
    uint64_t lookup;

    /* The specific hash function used to generate the table. */
    hash = ((uint32_t) (n * 314159265)) >> 15;

    FLINT_ASSERT(hash < WITNESS_BASE_HASH_SIZE);
    FLINT_ASSERT(hash / 3 < WITNESS_BASE_COMPRESSED_LEN);

    /* Unpack bits from the hash table word. */
    lookup = witness_base_tab[hash / 3];
    base = lookup >> ((hash % 3) * 21);
    base &= ((UWORD(1) << (21 + ((hash % 3) == 2))) - 1);

    if (base != 0)
        return base;

    /* If needed, binary search lookup in the oversize table. */
    int a = 0, b = NUM_OVERSIZE_BASES - 1, mid;

    while (a < b)
    {
        mid = (a + b) >> 1;
        if (oversize_index[mid] < hash)
            a = mid + 1;
        else
            b = mid;
    }

    return oversize_base[a];
}

#if FLINT_BITS == 64
#define LG_FLINT_BITS 6
#else
#define LG_FLINT_BITS 5
#endif

static ulong
_nmod_2_pow_fullword(ulong exp, nmod_t mod)
{
    ulong x, bit;
    unsigned int ebits;

    if (exp < FLINT_BITS)
        return nmod_set_ui(UWORD(1) << exp, mod);

    ebits = FLINT_BITS - flint_clz(exp);
    bit = UWORD(1) << (ebits - LG_FLINT_BITS);

    x = nmod_set_ui(UWORD(1) << (exp >> (ebits - LG_FLINT_BITS)), mod);

    while (bit >>= 1)
    {
        x = _nmod_mul_fullword(x, x, mod);
        if (bit & exp)
            x = nmod_add(x, x, mod);
    }

    return x;
}

static ulong
_nmod_2_pow_nonfullword(ulong exp, nmod_t mod)
{
    ulong x, bit;
    unsigned int ebits;

    if (exp < FLINT_BITS)
        return nmod_set_ui(UWORD(1) << exp, mod);

    ebits = FLINT_BITS - flint_clz(exp);
    bit = UWORD(1) << (ebits - LG_FLINT_BITS);

    x = nmod_set_ui(UWORD(1) << (exp >> (ebits - LG_FLINT_BITS)), mod);

    while (bit >>= 1)
    {
        x = nmod_mul(x, x, mod);
        if (bit & exp)
            x = nmod_add(x, x, mod);
    }

    return x;
}

static int
_n_is_strong_probabprime_nmod(ulong a, ulong d, nmod_t mod)
{
    ulong t = d;
    ulong n = mod.n;
    ulong y;

    FLINT_ASSERT(a < n);
    FLINT_ASSERT(a >= 2);
    FLINT_ASSERT(a != n - 1);

    if (mod.norm == 0)
    {
        if (a == 2)
            y = _nmod_2_pow_fullword(t, mod);
        else
            y = n_powmod2_ui_preinv(a, t, mod.n, mod.ninv);

        if (y == 1)
            return 1;
        t <<= 1;

        while ((t != n - 1) && (y != n - 1))
        {
            y = _nmod_mul_fullword(y, y, mod);
            t <<= 1;
        }

        return (y == n - 1);
    }
    else if (a == 2)
    {
        y = _nmod_2_pow_nonfullword(t, mod);

        if (y == 1)
            return 1;
        t <<= 1;

        while ((t != n - 1) && (y != n - 1))
        {
            y = nmod_mul(y, y, mod);
            t <<= 1;
        }

        return (y == n - 1);
    }
    else
    {
        return n_is_strong_probabprime2_preinv(mod.n, mod.ninv, a, d);
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
        ulong d, norm;
        nmod_t mod;

        d = n - 1;
        norm = flint_ctz(d);
        d >>= norm;
        nmod_init(&mod, n);

        if (!_n_is_strong_probabprime_nmod(2, d, mod))
            return 0;

        return _n_is_strong_probabprime_nmod(get_witness_base(n), d, mod);
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

