/*
    Copyright (C) 2008 Peter Shrimpton
    Copyright (C) 2009 William Hart
    Copyright (C) 2014, 2015 Dana Jacobsen
    Copyright (C) 2015 Kushagra Singh
    Copyright (C) 2025 Fredrik Johansson
    Copyright (C) 2026 Viorel Wegner

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdint.h>
#include <math.h>
#include "ulong_extras.h"
#include "nmod.h"

#define SMALL_ODDPRIME_LIMIT 32768
#define WITNESS_BASE_HASH_SIZE 32768
/* To keep this file readable, the lookup tables have been
   placed in a separate header file. */
#include "is_prime_tables.h"

// Detects if n is of the form pq where q = k*(p-1)+1
// Here we use a slightly simpler equivalence of p=(a+1) q=(k*a+1)
// k is in {2,3,4,5,6,7,8,9,12} 10 and 11 are skipped for efficiency
static int n_is_semiprime_k(ulong n)
{
// Precomputed multiplicative inverses of the sqrt of k \in [2,9]
   const double SQRTINV[8] = {
     // sqrt(2)^-1          sqrt(3)^-1
     0x1.6a09e667f3bccp-1 , 0x1.279a74590331dp-1 ,
     0x1p-1               , 0x1.c9f25c5bfedd9p-2 ,
     0x1.a20bd700c2c3fp-2 , 0x1.83091e6a7f7e6p-2 ,
     0x1.6a09e667f3bccp-2 , 0x1.5555555555555p-2 ,
   };

    // 1/sqrt(12)
  const double SQRTINV12 = 0x1.279a74590331dp-2;
  // We only need a single sqrt
  double sqrtn = sqrt(n);  
  // This is equivalent to sqrt(n/12)
  // The result of sqrt(n/k) is always less than 2^32 so it fits into a signed int64
  // Additionally float64 to signed int is generally faster than  float64 to unsigned int
  int64_t ki = sqrtn*SQRTINV12;
  // A simple interpretation change that should be zerocost
  uint64_t k = ki;  
  // check that N is of the form 
  if (((k+1)*(12*k+1))==n){
      return 1;
  }
  // Loop over the rest of the elements
   for (int idx=0;idx<8;idx++){
    // This is equivalent to sqrt(n/k)
      int64_t ai = sqrtn*SQRTINV[idx];

      uint64_t a = ai; 
      // increment the index to the offset
      uint64_t k = idx+2;
      // This is a composite so return true
      if ((a+1)*(k*a+1)==n){
         return 1;
     }
   }
  return 0;  
}

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
    
   uint32_t hash, multiplier=541578969,nsmall=n;
    
    hash = (nsmall*multiplier)>>17;
    
    return witness_base_tab[hash];
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

static int
n_is_mod30prime_small(uint32_t n)
{
    uint32_t div30 = ((uint32_t) n) / 30;
    uint32_t mod30 = ((uint32_t) n) % 30;

    static const uint8_t wheel_index[30] = { 0, 1, 0, 0, 0, 0, 0, 2, 0,
        0, 0, 3, 0, 4, 0, 0, 0, 5, 0, 6, 0, 0, 0, 7, 0, 0, 0, 0, 0, 8 };

    mod30 = wheel_index[mod30];
    if (mod30 == 0)
        return 0;

    return (flint_mod30_prime_lookup[div30] >> (mod30 - 1)) & 1;
}

int
n_is_prime_odd_no_trial(ulong n)
{
    FLINT_ASSERT(n % 2 != 0);

    if (n <= UINT32_MAX)
    {
        if (n < SMALL_ODDPRIME_LIMIT)
            return n_is_oddprime_small2(n);

        if (n < (1 << 20))
            return n_is_mod30prime_small(n);

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
        
        if (n_is_semiprime_k(n)){
           return 0;
        }
        return _n_is_strong_probabprime_redc(get_witness_base(n), d, one_red, ctx);
    }
}

int
n_is_prime(ulong n)
{
    if (n < SMALL_ODDPRIME_LIMIT)
        return (n % 2 ? n_is_oddprime_small2(n) : n == 2);

    if (!(n % 2) || !(n % 3) || !(n % 5))
        return 0;

    if (n < (1 << 20))
        return n_is_mod30prime_small(n);

    if (n > (1 << 20) &&
        (!(n % 7) || !(n % 11) || !(n % 13) || !(n % 17) || !(n % 19) ||
        !(n % 23) || !(n % 29) || !(n % 31) || !(n % 37) ||
        !(n % 41) || !(n % 43) || !(n % 47) || !(n % 53)))
        return 0;

    if (n > UINT32_MAX &&
        (!(n % 59) || !(n % 61) || !(n % 67) || !(n % 71) || !(n % 73) ||
         !(n % 79) || !(n % 83) || !(n % 89) || !(n % 97) || !(n %101) ||
         !(n %103) || !(n % 107) || !(n %109) || !(n %113) || !(n % 127) ||
         !(n %131) || !(n % 137) || !(n %139) || !(n %149)))
        return 0;

    return n_is_prime_odd_no_trial(n);
}

/* first 64 primes used for modular arithmetic */
#define N_MODULUS (UWORD(1) << (FLINT_BITS - 1))
#define N_MOD_TAB 64
static const unsigned short n_modular_primes_tab[N_MOD_TAB] = {
#if FLINT_BITS == 64
  29, 99, 123, 131, 155, 255, 269, 359, 435, 449, 453, 485, 491, 543, 585,
  599, 753, 849, 879, 885, 903, 995, 1209, 1251, 1311, 1373, 1403, 1485, 1533,
  1535, 1545, 1551, 1575, 1601, 1625, 1655, 1701, 1709, 1845, 1859, 1913,
  1995, 2045, 2219, 2229, 2321, 2363, 2385, 2483, 2499, 2523, 2543, 2613,
  2639, 2679, 2829, 2931, 3089, 3165, 3189, 3245, 3273, 3291, 3341
#else
  11, 45, 65, 95, 129, 135, 165, 209, 219, 221, 239, 245, 281, 303, 345, 351,
  359, 389, 393, 395, 413, 435, 461, 513, 519, 549, 555, 573, 575, 585, 591,
  611, 623, 629, 683, 689, 701, 729, 785, 791, 813, 843, 851, 869, 879, 893,
  905, 921, 953, 963, 965, 969, 993, 1031, 1049, 1073, 1085, 1103, 1143, 1173,
  1203, 1221, 1229, 1271
#endif
};

#if FLINT_BITS == 64
#define u64_ctz flint_ctz
#else
static unsigned int u64_ctz(uint64_t x)
{
    uint32_t lo = x;
    uint32_t hi = x >> 32;
    return (lo != 0) ? flint_ctz(lo) : flint_ctz(hi) + 32;
}
#endif

ulong n_nextprime(ulong n, int FLINT_UNUSED(proved))
{
    if (n < 32749)
    {
        if (n <= 1)
            return 2;

        int i = ((n + 1) / 2) / 64;
        int j = ((n + 1) / 2) % 64;
        uint64_t b = flint_odd_prime_lookup[i];

        if (b >> j)
            return i * 128 + 2 * u64_ctz(b >> j) + 2 * j + 1;
        else
            return (i + 1) * 128 + 2 * u64_ctz(flint_odd_prime_lookup[i + 1]) + 1;
    }
    else if (n < 1048573)
    {
        static const uint8_t wheel[8] = { 1, 7, 11, 13, 17, 19, 23, 29 };
        static const uint8_t wheel_index[30] = { 0, 0, 1, 1, 1, 1, 1, 1,
            2, 2, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7 };

        uint32_t i = ((uint32_t) (n + 1)) / 30;
        uint32_t j = wheel_index[((uint32_t) (n + 1)) % 30];
        uint32_t b = flint_mod30_prime_lookup[i];

        while ((b >> j) == 0)
        {
            i++;
            b = flint_mod30_prime_lookup[i];
            j = 0;
        }

        return 30 * i + wheel[flint_ctz(b >> j) + j];
    }
    else
    {
        ulong i, index;

        if (n >= N_MODULUS && n < N_MODULUS + n_modular_primes_tab[N_MOD_TAB-1])
        {
            for (i = 0; i < N_MOD_TAB; i++)
                if (N_MODULUS + n_modular_primes_tab[i] > n)
                    return N_MODULUS + n_modular_primes_tab[i];
        }

        if (n >= UWORD_MAX_PRIME)
        {
            flint_throw(FLINT_ERROR, "Exception (n_nextprime). No larger single-limb prime exists.\n");
        }

        static const unsigned int nextmod30[30] = {
           1, 6, 5, 4, 3, 2, 1, 4, 3, 2, 1, 2, 1, 4, 3, 2, 1, 2, 1,
           4, 3, 2, 1, 6, 5, 4, 3, 2, 1, 2
        };
        static const unsigned int nextindex[30] = {
           1, 7, 7, 7, 7, 7, 7, 11, 11, 11, 11, 13, 13, 17, 17, 17, 17, 19, 19,
           23, 23, 23, 23, 29, 29, 29, 29, 29, 29, 1
        };

        index = n % 30;
        do {
            n += nextmod30[index];
            index = nextindex[index];
        } while (!n_is_prime(n));

        return n;
    }
}

