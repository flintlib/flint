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

/* Use the division routines currently defined internally in radix.h;
   todo - move to ulong_extras */
#include "radix.h"

static ulong
n_2_powmod_c0(uint32_t exp, ulong n, const n_div_precomp_t pre)
{
    ulong x;
    int k, ebits;

    if (exp < FLINT_BITS)
        return n_rem_precomp_c0(UWORD(1) << exp, n, pre);

    ebits = FLINT_BITS - flint_clz(exp);
    x = n_rem_precomp_c0(UWORD(1) << (exp >> (ebits - 6)), n, pre);

    if (n < UWORD(1) << (FLINT_BITS / 2 - 1))
    {
        for (k = ebits - 7; k >= 0; k--)
        {
            x = n_rem_precomp_c0((ulong) x * x, n, pre);
            x <<= ((exp >> k) & 1);
        }

        if (x >= n)
            x -= n;
    }
    else
    {
        for (k = ebits - 7; k >= 0; k--)
        {
            x = n_rem_precomp_c0((ulong) x * x, n, pre);
            x <<= ((exp >> k) & 1);
            if (x >= n)
                x -= n;
        }
    }

    return x;
}

static ulong
n_2_powmod_c1(uint32_t exp, ulong n, const n_div_precomp_t pre)
{
    ulong x;
    int k, ebits;

    if (exp < FLINT_BITS)
        return n_rem_precomp_c1_bounded(UWORD(1) << exp, n, pre);

    ebits = FLINT_BITS - flint_clz(exp);
    x = n_rem_precomp_c1_bounded(UWORD(1) << (exp >> (ebits - 6)), n, pre);

    if (n < UWORD(1) << (FLINT_BITS / 2 - 1))
    {
        for (k = ebits - 7; k >= 0; k--)
        {
            x = n_rem_precomp_c1_bounded(x * x, n, pre);
            x <<= ((exp >> k) & 1);
        }

        if (x >= n)
            x -= n;
    }
    else
    {
        for (k = ebits - 7; k >= 0; k--)
        {
            x = n_rem_precomp_c1_bounded(x * x, n, pre);
            x <<= ((exp >> k) & 1);
            if (x >= n)
                x -= n;
        }
    }

    return x;
}

static int u32_is_base2_probabprime(uint32_t n)
{
    ulong d = n-1, s = 0;

    n_div_precomp_t pre;
    n_div_precomp_init(pre, n);

    s = flint_ctz(d);
    d >>= s;

    ulong t = d;
    ulong y;

    if (pre->c == 0)
    {
        y = n_2_powmod_c0(t, n, pre);
        if (y == 1)
            return 1;

        t <<= 1;
        while ((t != n - 1) && (y != n - 1))
        {
            y = n_rem_precomp_c0(y * y, n, pre);
            t <<= 1;
        }
    }
    else
    {
        y = n_2_powmod_c1(t, n, pre);
        if (y == 1)
            return 1;

        t <<= 1;
        while ((t != n - 1) && (y != n - 1))
        {
            y = n_rem_precomp_c1_bounded(y * y, n, pre);
            t <<= 1;
        }
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

/* 3, 5 (checked separately), and 32 more primes which may be vectorized */
#define NUM_ODD_TRIAL_PRIMES 34

static const uint32_t trial_inv_32[2 * NUM_ODD_TRIAL_PRIMES] = {
    UWORD(0xaaaaaaab), UWORD(0x55555555), /* 3 */
    UWORD(0xcccccccd), UWORD(0x33333333), /* 5 */
    UWORD(0xb6db6db7), UWORD(0x24924924), /* 7 */
    UWORD(0xba2e8ba3), UWORD(0x1745d174), /* 11 */
    UWORD(0xc4ec4ec5), UWORD(0x13b13b13), /* 13 */
    UWORD(0xf0f0f0f1), UWORD(0xf0f0f0f), /* 17 */
    UWORD(0x286bca1b), UWORD(0xd79435e), /* 19 */
    UWORD(0xe9bd37a7), UWORD(0xb21642c), /* 23 */
    UWORD(0x4f72c235), UWORD(0x8d3dcb0), /* 29 */
    UWORD(0xbdef7bdf), UWORD(0x8421084), /* 31 */
    UWORD(0x914c1bad), UWORD(0x6eb3e45), /* 37 */
    UWORD(0xc18f9c19), UWORD(0x63e7063), /* 41 */
    UWORD(0x2fa0be83), UWORD(0x5f417d0), /* 43 */
    UWORD(0x677d46cf), UWORD(0x572620a), /* 47 */
    UWORD(0x8c13521d), UWORD(0x4d4873e), /* 53 */
    UWORD(0xa08ad8f3), UWORD(0x456c797), /* 59 */
    UWORD(0xc10c9715), UWORD(0x4325c53), /* 61 */
    UWORD(0x7a44c6b), UWORD(0x3d22635), /* 67 */
    UWORD(0xe327a977), UWORD(0x39b0ad1), /* 71 */
    UWORD(0xc7e3f1f9), UWORD(0x381c0e0), /* 73 */
    UWORD(0x613716af), UWORD(0x33d91d2), /* 79 */
    UWORD(0x2b2e43db), UWORD(0x3159721), /* 83 */
    UWORD(0xfa3f47e9), UWORD(0x2e05c0b), /* 89 */
    UWORD(0x5f02a3a1), UWORD(0x2a3a0fd), /* 97 */
    UWORD(0x7c32b16d), UWORD(0x288df0c), /* 101 */
    UWORD(0xd3431b57), UWORD(0x27c4597), /* 103 */
    UWORD(0x8d28ac43), UWORD(0x2647c69), /* 107 */
    UWORD(0xda6c0965), UWORD(0x2593f69), /* 109 */
    UWORD(0xfdbc091), UWORD(0x243f6f0), /* 113 */
    UWORD(0xefdfbf7f), UWORD(0x2040810), /* 127 */
    UWORD(0xc9484e2b), UWORD(0x1f44659), /* 131 */
    UWORD(0x77975b9), UWORD(0x1de5d6e), /* 137 */
    UWORD(0x70586723), UWORD(0x1d77b65), /* 139 */
    UWORD(0x8ce2cabd), UWORD(0x1b7d6c3), /* 149 */
};

static const uint64_t trial_inv_64[2 * NUM_ODD_TRIAL_PRIMES] = {
    UWORD(0xaaaaaaaaaaaaaaab), UWORD(0x5555555555555555), /* 3 */
    UWORD(0xcccccccccccccccd), UWORD(0x3333333333333333), /* 5 */
    UWORD(0x6db6db6db6db6db7), UWORD(0x2492492492492492), /* 7 */
    UWORD(0x2e8ba2e8ba2e8ba3), UWORD(0x1745d1745d1745d1), /* 11 */
    UWORD(0x4ec4ec4ec4ec4ec5), UWORD(0x13b13b13b13b13b1), /* 13 */
    UWORD(0xf0f0f0f0f0f0f0f1), UWORD(0xf0f0f0f0f0f0f0f), /* 17 */
    UWORD(0x86bca1af286bca1b), UWORD(0xd79435e50d79435), /* 19 */
    UWORD(0xd37a6f4de9bd37a7), UWORD(0xb21642c8590b216), /* 23 */
    UWORD(0x34f72c234f72c235), UWORD(0x8d3dcb08d3dcb08), /* 29 */
    UWORD(0xef7bdef7bdef7bdf), UWORD(0x842108421084210), /* 31 */
    UWORD(0x14c1bacf914c1bad), UWORD(0x6eb3e45306eb3e4), /* 37 */
    UWORD(0x8f9c18f9c18f9c19), UWORD(0x63e7063e7063e70), /* 41 */
    UWORD(0x82fa0be82fa0be83), UWORD(0x5f417d05f417d05), /* 43 */
    UWORD(0x51b3bea3677d46cf), UWORD(0x572620ae4c415c9), /* 47 */
    UWORD(0x21cfb2b78c13521d), UWORD(0x4d4873ecade304d), /* 53 */
    UWORD(0xcbeea4e1a08ad8f3), UWORD(0x456c797dd49c341), /* 59 */
    UWORD(0x4fbcda3ac10c9715), UWORD(0x4325c53ef368eb0), /* 61 */
    UWORD(0xf0b7672a07a44c6b), UWORD(0x3d226357e16ece5), /* 67 */
    UWORD(0x193d4bb7e327a977), UWORD(0x39b0ad12073615a), /* 71 */
    UWORD(0x7e3f1f8fc7e3f1f9), UWORD(0x381c0e070381c0e), /* 73 */
    UWORD(0x9b8b577e613716af), UWORD(0x33d91d2a2067b23), /* 79 */
    UWORD(0xa3784a062b2e43db), UWORD(0x3159721ed7e7534), /* 83 */
    UWORD(0xf47e8fd1fa3f47e9), UWORD(0x2e05c0b81702e05), /* 89 */
    UWORD(0xa3a0fd5c5f02a3a1), UWORD(0x2a3a0fd5c5f02a3), /* 97 */
    UWORD(0x3a4c0a237c32b16d), UWORD(0x288df0cac5b3f5d), /* 101 */
    UWORD(0xdab7ec1dd3431b57), UWORD(0x27c45979c95204f), /* 103 */
    UWORD(0x77a04c8f8d28ac43), UWORD(0x2647c69456217ec), /* 107 */
    UWORD(0xa6c0964fda6c0965), UWORD(0x2593f69b02593f6), /* 109 */
    UWORD(0x90fdbc090fdbc091), UWORD(0x243f6f0243f6f02), /* 113 */
    UWORD(0x7efdfbf7efdfbf7f), UWORD(0x204081020408102), /* 127 */
    UWORD(0x3e88cb3c9484e2b), UWORD(0x1f44659e4a42715), /* 131 */
    UWORD(0xe21a291c077975b9), UWORD(0x1de5d6e3f8868a4), /* 137 */
    UWORD(0x3aef6ca970586723), UWORD(0x1d77b654b82c339), /* 139 */
    UWORD(0xdf5b0f768ce2cabd), UWORD(0x1b7d6c3dda338b2), /* 149 */
};

#define TRIAL32(n, ii) (((uint32_t) (n) * trial_inv_32[2 * ii]) <= trial_inv_32[2 * ii + 1])
#define TRIAL64(n, ii) (((uint64_t) (n) * trial_inv_64[2 * ii]) <= trial_inv_64[2 * ii + 1])

int
n_is_prime(ulong n)
{
    if (n < SMALL_ODDPRIME_LIMIT)
        return (n % 2 ? n_is_oddprime_small2(n) : n == 2);

    /* Trial division by 2, 3 and 5 */
#if FLINT_BITS == 64
    if (!(n % 2) || TRIAL64(n, 0) || TRIAL64(n, 1))
#else
    if (!(n % 2) || TRIAL32(n, 0) || TRIAL32(n, 1))
#endif
        return 0;

    if (n < (1 << 20))
        return n_is_mod30prime_small(n);

    if (n <= UINT32_MAX)
    {
        int i, c;
        /* The compiler should unroll and convert this to a few vector
           instructions when e.g. AVX2 is available; branchy return isn't
           worth it. */
        c = 0;
        for (i = 2; i < NUM_ODD_TRIAL_PRIMES; i++)
            c |= TRIAL32(n, i);

        if (c != 0)
            return 0;
    }
    else
    {
        int i;
        for (i = 2; i < NUM_ODD_TRIAL_PRIMES; i++)
            if (TRIAL64(n, i))
                return 0;
    }

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

