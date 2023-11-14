/*
    Copyright (C) 2008, 2009, David Harvey
    Copyright (C) 2021 Fredrik Johansson

This file has been adapted from the BSD-licensed bernmm package by
David Harvey (see original copyright text below).

Changes in the Arb version:
* Use FLINT functions instead of NTL functions. (Some helper functions
  have been added -- these should probably be moved to FLINT.)
* C instead of C++; some renaming and reformatting for better consistency
  with FLINT/Arb coding conventions.

Important note on performance:
* The modular arithmetic in FLINT is quite a bit slower than NTL
  since it is designed to support full 64-bit moduli which we don't need
  here. This mainly affects bernoulli_sum_powg and bernsum_pow2, which
  fortunately don't get called often for multimodular computation.
  However, we'd want to re-optimize these in case they find some other use.

*/

/*
===============================================================================

bernmm: an implementation of the algorithm described in "A multimodular
        algorithm for computing Bernoulli numbers", by David Harvey, 2008.

version 1.1

Copyright (C) 2008, 2009, David Harvey

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ``AS IS''
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

===============================================================================
*/

#include "ulong_extras.h"
#include "bernoulli.h"

#ifdef __GNUC__
# define memset __builtin_memset
#else
# include <string.h>
#endif

#define DEBUG 0
#define TIMING 1

/******************************************************************************

   Computing the main sum (general case)

******************************************************************************/

/* todo: this should be a function in ulong_extras */
ulong _bernoulli_n_muldivrem_precomp(ulong * q, ulong a, ulong b, ulong n, double bnpre)
{
    ulong qq, r;

    qq = (double) a * bnpre;
    r = a * b - qq * n;

    if ((slong) r < 0)
    {
        qq--;
        r += n;
    }

    if (r >= n)
    {
        qq++;
        r -= n;
    }

    *q = qq;
    return r;
}

/*
   Returns (1 - g^k) B_k / 2k mod p.

   PRECONDITIONS:
      5 <= p < 2^FLINT_D_BITS, p prime
      2 <= k <= p-3, k even
      pinv = n_preinvert_limb(p)
      g = a multiplicative generator of GF(p), in [0, p)
*/
static ulong
bernoulli_sum_powg(ulong p, ulong pinv, ulong k, ulong g)
{
    ulong half_gm1, sum, g_to_km1, g_to_jm1, g_to_km1_to_j, q, h;
    slong j;
    double g_pinv;

    g_pinv = (double) g / (double) p;
    half_gm1 = (g + ((g & 1) ? 0 : p) - 1) / 2;    /* (g-1)/2 mod p */
    g_to_km1 = n_powmod2_preinv(g, k-1, p, pinv);
    g_to_jm1 = 1;
    g_to_km1_to_j = g_to_km1;
    sum = 0;

    for (j = 1; j <= (p - 1) / 2; j++)
    {
        g_to_jm1 = _bernoulli_n_muldivrem_precomp(&q, g_to_jm1, g, p, g_pinv);
        h = n_submod(q, half_gm1, p);
        sum = n_submod(sum, n_mulmod2_preinv(h, g_to_km1_to_j, p, pinv), p);
        g_to_km1_to_j = n_mulmod2_preinv(g_to_km1_to_j, g_to_km1, p, pinv);
    }

    return sum;
}


/******************************************************************************

   Computing the main sum (c = 1/2 case)

******************************************************************************/

/*
    The Expander class stores precomputed information for a fixed integer p,
    that subsequently permits fast computation of the binary expansion of s/p
    for any 0 < s < p.

    The constructor takes p and max_words as input. Must have 1 <= max_words <=
    MAX_INV. It computes an approximation to 1/p.

    The function expand(mp_ptr res, long s, long n) computes n limbs of s/p.
    Must have 0 < s < p and 1 <= n <= max_words. The output is written to res.
    The first word of output is junk. The next n words are the digits of s/p,
    from least to most significant. The buffer must be at least n+2 words long
    (even though the first and last words are never used for output).
*/

#define MAX_INV 256

typedef struct
{
   /* Approximation to 1/p. We store (max_words + 1) limbs. */
   mp_limb_t pinv[MAX_INV + 2];
   mp_limb_t p;
   int max_words;
}
expander_t;

static void
expander_init(expander_t * this, ulong p, int max_words)
{
    mp_limb_t one;

    FLINT_ASSERT(max_words >= 1);
    FLINT_ASSERT(max_words <= MAX_INV);

    this->max_words = max_words;
    this->p = p;
    one = 1;
    mpn_divrem_1(this->pinv, max_words + 1, &one, 1, p);
}

static void
expander_expand(mp_ptr res, expander_t * this, ulong s, ulong n)
{
    slong i;

    FLINT_ASSERT(s > 0 && s < this->p);
    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(n <= this->max_words);

    if (s == 1)
    {
        /* already have 1/p; just copy it */
        for (i = 1; i <= n; i++)
            res[i] = this->pinv[this->max_words - n + i];
    }
    else
    {
        mpn_mul_1(res, this->pinv + this->max_words - n, n + 1, (mp_limb_t) s);

        /* If the first output limb is really close to 0xFFFF..., then there's
           a possibility of overflow, so fall back on doing division directly.
           This should happen extremely rarely --- essentially never on a
           64-bit system, and very occasionally on a 32-bit system. */
        if (res[0] > -((mp_limb_t) s))
        {
            mp_limb_t ss = s;
            mpn_divrem_1(res, n + 1, &ss, 1, this->p);
        }
    }
}


/*
   Returns (2^(-k) - 1) 2 B_k / k  mod p.

   (Note: this is useless if 2^k = 1 mod p.)

   PRECONDITIONS:
      5 <= p, p prime
      2 <= k <= p-3, k even
      pinv = n_preinvert_limb(p)
      g = a multiplicative generator of GF(p), in [0, p)
      n = multiplicative n_multiplicative_order of 2 in GF(p)
*/

#define TABLE_LG_SIZE 8
#define TABLE_SIZE (WORD(1) << TABLE_LG_SIZE)
#define TABLE_MASK (TABLE_SIZE - 1)
#define NUM_TABLES (FLINT_BITS / TABLE_LG_SIZE)

#if FLINT_BITS % TABLE_LG_SIZE != 0
#error Number of bits in a ulong must be divisible by TABLE_LG_SIZE
#endif

ulong bernsum_pow2(ulong p, ulong pinv, ulong k, ulong g, ulong n)
{
    slong i, m;
    ulong g_to_km1, two_to_km1, B_to_km1, s_jump;
    ulong tables[NUM_TABLES][TABLE_SIZE];
    ulong g_to_km1_to_i;
    ulong g_to_i;
    ulong sum;
    expander_t expander;

    slong h;
    ulong x;
    ulong weights[TABLE_SIZE];
    ulong x_jump;

    /* In the main summation loop we accumulate data into the _tables_ array;
      tables[y][z] contributes to the final answer with a weight of

      sum(-(-1)^z[t] * (2^(k-1))^(FLINT_BITS - 1 - y * TABLE_LG_SIZE - t) :
                                                       0 <= t < TABLE_LG_SIZE),

      where z[t] denotes the t-th binary digit of z (LSB is t = 0).
      The memory footprint for _tables_ is 4KB on a 32-bit machine, or 16KB
      on a 64-bit machine, so should fit easily into L1 cache. */
    memset(tables, 0, sizeof(ulong) * NUM_TABLES * TABLE_SIZE);

    m = (p-1) / n;

    /* take advantage of symmetry (n' and m' from the paper) */
    if (n & 1)
        m >>= 1;
    else
        n >>= 1;

    /* g^(k-1) */
    g_to_km1 = n_powmod2_preinv(g, k - 1, p, pinv);
    /* 2^(k-1) */
    two_to_km1 = n_powmod2_preinv(2, k - 1, p, pinv);
    /* B^(k-1), where B = 2^FLINT_BITS */
    B_to_km1 = n_powmod2_preinv(two_to_km1, FLINT_BITS, p, pinv);
    /* B^(MAX_INV) */
    s_jump = n_powmod2_preinv(2, MAX_INV * FLINT_BITS, p, pinv);

    /* todo - help speed up modmuls
    mulmod_precon_t g_pinv = PrepMulModPrecon(g, p, pinv);
    mulmod_precon_t g_to_km1_pinv = PrepMulModPrecon(g_to_km1, p, pinv);
    mulmod_precon_t two_to_km1_pinv = PrepMulModPrecon(two_to_km1, p, pinv);
    mulmod_precon_t B_to_km1_pinv = PrepMulModPrecon(B_to_km1, p, pinv);
    mulmod_precon_t s_jump_pinv = PrepMulModPrecon(s_jump, p, pinv);
    */

    g_to_km1_to_i = 1;
    g_to_i = 1;
    sum = 0;

    /* Precompute some of the binary expansion of 1/p; at most MAX_INV words,
      or possibly less if n is sufficiently small */
    expander_init(&expander, p, (n >= MAX_INV * FLINT_BITS)
                                       ? MAX_INV : ((n - 1) / FLINT_BITS + 1));

    /* =========== phase 1: main summation loop */

    /* loop over outer sum */
    for (i = 0; i < m; i++)
    {
        ulong s, x, y;
        slong nn;

        /* s keeps track of g^i*2^j mod p */
        s = g_to_i;
        /* x keeps track of (g^i*2^j)^(k-1) mod p */
        x = g_to_km1_to_i;

        /* loop over inner sum; break it up into chunks of length at most */
        /* MAX_INV * FLINT_BITS. If n is large, this allows us to do most of */
        /* the work with mpn_mul_1 instead of mpn_divrem_1, and also improves */
        /* memory locality. */
        for (nn = n; nn > 0; nn -= MAX_INV * FLINT_BITS)
        {
            mp_limb_t s_over_p[MAX_INV + 2];
            slong bits, words;
            mp_ptr next;

            if (nn >= MAX_INV * FLINT_BITS)
            {
                /* do one chunk of length exactly MAX_INV * FLINT_BITS */
                bits = MAX_INV * FLINT_BITS;
                words = MAX_INV;
            }
            else
            {
                /* last chunk of length less than MAX_INV * FLINT_BITS */
                bits = nn;
                words = (nn - 1) / FLINT_BITS + 1;
            }

            /* compute some bits of the binary expansion of s/p */
            expander_expand(s_over_p, &expander, s, words);

            next = s_over_p + words;

            /* loop over whole words */
            for (; bits >= FLINT_BITS; bits -= FLINT_BITS, next--)
            {
                mp_limb_t y;
#if NUM_TABLES != 8 && NUM_TABLES != 4
                mp_ptr target;
#else
                mp_ptr target0, target1, target2, target3, target4, target5, target6, target7;
#endif

                y = *next;

#if NUM_TABLES != 8 && NUM_TABLES != 4
                /* generic version */
                for (h = 0; h < NUM_TABLES; h++)
                {
                    target = &(tables[h][y & TABLE_MASK]);
                    *target = n_submod(*target, x, p);
                    y >>= TABLE_LG_SIZE;
                }
#else
                /* unrolled versions for 32-bit/64-bit machines */
                target0 = &(tables[0][y & TABLE_MASK]);
                *target0 = n_submod(*target0, x, p);

                target1 = &(tables[1][(y >> TABLE_LG_SIZE) & TABLE_MASK]);
                *target1 = n_submod(*target1, x, p);

                target2 = &(tables[2][(y >> (2*TABLE_LG_SIZE)) & TABLE_MASK]);
                *target2 = n_submod(*target2, x, p);

                target3 = &(tables[3][(y >> (3*TABLE_LG_SIZE)) & TABLE_MASK]);
                *target3 = n_submod(*target3, x, p);
#if NUM_TABLES == 8
                target4 = &(tables[4][(y >> (4*TABLE_LG_SIZE)) & TABLE_MASK]);
                *target4 = n_submod(*target4, x, p);

                target5 = &(tables[5][(y >> (5*TABLE_LG_SIZE)) & TABLE_MASK]);
                *target5 = n_submod(*target5, x, p);

                target6 = &(tables[6][(y >> (6*TABLE_LG_SIZE)) & TABLE_MASK]);
                *target6 = n_submod(*target6, x, p);

                target7 = &(tables[7][(y >> (7*TABLE_LG_SIZE)) & TABLE_MASK]);
                *target7 = n_submod(*target7, x, p);
#endif
#endif

                x = n_mulmod2_preinv(x, B_to_km1, p, pinv);  /* todo: precond */
            }

            /* loop over remaining bits in the last word */
            y = *next;
            for (; bits > 0; bits--)
            {
                if (y & (UWORD(1) << (FLINT_BITS - 1)))
                    sum = n_submod(sum, x, p);
                else
                    sum = n_addmod(sum, x, p);

                x = n_mulmod2_preinv(x, two_to_km1, p, pinv);   /* todo: precond */
                y <<= 1;
            }

            s = n_mulmod2_preinv(s, s_jump, p, pinv);   /* todo: precond */
        }

        /* update g^i and (g^(k-1))^i */
        g_to_i = n_mulmod2_preinv(g_to_i, g, p, pinv);
        g_to_km1_to_i = n_mulmod2_preinv(g_to_km1_to_i, g_to_km1, p, pinv);
    }

#if DEBUG
    {
        slong i, j;
        for (i = 0; i < NUM_TABLES; i++)
        {
            printf("tab[%lu] = ", i);
            for (j = 0; j < TABLE_SIZE; j++)
                printf("%lu ", tables[i][j]);
            printf("\n");
        }
    }
#endif


    /* =========== phase 2: consolidate table data */

    /* compute weights[z] = sum((-1)^z[t] * (2^(k-1))^(TABLE_LG_SIZE - 1 - t) :
                                                       0 <= t < TABLE_LG_SIZE). */

    weights[0] = 0;
    for (h = 0, x = 1; h < TABLE_LG_SIZE;
         h++, x = n_mulmod2_preinv(x, two_to_km1, p, pinv))
    {
        for (i = (WORD(1) << h) - 1; i >= 0; i--)
        {
            weights[2*i+1] = n_submod(weights[i], x, p);
            weights[2*i]   = n_addmod(weights[i], x, p);
        }
    }

    /* combine table data with weights */

    x_jump = n_powmod2_preinv(two_to_km1, TABLE_LG_SIZE, p, pinv);

    for (h = NUM_TABLES - 1, x = 1; h >= 0; h--)
    {
        for (i = 0; i < TABLE_SIZE; i++)
        {
            ulong y = n_mulmod2_preinv(tables[h][i], weights[i], p, pinv);
            y = n_mulmod2_preinv(y, x, p, pinv);
            sum = n_submod(sum, y, p);
        }

        x = n_mulmod2_preinv(x_jump, x, p, pinv);
    }

    return sum;
}

/******************************************************************************

   Computing the main sum (c = 1/2 case, with REDC arithmetic)

   Throughout this section F denotes 2^(FLINT_BITS / 2).

******************************************************************************/


/*
   Returns x/F mod n. Output is in [0, 2n), i.e. *not* reduced completely
   into [0, n).

   PRECONDITIONS:
      3 <= n < F, n odd
      0 <= x < nF    (if n < F/2)
      0 <= x < nF/2  (if n > F/2)
      ninv2 = -1/n mod F
*/
#define LOW_MASK ((UWORD(1) << (FLINT_BITS / 2)) - 1)
static inline ulong RedcFast(ulong x, ulong n, ulong ninv2)
{
    ulong y = (x * ninv2) & LOW_MASK;
    ulong z = x + (n * y);
    return z >> (FLINT_BITS / 2);
}

/*
   Same as RedcFast(), but reduces output into [0, n).
*/
static inline ulong Redc(ulong x, ulong n, ulong ninv2)
{
    ulong y = RedcFast(x, n, ninv2);
    if (y >= n)
       y -= n;
    return y;
}

/*
   Computes -1/n mod F, in [0, F).

   PRECONDITIONS:
      3 <= n < F, n odd
*/
static ulong PrepRedc(ulong n)
{
    ulong bits;
    ulong ninv2 = -n;   /* already correct mod 8 */

    /* Newton's method for 2-adic inversion */
    for (bits = 3; bits < FLINT_BITS/2; bits *= 2)
       ninv2 = 2*ninv2 + n * ninv2 * ninv2;

    return ninv2 & LOW_MASK;
}

/*
   Same as bernsum_pow2(), but uses REDC arithmetic, and various delayed
   reduction strategies.

   PRECONDITIONS:
      Same as bernsum_pow2(), and in addition:
      p < 2^(FLINT_BITS/2 - 1)

   (See bernsum_pow2() for code comments; we only add comments here where
   something is different from bernsum_pow2())
*/
ulong bernsum_pow2_redc(ulong p, ulong pinv, ulong k, ulong g, ulong n)
{
    ulong pinv2 = PrepRedc(p);
    ulong F = (UWORD(1) << (FLINT_BITS/2)) % p;
    ulong x;
    slong h, i, m;
    ulong weights[TABLE_SIZE];
    ulong x_jump;
    ulong x_jump_redc;

    ulong g_to_km1;
    ulong two_to_km1;
    ulong B_to_km1;
    ulong s_jump;

    ulong g_redc;
    ulong g_to_km1_redc;
    ulong two_to_km1_redc;
    ulong B_to_km1_redc;
    ulong s_jump_redc;

    ulong g_to_km1_to_i;
    ulong g_to_i;
    ulong sum;

    ulong tables[NUM_TABLES][TABLE_SIZE];

    expander_t expander;

    memset(tables, 0, sizeof(ulong) * NUM_TABLES * TABLE_SIZE);

    m = (p-1) / n;

    if (n & 1)
        m >>= 1;
    else
        n >>= 1;

    g_to_km1 = n_powmod2_preinv(g, k-1, p, pinv);
    two_to_km1 = n_powmod2_preinv(2, k-1, p, pinv);
    B_to_km1 = n_powmod2_preinv(two_to_km1, FLINT_BITS, p, pinv);
    s_jump = n_powmod2_preinv(2, MAX_INV * FLINT_BITS, p, pinv);

    g_redc = n_mulmod2_preinv(g, F, p, pinv);
    g_to_km1_redc = n_mulmod2_preinv(g_to_km1, F, p, pinv);
    two_to_km1_redc = n_mulmod2_preinv(two_to_km1, F, p, pinv);
    B_to_km1_redc = n_mulmod2_preinv(B_to_km1, F, p, pinv);
    s_jump_redc = n_mulmod2_preinv(s_jump, F, p, pinv);

    g_to_km1_to_i = 1;    /* always in [0, 2p)  */
    g_to_i = 1;           /* always in [0, 2p)  */
    sum = 0;

#if DEBUG
    printf("%lu %lu %lu %lu %lu  %lu %lu %lu %lu\n", F, g_to_km1, two_to_km1, B_to_km1, s_jump, g_redc, g_to_km1_redc, B_to_km1_redc, s_jump_redc);
#endif

    expander_init(&expander, p, (n >= MAX_INV * FLINT_BITS)
                                       ? MAX_INV : ((n - 1) / FLINT_BITS + 1));

    /* =========== phase 1: main summation loop */

    for (i = 0; i < m; i++)
    {
        ulong s, x, y;
        slong nn, bits, words;
        mp_ptr next;

        s = g_to_i;           /* always in [0, p) */
        if (s >= p)
            s -= p;

        x = g_to_km1_to_i;    /* always in [0, 2p) */

        for (nn = n; nn > 0; nn -= MAX_INV * FLINT_BITS)
        {
            ulong s_over_p[MAX_INV + 2];

            if (nn >= MAX_INV * FLINT_BITS)
            {
                bits = MAX_INV * FLINT_BITS;
                words = MAX_INV;
            }
            else
            {
                bits = nn;
                words = (nn - 1) / FLINT_BITS + 1;
            }

            expander_expand(s_over_p, &expander, s, words);
            next = s_over_p + words;

            for (; bits >= FLINT_BITS; bits -= FLINT_BITS, next--)
            {
                y = *next;

#if DEBUG
                printf("i = %lu  nn = %lu  words = %lu  bits = %lu  y = %lu\n", i, nn, words, bits, y);
#endif

                /* note: we add the values into tables *without* reduction mod p */

#if NUM_TABLES != 8 && NUM_TABLES != 4
                /* generic version */
                for (h = 0; h < NUM_TABLES; h++)
                {
                    tables[h][y & TABLE_MASK] += x;
                    y >>= TABLE_LG_SIZE;
                }
#else
                /* unrolled versions for 32-bit/64-bit machines */
                tables[0][ y                       & TABLE_MASK] += x;
                tables[1][(y >>    TABLE_LG_SIZE ) & TABLE_MASK] += x;
                tables[2][(y >> (2*TABLE_LG_SIZE)) & TABLE_MASK] += x;
                tables[3][(y >> (3*TABLE_LG_SIZE)) & TABLE_MASK] += x;
#if NUM_TABLES == 8
                tables[4][(y >> (4*TABLE_LG_SIZE)) & TABLE_MASK] += x;
                tables[5][(y >> (5*TABLE_LG_SIZE)) & TABLE_MASK] += x;
                tables[6][(y >> (6*TABLE_LG_SIZE)) & TABLE_MASK] += x;
                tables[7][(y >> (7*TABLE_LG_SIZE)) & TABLE_MASK] += x;
#endif
#endif

                x = RedcFast(x * B_to_km1_redc, p, pinv2);
            }

            /* bring x into [0, p) for next loop */
            if (x >= p)
                x -= p;

            y = *next;
            for (; bits > 0; bits--)
            {
                if (y & (UWORD(1) << (FLINT_BITS - 1)))
                    sum = n_submod(sum, x, p);
                else
                    sum = n_addmod(sum, x, p);

                x = Redc(x * two_to_km1_redc, p, pinv2);
                y <<= 1;
            }

            s = Redc(s * s_jump_redc, p, pinv2);
        }

        g_to_i = RedcFast(g_to_i * g_redc, p, pinv2);
        g_to_km1_to_i = RedcFast(g_to_km1_to_i * g_to_km1_redc, p, pinv2);
    }

    /* At this point, each table entry is at most p^2 (since x was always */
    /* in [0, 2p), and the inner loop was called at most (p/2) / FLINT_BITS */
    /* times, and 2p * p/2 / FLINT_BITS * TABLE_LG_SIZE <= p^2). */

#if DEBUG
    {
        slong i, j;
        for (i = 0; i < NUM_TABLES; i++)
        {
            printf("tab[%lu] = ", i);
            for (j = 0; j < TABLE_SIZE; j++)
                printf("%lu ", tables[i][j]);
            printf("\n");
        }
    }
#endif

    /* =========== phase 2: consolidate table data */

    weights[0] = 0;
    /* we store the weights multiplied by a factor of 2^(3*FLINT_BITS/2) to */
    /* compensate for the three rounds of REDC reduction in the loop below */
    for (h = 0, x = n_powmod2_preinv(2, 3*FLINT_BITS/2, p, pinv);
         h < TABLE_LG_SIZE; h++, x = Redc(x * two_to_km1_redc, p, pinv2))
    {
        for (i = (WORD(1) << h) - 1; i >= 0; i--)
        {
            weights[2*i+1] = n_submod(weights[i], x, p);
            weights[2*i]   = n_addmod(weights[i], x, p);
        }
    }

    x_jump = n_powmod2_preinv(two_to_km1, TABLE_LG_SIZE, p, pinv);
    x_jump_redc = n_mulmod2_preinv(x_jump, F, p, pinv);

    for (h = NUM_TABLES - 1, x = 1; h >= 0; h--)
    {
        for (i = 0; i < TABLE_SIZE; i++)
        {
            ulong y;
            y = RedcFast(tables[h][i], p, pinv2);
            y = RedcFast(y * weights[i], p, pinv2);
            y = RedcFast(y * x, p, pinv2);
            sum += y;
        }

        x = Redc(x * x_jump_redc, p, pinv2);
    }

    return sum % p;
}


/******************************************************************************

   Wrappers for bernsum_*

******************************************************************************/

/*
    Returns B_k/k mod p, in the range [0, p).

    PRECONDITIONS:
        5 <= p < NTL_SP_BOUND, p prime
        2 <= k <= p-3, k even
        pinv = PrepMulMod(p)

    Algorithm: uses bernoulli_sum_powg() to compute the main sum.
*/
ulong _bernoulli_mod_p_harvey_powg(ulong p, ulong pinv, ulong k)
{
    ulong x, g, t, g_to_k;

    g = n_primitive_root_prime(p);

    /* compute main sum */
    x = bernoulli_sum_powg(p, pinv, k, g);

    /* divide by (1 - g^k) and multiply by 2 */
    g_to_k = n_powmod2_preinv(g, k, p, pinv);
    t = n_invmod(p + 1 - g_to_k, p);
    x = n_mulmod2_preinv(x, t, p, pinv);
    x = n_addmod(x, x, p);

    return x;
}

static ulong
n_multiplicative_order(ulong x, ulong p, ulong pinv, n_factor_t * F)
{
    ulong m, q, mm;
    slong i;

    /* in the loop below, m is always some multiple of the n_multiplicative_order of x */
    m = p - 1;

    /* try to remove factors from m until we can't remove any more */
    for (i = 0; i < F->num; i++)
    {
        q = F->p[i];

        while (m % q == 0)
        {
            mm = m / q;
            if (n_powmod2_preinv(x, mm, p, pinv) != 1)
                break;
            m = mm;
        }
    }

    return m;
}

/*
   Returns B_k/k mod p, in the range [0, p).

   PRECONDITIONS:
      5 <= p < NTL_SP_BOUND, p prime
      2 <= k <= p-3, k even
      pinv = PrepMulMod(p)
      2^k != 1 mod p

   Algorithm: uses bernsum_pow2() (or bernsum_pow2_redc() if p is small
   enough) to compute the main sum.
*/
ulong _bernoulli_mod_p_harvey_pow2(ulong p, ulong pinv, ulong k)
{
    n_factor_t F;
    ulong g, n, x, t;

    n_factor_init(&F);
    n_factor(&F, p - 1, 1);

    g = n_primitive_root_prime_prefactor(p, &F);
    n = n_multiplicative_order(2, p, pinv, &F);

#if DEBUG
    printf("g = %lu, n = %lu\n", g, n);
#endif

    if (p < (UWORD(1) << (FLINT_BITS/2 - 1)))
        x = bernsum_pow2_redc(p, pinv, k, g, n);
    else
        x = bernsum_pow2(p, pinv, k, g, n);

    /* divide by 2*(2^(-k) - 1) */
    t = n_submod(n_invmod(n_powmod2_preinv(2, k, p, pinv), p), 1, p);
    t = n_addmod(t, t, p);
    t = n_invmod(t, p);
    x = n_mulmod2_preinv(x, t, p, pinv);

    return x;
}

/*
   Returns B_k/k mod p, in the range [0, p).

   PRECONDITIONS:
      5 <= p < NTL_SP_BOUND, p prime
      2 <= k <= p-3, k even
      pinv = PrepMulMod(p)
*/
ulong _bernoulli_mod_p_harvey(ulong p, ulong pinv, ulong k)
{
    if (n_powmod2_preinv(2, k, p, pinv) != 1)
    {
        /* 2^k != 1 mod p, so we use the faster version */
        return _bernoulli_mod_p_harvey_pow2(p, pinv, k);
    }
    else
    {
          /* forced to use slower version */
        return _bernoulli_mod_p_harvey_powg(p, pinv, k);
    }
}

/******************************************************************************

   Main bernoulli_mod_p routine

******************************************************************************/

/* beware: opposite argument order to bernmm */
ulong bernoulli_mod_p_harvey(ulong k, ulong p)
{
    ulong m, x, pinv;

    FLINT_ASSERT(k >= 0);
    FLINT_ASSERT(2 <= p && p < (UWORD(1) << FLINT_D_BITS));

    if (k == 0)
        return 1;

    if (k == 1)
    {
        if (p == 2)
            return -1;
        return (p - 1) / 2;
    }

    if (k & 1)
        return 0;

    /* denominator of B_k is always divisible by 6 for k >= 2 */
    if (p <= 3)
        return UWORD_MAX;

    /* use Kummer's congruence (k = m mod p-1  =>  B_k/k = B_m/m mod p) */
    m = k % (p - 1);
    if (m == 0)
        return UWORD_MAX;

    pinv = n_preinvert_limb(p);
    x = _bernoulli_mod_p_harvey(p, pinv, m);    /* = B_m/m mod p  */

    return n_mulmod2_preinv(x, k % p, p, pinv);
}

