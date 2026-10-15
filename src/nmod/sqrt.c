/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdint.h>
#include "ulong_extras.h"
#include "nmod.h"

/*
    Square roots modulo an odd prime n, computed entirely in Montgomery form.

    The arithmetic is selected by a compile-time constant so that each
    variant gets its own fully inlined code:

      SQRT_HALF  n < 2^30: half-limb Montgomery, lazy reduction to [0, 2n)
      SQRT_FAST  n < 2^62: full-limb Montgomery, lazy reduction to [0, 2n)
      SQRT_FULL  otherwise: full-limb Montgomery, reduced to [0, n)

    No Jacobi symbol is computed for the input. Instead, a candidate root
    r is always produced and the result is accepted only if r^2 = a. This
    costs one multiplication, is a complete residuosity test when n is
    prime, and guarantees that a returned root is genuine even when n is
    not prime.
*/

#define SQRT_FULL 0
#define SQRT_FAST 1
#define SQRT_HALF 2

/*
    Modular subtraction x - y mod m for x, y < m, written with a mask
    rather than a comparison so that the compiler emits no branch. Plain
    conditional expressions can end up as branches at -O3 (e.g. through
    path splitting), and in the loops below the conditions are random,
    so mispredictions would dominate the running time.
*/
FLINT_FORCE_INLINE ulong
_sqrt_submod(ulong x, ulong y, ulong m)
{
    return (x - y) + (m & (-(ulong) (x < y)));
}

FLINT_FORCE_INLINE ulong
_sqrt_addmod(ulong x, ulong y, ulong m)
{
    return _sqrt_submod(x, m - y, m);
}

FLINT_FORCE_INLINE ulong
_sqrt_mul(ulong x, ulong y, const nmod_redc_ctx_t ctx, int v)
{
    if (v == SQRT_HALF)
    {
        return nmod_redc_half_fast_mul(x, y, ctx);
    }
    else if (v == SQRT_FAST)
    {
        return nmod_redc_fast_mul(x, y, ctx);
    }
    else
    {
        /* nmod_redc_mul, without a branch */
        ulong n = ctx->mod.n;
        ull_t t = ull_u_mul_u(x, y);
        ulong lo = n_mulhi(ull_lo(t) * (-ctx->nred), n);
        return _sqrt_submod(ull_hi(t), lo, n);
    }
}

FLINT_FORCE_INLINE ulong
_sqrt_add(ulong x, ulong y, const nmod_redc_ctx_t ctx, int v)
{
    if (v == SQRT_FULL)
        return _sqrt_addmod(x, y, ctx->mod.n);
    else
        return _sqrt_addmod(x, y, 2 * ctx->mod.n);
}

FLINT_FORCE_INLINE ulong
_sqrt_sub(ulong x, ulong y, const nmod_redc_ctx_t ctx, int v)
{
    if (v == SQRT_FULL)
        return _sqrt_submod(x, y, ctx->mod.n);
    else
        return _sqrt_submod(x, y, 2 * ctx->mod.n);
}

/* to [0, n) */
FLINT_FORCE_INLINE ulong
_sqrt_normalise(ulong x, const nmod_redc_ctx_t ctx, int v)
{
    ulong n = ctx->mod.n;

    if (v == SQRT_FULL)
        return x;
    else
        return x - (n & (-(ulong) (x >= n)));
}

/* Montgomery form of x in [0, n), as an element of [0, n) */
FLINT_FORCE_INLINE ulong
_sqrt_set_nmod(ulong x, const nmod_redc_ctx_t ctx, int v)
{
    if (v == SQRT_HALF)
        return nmod_redc_half_set_nmod(x, ctx);
    else
        return nmod_redc_set_nmod(x, ctx);
}

FLINT_FORCE_INLINE ulong
_sqrt_get_nmod(ulong x, const nmod_redc_ctx_t ctx, int v)
{
    if (v == SQRT_HALF)
        return nmod_redc_half_get_nmod(x, ctx);
    else
        return nmod_redc_get_nmod(x, ctx);
}

/* x if bit is set, otherwise y; without a branch, as exponent bits
   are unpredictable when the modulus varies between calls */
FLINT_FORCE_INLINE ulong
_sqrt_select(ulong bit, ulong x, ulong y)
{
    return y ^ ((x ^ y) & (-(bit & 1)));
}

/* a^e for e >= 1. The multiplications into x are off the critical path
   (the squaring chain), so multiplying by one instead of branching
   is essentially free. */
FLINT_FORCE_INLINE ulong
_sqrt_pow(ulong a, ulong e, ulong one, const nmod_redc_ctx_t ctx, int v)
{
    ulong x;

    FLINT_ASSERT(e != 0);

    x = _sqrt_select(e, a, one);

    while (e >>= 1)
    {
        a = _sqrt_mul(a, a, ctx, v);
        x = _sqrt_mul(x, _sqrt_select(e, a, one), ctx, v);
    }

    return x;
}

/*
    a^e and b^e for e >= 1. The two squaring chains are independent,
    so interleaving them hides most of the multiplication latency.
*/
FLINT_FORCE_INLINE void
_sqrt_pow2(ulong * ra, ulong * rb, ulong a, ulong b, ulong e, ulong one,
                            const nmod_redc_ctx_t ctx, int v)
{
    ulong x, y;

    FLINT_ASSERT(e != 0);

    x = _sqrt_select(e, a, one);
    y = _sqrt_select(e, b, one);

    while (e >>= 1)
    {
        a = _sqrt_mul(a, a, ctx, v);
        b = _sqrt_mul(b, b, ctx, v);
        x = _sqrt_mul(x, _sqrt_select(e, a, one), ctx, v);
        y = _sqrt_mul(y, _sqrt_select(e, b, one), ctx, v);
    }

    *ra = x;
    *rb = y;
}

/*
    For the odd primes k = 3, ..., 61: kinv = 1/k mod B (B = 2^FLINT_BITS),
    and a bitmask (64 bits, also on 32-bit machines) with bit h set iff
    -h is a nonzero square modulo k.

    For any n, with m = n kinv mod B, we have m k = n + h B where h is the
    high word of m k, so h = -n / B mod k with 0 <= h < k. Since B is a
    square, (n | k) = (-h | k), which is read off from the mask. This
    tests n modulo k using two multiplications instead of a division.
*/
static const unsigned char _sqrt_small_primes[] =
    { 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61 };

static const ulong _sqrt_small_inv[] =
{
#if FLINT_BITS == 64
    UWORD(0xaaaaaaaaaaaaaaab), UWORD(0xcccccccccccccccd), UWORD(0x6db6db6db6db6db7),
    UWORD(0x2e8ba2e8ba2e8ba3), UWORD(0x4ec4ec4ec4ec4ec5), UWORD(0xf0f0f0f0f0f0f0f1),
    UWORD(0x86bca1af286bca1b), UWORD(0xd37a6f4de9bd37a7), UWORD(0x34f72c234f72c235),
    UWORD(0xef7bdef7bdef7bdf), UWORD(0x14c1bacf914c1bad), UWORD(0x8f9c18f9c18f9c19),
    UWORD(0x82fa0be82fa0be83), UWORD(0x51b3bea3677d46cf), UWORD(0x21cfb2b78c13521d),
    UWORD(0xcbeea4e1a08ad8f3), UWORD(0x4fbcda3ac10c9715)
#else
    UWORD(0xaaaaaaab), UWORD(0xcccccccd), UWORD(0xb6db6db7), UWORD(0xba2e8ba3),
    UWORD(0xc4ec4ec5), UWORD(0xf0f0f0f1), UWORD(0x286bca1b), UWORD(0xe9bd37a7),
    UWORD(0x4f72c235), UWORD(0xbdef7bdf), UWORD(0x914c1bad), UWORD(0xc18f9c19),
    UWORD(0x2fa0be83), UWORD(0x677d46cf), UWORD(0x8c13521d), UWORD(0xa08ad8f3),
    UWORD(0xc10c9715)
#endif
};

static const uint64_t _sqrt_small_qr_mask[] =
{
    UINT64_C(0x4), UINT64_C(0x12), UINT64_C(0x68), UINT64_C(0x5c4),
    UINT64_C(0x161a), UINT64_C(0x1a316), UINT64_C(0x4f50c), UINT64_C(0x7acca0),
    UINT64_C(0x13d122f2), UINT64_C(0x6de2b848), UINT64_C(0x165e211e9a), UINT64_C(0x1b382b50736),
    UINT64_C(0x4a77c5c11ac), UINT64_C(0x7bcae4d8ac20), UINT64_C(0x12dd703303aed2), UINT64_C(0x5d49de7c1846d44),
    UINT64_C(0x1713e6940a59f23a)
};

/*
    Some k with Jacobi symbol (k | n) = -1, for odd n = 1 mod 8; returns 0
    if none is found (which can only happen if n is not prime). Since n =
    1 mod 4, reciprocity gives (k | n) = (n | k) for odd prime k. 2 is a
    residue for n = 1 mod 8, and the least nonresidue is prime, so only
    odd primes need testing. This runs once per call, so it is kept out of
    line to avoid disturbing code generation in the callers.
*/
FLINT_STATIC_NOINLINE ulong
_sqrt_nonresidue_1mod8(ulong n)
{
    ulong k, h;
    slong i;

    FLINT_ASSERT(n % 8 == 1);

    for (i = 0; i < (slong) (sizeof(_sqrt_small_primes) / sizeof(_sqrt_small_primes[0])); i++)
    {
        k = _sqrt_small_primes[i];
        h = n_mulhi(n * _sqrt_small_inv[i], k);
        if (!((_sqrt_small_qr_mask[i] >> h) & 1))
            return k;
    }

    for (k = 67; ; k += 2)
    {
        if (n_jacobi_unsigned(k, n) == -1)
            return k;

        /* For prime n < 2^64 the least nonresidue is far smaller than
           this; only a square n has no nonresidue at all. */
        if (k == 1025 && n_is_square(n))
            return 0;

        if (k >= (UWORD(1) << 20))
            return 0;
    }
}

/*
    Mueller's variant of Cipolla-Lehmer, for n = 1 mod 4. With t chosen so
    that (a t^2 - 4 | n) = -1, a square root of a is V_{(n-1)/4}(P) / t
    where P = a t^2 - 2 and V_k is the Lucas sequence V_0 = 2, V_1 = P,
    V_{k+1} = P V_k - V_{k-1}. The ladder computing V_k performs one
    squaring and one independent multiplication per bit, so its latency is
    that of a single exponentiation regardless of the 2-adic valuation of
    n - 1, at the price of a few Jacobi symbols to find t.

    We restrict to t = 2^j so that dividing by t is just halving. Candidate
    discriminants are a - 4 (j = 0) and 4 c_j with c_j = 4^(j-1) a - 1
    (j >= 1), where c_{j+1} = 4 c_j + 3.

    Here a is in standard form and one is in Montgomery form. Returns 0
    if no suitable t was found (the caller should fall back on something
    else), otherwise sets *res to a candidate root in Montgomery form.
*/
FLINT_FORCE_INLINE int
_sqrt_mueller(ulong * res, ulong a, ulong one, const nmod_redc_ctx_t ctx, int v)
{
    nmod_t mod = ctx->mod;
    ulong n = mod.n;
    ulong c, P, two, v0, v1, x, y, w, k, half;
    slong j;
    int bit;

    FLINT_ASSERT(n % 4 == 1 && n > 5);

    c = nmod_sub(a, 4, mod);

    for (j = 0; ; j++)
    {
        if (n_jacobi_unsigned(c, n) == -1)
            break;

        /* each candidate fails with probability about 1/2 */
        if (j >= 32)
            return 0;

        if (j == 0)
            c = nmod_sub(a, 1, mod);
        else
            c = nmod_add(nmod_add(nmod_add(c, c, mod), nmod_add(c, c, mod), mod), 3, mod);
    }

    /* P = a t^2 - 2 */
    if (j == 0)
        P = nmod_add(c, 2, mod);
    else
        P = nmod_add(nmod_add(nmod_add(c, c, mod), nmod_add(c, c, mod), mod), 2, mod);

    P = _sqrt_set_nmod(P, ctx, v);
    two = _sqrt_add(one, one, ctx, v);

    /* Lucas ladder for V_k, k = (n - 1) / 4, maintaining (V_m, V_{m+1}) */
    k = n >> 2;
    v0 = two;
    v1 = P;

    for (bit = FLINT_BITS - 1 - flint_clz(k); bit >= 0; bit--)
    {
        ulong b = k >> bit;

        /* bit set: (V_m, V_{m+1}) -> (V_{2m+1}, V_{2m+2}),
           otherwise -> (V_{2m}, V_{2m+1}) */
        x = _sqrt_sub(_sqrt_mul(v0, v1, ctx, v), P, ctx, v);
        w = _sqrt_select(b, v1, v0);
        y = _sqrt_sub(_sqrt_mul(w, w, ctx, v), two, ctx, v);
        v0 = _sqrt_select(b, x, y);
        v1 = _sqrt_select(b, y, x);
    }

    /* divide by t = 2^j */
    v0 = _sqrt_normalise(v0, ctx, v);
    half = (n >> 1) + 1;   /* 1/2 mod n */
    for ( ; j > 0; j--)
        v0 = (v0 >> 1) + ((v0 & 1) ? half : 0);

    *res = v0;
    return 1;
}

/*
    Tonelli-Shanks for n - 1 = q 2^s with q odd, s >= 3. With t = a^q and
    g = k^q for a nonresidue k, g generates the 2-Sylow subgroup (of order
    2^s) and we find x with t g^x = 1 (x is even if a is a square), giving
    the root r = a^((q+1)/2) g^(x/2).

    Here a and one are in Montgomery form.
*/
FLINT_FORCE_INLINE int
_sqrt_tonelli_shanks(ulong * res, ulong a, ulong one, slong s,
                                    const nmod_redc_ctx_t ctx, int v)
{
    ulong n = ctx->mod.n;
    ulong e, r, t, u, g, b, k;
    slong m, i, j;

    k = _sqrt_nonresidue_1mod8(n);
    if (k == 0)
        return 0;

    k = _sqrt_set_nmod(k, ctx, v);

    /* u = a^((q-1)/2), g = k^((q-1)/2) */
    e = (n - 1) >> (s + 1);
    if (e == 0)
    {
        u = one;
        g = one;
    }
    else
    {
        _sqrt_pow2(&u, &g, a, k, e, one, ctx, v);
    }

    r = _sqrt_mul(u, a, ctx, v);                    /* a^((q+1)/2) */
    t = _sqrt_mul(u, r, ctx, v);                    /* a^q */
    g = _sqrt_mul(_sqrt_mul(g, g, ctx, v), k, ctx, v);  /* k^q */

    if (s < 12)
    {
        /*
            Classical descent: each step finds the order 2^i of t by
            repeated squaring. This is a chain of O(s^2) dependent
            multiplications, which is cheap for small s.
        */
        m = s;

        while (_sqrt_normalise(t, ctx, v) != one)
        {
            /* least i with 0 < i < m and t^(2^i) = 1 */
            u = t;
            for (i = 1; i < m; i++)
            {
                u = _sqrt_mul(u, u, ctx, v);
                if (_sqrt_normalise(u, ctx, v) == one)
                    break;
            }

            /* a is not a square (or n is not prime); this also
               guarantees termination since m decreases */
            if (i >= m)
                return 0;

            /* b = g^(2^(m - i - 1)) */
            b = g;
            for (j = 0; j < m - i - 1; j++)
                b = _sqrt_mul(b, b, ctx, v);

            m = i;
            g = _sqrt_mul(b, b, ctx, v);
            t = _sqrt_mul(t, g, ctx, v);
            r = _sqrt_mul(r, b, ctx, v);
        }
    }
    else
    {
        /*
            Find the bits x_j of x from the bottom, maintaining
            D[i] = (t g^X)^(2^i) where X holds the bits found so far.
            Bit j is set iff D[s-1-j] != 1, and setting it multiplies each
            D[i] by G[i+j] = g^(2^(i+j)). These O(s^2) multiplications are
            mutually independent, so unlike the classical descent this is
            bounded by multiplication throughput rather than latency.
        */
        ulong D[FLINT_BITS], G[FLINT_BITS];

        D[0] = t;
        G[0] = g;
        for (i = 1; i < s; i++)
        {
            D[i] = _sqrt_mul(D[i - 1], D[i - 1], ctx, v);
            G[i] = _sqrt_mul(G[i - 1], G[i - 1], ctx, v);
        }

        /* x_0 = 1 means a is not a square (or n is not prime) */
        if (_sqrt_normalise(D[s - 1], ctx, v) != one)
            return 0;

        for (j = 1; j < s; j++)
        {
            if (_sqrt_normalise(D[s - 1 - j], ctx, v) != one)
            {
                /* only D[i] for i < s - 1 - j are needed later */
                for (i = 0; i < s - 1 - j; i++)
                    D[i] = _sqrt_mul(D[i], G[i + j], ctx, v);

                /* r *= g^(2^j / 2) */
                r = _sqrt_mul(r, G[j - 1], ctx, v);
            }
        }
    }

    *res = r;
    return 1;
}

/*
    Computes a candidate square root of a not in {0, 1} (given both in
    standard form, a, and in Montgomery form, am) modulo the odd number n
    > 2. Returns 0 if the algorithm detects that a is not a square or that
    n is not prime; otherwise sets *res to the candidate root (in
    Montgomery form, possibly unnormalised) which the caller must verify.
*/
FLINT_FORCE_INLINE int
_nmod_sqrt_redc_candidate(ulong * res, ulong a, ulong am, ulong one,
    const nmod_redc_ctx_t ctx, int v)
{
    ulong n = ctx->mod.n;
    ulong e, t, u, b;
    slong s;

    if ((n & 3) == 3)
    {
        /* r = a^((n+1)/4) */
        *res = _sqrt_pow(am, (n >> 2) + 1, one, ctx, v);
        return 1;
    }

    if ((n & 7) == 5)
    {
        /*
            Atkin: t = 2a, b = t^((n-5)/8), y = t b^2 (y = +/- sqrt(-1) if
            a is a square), r = a b (y - 1).
        */
        t = _sqrt_add(am, am, ctx, v);
        e = n >> 3;

        if (e == 0)
            b = one;
        else
            b = _sqrt_pow(t, e, one, ctx, v);

        u = _sqrt_mul(_sqrt_mul(b, b, ctx, v), t, ctx, v);
        u = _sqrt_sub(u, one, ctx, v);
        *res = _sqrt_mul(_sqrt_mul(am, b, ctx, v), u, ctx, v);
        return 1;
    }

    s = flint_ctz(n - 1);

    /*
        The cost of Tonelli-Shanks grows with s, while Mueller's algorithm
        has a fixed overhead (mainly, Jacobi symbols) but is otherwise
        independent of s. The latter wins when n - 1 is divisible by a
        large power of two; the crossover is around s = 19 for 24-bit n
        and s = 30 for 64-bit n.
    */
    if (6 * s >= (FLINT_BITS - (slong) flint_clz(n)) + 120 &&
            _sqrt_mueller(res, a, one, ctx, v))
        return 1;

    return _sqrt_tonelli_shanks(res, am, one, s, ctx, v);
}

FLINT_FORCE_INLINE int
_nmod_sqrt_redc(ulong * res, ulong a, const nmod_redc_ctx_t ctx, int v)
{
    ulong am, one, r, n = ctx->mod.n;

    am = _sqrt_set_nmod(a, ctx, v);
    one = _sqrt_set_nmod(1, ctx, v);

    if (!_nmod_sqrt_redc_candidate(&r, a, am, one, ctx, v) ||
        _sqrt_normalise(_sqrt_mul(r, r, ctx, v), ctx, v) != am)
    {
        *res = 0;
        return 0;
    }

    r = _sqrt_get_nmod(r, ctx, v);

    /* canonical choice: the root in [0, n/2] */
    if (r > n / 2)
        r = n - r;

    *res = r;
    return 1;
}

FLINT_STATIC_NOINLINE int
_nmod_sqrt_half(ulong * res, ulong a, const nmod_redc_ctx_t ctx)
{
    return _nmod_sqrt_redc(res, a, ctx, SQRT_HALF);
}

FLINT_STATIC_NOINLINE int
_nmod_sqrt_fast(ulong * res, ulong a, const nmod_redc_ctx_t ctx)
{
    return _nmod_sqrt_redc(res, a, ctx, SQRT_FAST);
}

FLINT_STATIC_NOINLINE int
_nmod_sqrt_full(ulong * res, ulong a, const nmod_redc_ctx_t ctx)
{
    return _nmod_sqrt_redc(res, a, ctx, SQRT_FULL);
}

/*
    Exhaustive search over t <= (n-1)/2 with t^2 computed incrementally.
    Every root r has r <= n/2 or n - r <= n/2, so this is exact for any
    odd n, and it finds the canonical root. For tiny n it is faster than
    the setup of the general algorithm.
*/
int
_nmod_sqrt_small(ulong * res, ulong a, ulong n)
{
    ulong t, t2;

    FLINT_ASSERT(n & 1);
    FLINT_ASSERT(a < n);

    if (a <= 1)
    {
        *res = a;
        return 1;
    }

    t2 = 0;
    for (t = 1; t <= n / 2; t++)
    {
        /* t^2 = (t-1)^2 + 2t - 1, where 2t - 1 < n */
        t2 += 2 * t - 1;
        t2 -= (t2 >= n) ? n : 0;

        if (t2 == a)
        {
            *res = t;
            return 1;
        }
    }

    *res = 0;
    return 0;
}

int
nmod_sqrt(ulong * res, ulong a, nmod_t mod)
{
    nmod_redc_ctx_t ctx;

    FLINT_ASSERT(a < mod.n);

    if (a <= 1)
    {
        *res = a;
        return 1;
    }

    /* n = 2 is covered above */
    if ((mod.n & 1) == 0)
    {
        *res = 0;
        return 0;
    }

    if (mod.n < NMOD_SQRT_SMALL_CUTOFF)
        return _nmod_sqrt_small(res, a, mod.n);

    if (mod.norm >= FLINT_BITS / 2 + 2)
    {
        nmod_redc_half_ctx_init_nmod(ctx, mod);
        return _nmod_sqrt_half(res, a, ctx);
    }

    nmod_redc_ctx_init_nmod(ctx, mod);

    if (nmod_redc_can_use_fast(ctx))
        return _nmod_sqrt_fast(res, a, ctx);
    else
        return _nmod_sqrt_full(res, a, ctx);
}
