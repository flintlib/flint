/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "mpn_extras.h"

/*
    GCD of many integers (the content of a vector).

    The running gcd g = 2^k h (h odd) only ever decreases, and after a few
    elements it is usually the final result. So instead of chained gcds,
    each element x is first tested for divisibility by g, which costs a
    multiplication for single limbs (x h^{-1} mod B <= (B - 1) / h, and the
    low k bits of x zero) and one Hensel remainder pass otherwise; these
    tests are independent of each other, so that they run at the throughput
    of the multiplier rather than at the latency of a gcd. Only when the
    test fails is a gcd computed, from a remainder that is already reduced
    modulo h; this happens at most once per prime factor of g (counted with
    multiplicity), plus once for each element where g becomes 1.
*/

/* gcd(t, h) for h odd and t nonzero */
FLINT_FORCE_INLINE ulong
_gcd_odd(ulong t, ulong h)
{
    t >>= flint_ctz(t);
    if (t == 1 || h == 1)
        return 1;
    return mpn_gcd_11(t, h);
}

/* Elements tested one at a time before the blocked loop (which needs a
   division to set up), and the block size of the blocked loop. The
   elements are taken alternately from both ends of the vector, as for
   polynomials the leading and trailing coefficients are the most likely
   to differ from the rest in their factors. */
#define GCD_VEC_1_WARMUP 4
#define GCD_VEC_1_BLOCK 32

/* test x against g = 2^k h and replace g by gcd(g, x) if needed */
#define GCD_VEC_1_STEP(xj)                                      \
    do {                                                        \
        if (FLINT_UNLIKELY(((xj) & mask) != 0))                 \
        {                                                       \
            k = flint_ctz(xj);                                  \
            mask = (UWORD(1) << k) - 1;                         \
        }                                                       \
        /* with q = x h^{-1} mod B, q h = x + B t where t < h,  \
           so that x = -B t mod h and gcd(x, h) = gcd(t, h),    \
           and h divides x iff t = 0 */                         \
        t = n_mulhi((xj) * hinv, h);                            \
        if (FLINT_UNLIKELY(t != 0))                             \
        {                                                       \
            h = _gcd_odd(t, h);                                 \
            hinv = n_binvert(h);                                \
        }                                                       \
        if (FLINT_UNLIKELY(h == 1 && k == 0))                   \
            return 1;                                           \
    } while (0)

ulong
flint_mpn_gcd_vec_1(mp_srcptr x, slong len, ulong g)
{
    slong lo = 0, hi = len, j, s, e, w;
    ulong h, hinv, lim, mask, bad, xj, t;
    unsigned int k;
    int side;

    if (g == 0)
    {
        while (lo < hi && x[lo] == 0)
            lo++;
        if (lo == hi)
            return 0;
        g = x[lo++];

        /* the first two elements by a direct gcd (the setup below costs
           as much as a gcd, which often gives 1) */
        if (lo < hi && g != 1)
        {
            xj = x[--hi];
            if (xj != 0)
                g = mpn_gcd_1(&xj, 1, g);
        }
    }

    if (g == 1 || lo == hi)
        return g;

    k = flint_ctz(g);
    h = g >> k;
    mask = (UWORD(1) << k) - 1;
    hinv = n_binvert(h);

    for (w = 0; w < GCD_VEC_1_WARMUP && lo < hi; w++)
    {
        xj = (w & 1) ? x[--hi] : x[lo++];
        GCD_VEC_1_STEP(xj);
    }

    if (lo == hi)
        return h << k;

    lim = UWORD_MAX / h;

    for (side = 0; lo < hi; side ^= 1)
    {
        if (side == 0)
        {
            s = lo;
            e = lo = FLINT_MIN(lo + GCD_VEC_1_BLOCK, hi);
        }
        else
        {
            e = hi;
            s = hi = FLINT_MAX(hi - GCD_VEC_1_BLOCK, lo);
        }

        /* x[j] is divisible by g iff its low k bits are zero and
           x[j] h^{-1} mod B <= (B - 1) / h */
        bad = 0;
        for (j = s; j < e; j++)
            bad |= (x[j] & mask) | (ulong) (x[j] * hinv > lim);

        if (FLINT_UNLIKELY(bad != 0))
        {
            for (j = s; j < e; j++)
            {
                xj = x[j];
                GCD_VEC_1_STEP(xj);
            }

            lim = UWORD_MAX / h;
        }
    }

    return h << k;
}

/* Entries of 3 to GCD_VEC_FOLD_MAX limbs are reduced modulo h by folding
   x = x_0 + x_1 B + sum_{i >= 2} x_i (B^i mod h) into three limbs with
   precomputed powers, which unlike the Hensel remainder (a chain of
   dependent multiplications) runs at the throughput of the multiplier;
   longer entries by mpn_mod_1, whose setup cost is then negligible. */
#define GCD_VEC_FOLD_MAX 64

/* gcd(t, h) = gcd(x, h) and t = 0 iff h divides x, for x = {x1, x0}, h odd */
FLINT_FORCE_INLINE ulong
_rem_odd_2(ulong x1, ulong x0, ulong h, ulong hinv)
{
    ulong q, c;

    /* with q = x0 h^{-1} mod B and c = (q h - x0) / B < h,
       x = q h + (x1 - c) B */
    q = x0 * hinv;
    c = n_mulhi(q, h);
    if (x1 < c)
        return c - x1;
    q = (x1 - c) * hinv;
    return n_mulhi(q, h);
}

/* the same for x = {x2, x1, x0} with x2 < B - 1 */
FLINT_FORCE_INLINE ulong
_rem_odd_3(ulong x2, ulong x1, ulong x0, ulong h, ulong hinv)
{
    ulong q, c;

    /* x = q h + (x2 B + x1 - c) B */
    q = x0 * hinv;
    c = n_mulhi(q, h);
    if (x1 < c)
    {
        if (x2 == 0)
            return c - x1;
        x2--;
    }
    return _rem_odd_2(x2, x1 - c, h, hinv);
}

/* powers B^i mod h for 1 <= i < nr, for the fold; kept per thread for
   the last h, as callers (_fmpz_vec_content) make many calls with the
   same h */
typedef struct
{
    ulong r[GCD_VEC_FOLD_MAX];
    slong nr;
    ulong h;
    ulong ninv;
}
_gcd_vec_pow_t;

#if FLINT_USES_TLS
static FLINT_TLS_PREFIX _gcd_vec_pow_t _gcd_vec_pow_cache = { {0}, 0, 0, 0 };
#endif

/* a fold needs the powers up to n, costing a modular multiplication for
   each new one; the Hensel remainder is used instead when fewer than this
   many values are left */
#define GCD_VEC_FOLD_MIN_LEFT 8

static void
_gcd_vec_pow_extend(_gcd_vec_pow_t * P, slong n, ulong h)
{
    slong i;

    if (P->nr == 0 || P->h != h)
    {
        P->h = h;
        P->ninv = n_preinvert_limb(h);
        P->r[1] = n_ll_mod_preinv(1, 0, h, P->ninv);
        P->nr = 2;
    }

    for (i = P->nr; i < n; i++)
        P->r[i] = n_mulmod2_preinv(P->r[i - 1], P->r[1], h, P->ninv);

    P->nr = FLINT_MAX(P->nr, n);
}

/* {x, n} reduced to three limbs modulo h, for 3 <= n <= GCD_VEC_FOLD_MAX */
FLINT_FORCE_INLINE ulong
_rem_odd_fold(mp_srcptr x, mp_size_t n, const ulong * r, ulong h, ulong hinv)
{
    ulong a2, a1, a0, b2, b1, b0, p1, p0;
    mp_size_t i;

    /* two accumulators for the even and odd terms, so that the carry
       chains overlap; each sum stays below n B^2 */
    a2 = 0;
    a1 = x[1];
    a0 = x[0];
    b2 = b1 = b0 = 0;
    for (i = 2; i + 1 < n; i += 2)
    {
        umul_ppmm(p1, p0, x[i], r[i]);
        add_sssaaaaaa(a2, a1, a0, a2, a1, a0, 0, p1, p0);
        umul_ppmm(p1, p0, x[i + 1], r[i + 1]);
        add_sssaaaaaa(b2, b1, b0, b2, b1, b0, 0, p1, p0);
    }
    if (i < n)
    {
        umul_ppmm(p1, p0, x[i], r[i]);
        add_sssaaaaaa(a2, a1, a0, a2, a1, a0, 0, p1, p0);
    }
    add_sssaaaaaa(a2, a1, a0, a2, a1, a0, b2, b1, b0);

    return _rem_odd_3(a2, a1, a0, h, hinv);
}

/* the same with x_1 also folded, into two limbs, for n h <= B */
FLINT_FORCE_INLINE ulong
_rem_odd_fold2(mp_srcptr x, mp_size_t n, const ulong * r, ulong h, ulong hinv)
{
    ulong a1, a0, b1, b0, p1, p0;
    mp_size_t i;

    a1 = 0;
    a0 = x[0];
    b1 = b0 = 0;
    for (i = 1; i + 1 < n; i += 2)
    {
        umul_ppmm(p1, p0, x[i], r[i]);
        add_ssaaaa(a1, a0, a1, a0, p1, p0);
        umul_ppmm(p1, p0, x[i + 1], r[i + 1]);
        add_ssaaaa(b1, b0, b1, b0, p1, p0);
    }
    if (i < n)
    {
        umul_ppmm(p1, p0, x[i], r[i]);
        add_ssaaaa(a1, a0, a1, a0, p1, p0);
    }
    add_ssaaaa(a1, a0, a1, a0, b1, b0);

    return _rem_odd_2(a1, a0, h, hinv);
}

/* for h below this, n h < B for all n <= GCD_VEC_FOLD_MAX = 64 */
#define GCD_VEC_FOLD2_HMAX (UWORD(1) << (FLINT_BITS - 7))

/* gcd of 2^k h (h odd, 2^k h > 1) and the nonzero values {xp[i], xn[i]} */
static ulong
_gcd_vec_ptr_1(const mp_srcptr * xp, const mp_size_t * xn, slong len, ulong h, unsigned int k)
{
    _gcd_vec_pow_t Pl, * P = &Pl;
    slong i, nr;
    ulong hinv, mask, t, x0;
    mp_size_t n;
    mp_srcptr x;
    int extended = 0;

    mask = (UWORD(1) << k) - 1;
    hinv = n_binvert(h);

    /* the powers are used from a local copy, which is faster to access */
    Pl.nr = 0;
    Pl.h = 0;
#if FLINT_USES_TLS
    if (_gcd_vec_pow_cache.h == h && _gcd_vec_pow_cache.nr != 0)
    {
        Pl.h = h;
        Pl.ninv = _gcd_vec_pow_cache.ninv;
        Pl.nr = _gcd_vec_pow_cache.nr;
        for (i = 1; i < Pl.nr; i++)
            Pl.r[i] = _gcd_vec_pow_cache.r[i];
    }
#endif
    nr = Pl.nr;

    for (i = 0; i < len; i++)
    {
        n = xn[i];
        if (n == 0)
            continue;

        x = xp[i];
        x0 = x[0];

        /* if x0 = 0, x has at least FLINT_BITS > k trailing zeros */
        if (FLINT_UNLIKELY((x0 & mask) != 0))
        {
            k = flint_ctz(x0);
            mask = (UWORD(1) << k) - 1;
        }

        /* h is odd: zero low limbs can be dropped */
        if (FLINT_UNLIKELY(x0 == 0))
        {
            do
            {
                x++;
                n--;
            }
            while (x[0] == 0);
            x0 = x[0];
        }

        if (h != 1)
        {
            if (n == 1)
                t = n_mulhi(x0 * hinv, h);
            else if (n == 2)
                t = _rem_odd_2(x[1], x0, h, hinv);
            else if (n <= GCD_VEC_FOLD_MAX && (n <= nr || len - i >= GCD_VEC_FOLD_MIN_LEFT))
            {
                if (FLINT_UNLIKELY(n > nr))
                {
                    _gcd_vec_pow_extend(P, n, h);
                    nr = P->nr;
                    extended = 1;
                }
                if (h < GCD_VEC_FOLD2_HMAX)
                    t = _rem_odd_fold2(x, n, P->r, h, hinv);
                else
                    t = _rem_odd_fold(x, n, P->r, h, hinv);
            }
#if defined(mpn_modexact_1_odd)
            else if (n <= GCD_VEC_FOLD_MAX)
                t = mpn_modexact_1_odd(x, n, h);
#endif
            else
                t = mpn_mod_1(x, n, h);

            if (FLINT_UNLIKELY(t != 0))
            {
                h = _gcd_odd(t, h);
                hinv = n_binvert(h);
                nr = 0;
            }
        }

        if (FLINT_UNLIKELY(h == 1 && k == 0))
            return 1;
    }

#if FLINT_USES_TLS
    if (extended && Pl.h == h)
    {
        _gcd_vec_pow_cache.h = h;
        _gcd_vec_pow_cache.ninv = Pl.ninv;
        _gcd_vec_pow_cache.nr = Pl.nr;
        for (i = 1; i < Pl.nr; i++)
            _gcd_vec_pow_cache.r[i] = Pl.r[i];
    }
#else
    (void) extended;
#endif

    return h << k;
}

/* number of trailing zero bits of the nonzero {x, n} */
FLINT_FORCE_INLINE flint_bitcnt_t
_mpn_ctz(mp_srcptr x)
{
    mp_size_t i = 0;
    while (x[i] == 0)
        i++;
    return (flint_bitcnt_t) i * FLINT_BITS + flint_ctz(x[i]);
}

/* {g, return value} = 2^v h; g has room for the result */
static mp_size_t
_set_2exp_1(mp_ptr g, ulong h, flint_bitcnt_t v)
{
    mp_size_t s = v / FLINT_BITS;
    unsigned int c = v % FLINT_BITS;

    flint_mpn_zero(g, s);
    g[s] = h << c;
    if (c != 0 && (h >> (FLINT_BITS - c)) != 0)
    {
        g[s + 1] = h >> (FLINT_BITS - c);
        return s + 2;
    }
    return s + 1;
}

/* The gcd 2^v h of the nonzero values {xp[i], xn[i]} with 2^v h, h odd,
   written to g which has room for it */
static mp_size_t
_gcd_vec_finish_1(mp_ptr g, ulong h, flint_bitcnt_t v,
    const mp_srcptr * xp, const mp_size_t * xn, slong len)
{
    slong i;

    if (v < FLINT_BITS && (h << v) >> v == h)
    {
        if (h != 1 || v != 0)
            h = _gcd_vec_ptr_1(xp, xn, len, h, v);
        g[0] = h;
        return 1;
    }

    /* the power of two does not fit with h in a limb: the odd part and
       the power of two separately */
    if (h != 1)
        h = _gcd_vec_ptr_1(xp, xn, len, h, 0);
    for (i = 0; i < len && v != 0; i++)
        if (xn[i] != 0)
            v = FLINT_MIN(v, _mpn_ctz(xp[i]));

    return _set_2exp_1(g, h, v);
}

/* Divisors of at least this many limbs are tested by computing the
   remainder (which is needed anyway if the test fails) rather than with
   flint_mpn_divisible, which is only faster for short divisors. */
#define GCD_VEC_DIVISIBLE_CUTOFF 100

/* state of the gcd with a multi-limb running gcd 2^v g, g odd */
typedef struct
{
    mp_ptr g;
    mp_size_t gn;
    flint_bitcnt_t v;
    int fresh;      /* g was set from a single value */
    mp_ptr tmp;     /* 3 max(gn, xn) + 3 limbs */
    ulong ginv;     /* g^{-1} mod B, or 0 if not computed */
}
_gcd_vec_state_t;

/* whether the odd {g1, g0} (g1 != 0) divides {x, n} for n = 2 or 3, with
   ginv = g0^{-1} mod B: the quotient q has n - 1 limbs and is the Hensel
   quotient x g^{-1} mod B^{n-1}, so that g divides x iff q g = x */
FLINT_FORCE_INLINE int
_divisible_2(mp_srcptr x, mp_size_t n, ulong g1, ulong g0, ulong ginv)
{
    ulong q, p1, p0, r1, r0, s, y1, y0;

    if (n == 2)
    {
        y1 = x[1];
        y0 = x[0];
    }
    else
    {
        /* {y1, y0} = (x - q0 g) / B, which must be nonnegative */
        q = x[0] * ginv;
        umul_ppmm(p1, p0, q, g0);
        umul_ppmm(r1, r0, q, g1);
        add_ssaaaa(r1, r0, r1, r0, 0, p1);
        if (x[2] < r1 || (x[2] == r1 && x[1] < r0))
            return 0;
        sub_ddmmss(y1, y0, x[2], x[1], r1, r0);
    }

    /* {y1, y0} = q g with q < B */
    q = y0 * ginv;
    umul_ppmm(p1, p0, q, g0);
    umul_ppmm(r1, r0, q, g1);
    (void) p0;
    s = r0 + p1;
    return r1 == 0 && s >= r0 && s == y1;
}

/* replace g by gcd(g, x) for x nonzero */
static void
_gcd_vec_step(_gcd_vec_state_t * S, mp_srcptr x, mp_size_t n)
{
    mp_ptr g = S->g, r, q;
    mp_size_t gn = S->gn, rn;

    if (S->v != 0)
        S->v = FLINT_MIN(S->v, _mpn_ctz(x));

    /* g is odd, so that dividing x by B changes neither whether g divides
       it nor gcd(g, x); this matters for entries with many factors of 2
       (factorials), which would otherwise give long quotients */
    while (x[0] == 0)
    {
        x++;
        n--;
    }

    r = S->tmp;

    if (n >= gn && !S->fresh)
    {
        /* g usually divides x once it has been through a gcd */
        if (gn == 2 && n <= 3)
        {
            if (S->ginv == 0)
                S->ginv = n_binvert(g[0]);
            if (_divisible_2(x, n, g[1], g[0], S->ginv))
                return;
        }
        else if (gn <= FLINT_MPN_DIVEXACT_SMALL_BN)
        {
            if (_flint_mpn_divisible_small(x, n, g, gn, 0))
                return;
        }
        else if (gn < GCD_VEC_DIVISIBLE_CUTOFF)
        {
            if (flint_mpn_divisible(x, n, g, gn))
                return;
        }
        else
        {
            /* the remainder, which the gcd needs if the test fails */
            q = S->tmp + gn;
            flint_mpn_tdiv_qr(q, r, x, n, g, gn);
            rn = gn;
            MPN_NORM(r, rn);
            if (rn == 0)
                return;
            x = r;
            n = rn;
            r = S->tmp + gn;
        }
    }

    /* g = gcd(g, x), g being odd (flint_mpn_gcd reduces the longer
       operand first) */
    if (x != r)
        flint_mpn_copyi(r, x, n);
    q = r + FLINT_MAX(n, gn) + 1;
    if (n >= gn)
        S->gn = _flint_mpn_gcd_inplace(q, r, n, g, gn);
    else
        S->gn = _flint_mpn_gcd_inplace(q, g, gn, r, n);
    flint_mpn_copyi(g, q, S->gn);
    S->fresh = 0;
    S->ginv = 0;
}

/* sets S->g to the odd part of {x, n} and returns the number of trailing
   zeros */
static flint_bitcnt_t
_odd_part(mp_ptr g, mp_size_t * gn, mp_srcptr x, mp_size_t n)
{
    flint_bitcnt_t v = _mpn_ctz(x);
    mp_size_t m = n - v / FLINT_BITS;

    if (v % FLINT_BITS != 0)
        mpn_rshift(g, x + v / FLINT_BITS, m, v % FLINT_BITS);
    else if (g != x + v / FLINT_BITS)
        flint_mpn_copyi(g, x + v / FLINT_BITS, m);
    MPN_NORM(g, m);
    *gn = m;
    return v;
}

mp_size_t
flint_mpn_gcd_vec(mp_ptr g, mp_size_t gn, const mp_srcptr * xp, const mp_size_t * xn, slong len)
{
    slong i, j, n1, imin, imin2;
    slong nmin, nmin2;      /* not mp_size_t, which may be 32 bits */
    mp_size_t n, xmax;
    ulong g1;

    /* a single-limb g: the values in place */
    if (gn == 1)
    {
        unsigned int k = flint_ctz(g[0]);
        if (g[0] != 1)
            g[0] = _gcd_vec_ptr_1(xp, xn, len, g[0] >> k, k);
        return 1;
    }

    /* number of single-limb values, the two shortest longer ones and the
       longest one */
    n1 = 0;
    nmin = nmin2 = WORD_MAX;
    imin = imin2 = -1;
    xmax = 0;
    for (i = 0; i < len; i++)
    {
        n = xn[i];
        if (n == 1)
            n1++;
        else if (n != 0)
        {
            if (n < nmin2)
            {
                if (n < nmin)
                {
                    nmin2 = nmin;
                    imin2 = imin;
                    nmin = n;
                    imin = i;
                }
                else
                {
                    nmin2 = n;
                    imin2 = i;
                }
            }
            xmax = FLINT_MAX(xmax, n);
        }
    }

    if (gn == 1 || n1 != 0)
    {
        /* the gcd has a single limb: first the gcd of the single-limb
           values, gathered in blocks, then the others */
        mp_limb_t buf[GCD_VEC_1_BLOCK];
        unsigned int k;

        g1 = (gn == 1) ? g[0] : 0;

        for (i = 0; i < len && g1 != 1; )
        {
            for (j = 0; j < GCD_VEC_1_BLOCK && i < len; i++)
                if (xn[i] == 1)
                    buf[j++] = xp[i][0];
            g1 = flint_mpn_gcd_vec_1(buf, j, g1);
        }

        if (g1 != 1 && gn > 1)
        {
            mp_srcptr gp = g;
            k = flint_ctz(g1);
            g1 = _gcd_vec_ptr_1(&gp, &gn, 1, g1 >> k, k);
        }

        if (g1 != 1 && imin != -1)
        {
            k = flint_ctz(g1);
            g1 = _gcd_vec_ptr_1(xp, xn, len, g1 >> k, k);
        }

        g[0] = g1;
        return 1;
    }

    if (gn == 0 && imin == -1)
        return 0;

    if (imin == -1)
        return gn;

    {
        _gcd_vec_state_t S[1];
        mp_ptr gc = NULL;
        mp_size_t gcn = gn;
        TMP_INIT;
        TMP_START;

        if (gn != 0 && gn <= nmin)
        {
            /* start from the given g */
            S->v = _odd_part(g, &S->gn, g, gn);
            S->fresh = 0;
            imin = imin2 = -1;
        }
        else
        {
            /* start from the shortest value, then the given g if any */
            if (gn != 0)
            {
                gc = TMP_ALLOC(gn * sizeof(mp_limb_t));
                flint_mpn_copyi(gc, g, gn);
                xmax = FLINT_MAX(xmax, gn);
            }
            S->v = _odd_part(g, &S->gn, xp[imin], xn[imin]);
            S->fresh = 1;
        }

        S->g = g;
        S->ginv = 0;
        S->tmp = TMP_ALLOC((3 * FLINT_MAX(S->gn, xmax) + 3) * sizeof(mp_limb_t));

        if (gc != NULL && S->gn > 1)
        {
            _gcd_vec_step(S, gc, gcn);
            gc = NULL;
        }

        /* then the second shortest value, for a short g early */
        if (imin2 != -1 && S->gn > 1)
            _gcd_vec_step(S, xp[imin2], xn[imin2]);
        else
            imin2 = -1;

        for (i = 0; i < len; i++)
        {
            if (S->gn == 1)
            {
                /* finish with single limbs (and the given g if it has
                   not been used) */
                gn = _gcd_vec_finish_1(g, g[0], S->v, xp + i, xn + i, len - i);
                if (gc != NULL)
                    gn = flint_mpn_gcd_vec(g, gn, (const mp_srcptr *) &gc, &gcn, 1);
                goto done;
            }

            if (xn[i] != 0 && i != imin && i != imin2)
                _gcd_vec_step(S, xp[i], xn[i]);
        }

        FLINT_ASSERT(gc == NULL);

        /* g times 2^v */
        gn = S->gn;
        if (S->v != 0)
        {
            mp_size_t s = S->v / FLINT_BITS;
            mp_limb_t cy;

            if (s != 0)
            {
                flint_mpn_copyd(g + s, g, gn);
                flint_mpn_zero(g, s);
            }

            /* the result is at most the original g, so that the carry limb
               (if any) fits */
            if (S->v % FLINT_BITS != 0)
            {
                cy = mpn_lshift(g + s, g + s, gn, S->v % FLINT_BITS);
                if (cy != 0)
                    g[gn++ + s] = cy;
            }

            gn += s;
        }

done:
        TMP_END;
    }

    return gn;
}
